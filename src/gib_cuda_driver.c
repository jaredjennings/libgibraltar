/* Author:  Matthew Curry
 * Email:   mlcurry@sandia.gov
 *
 * This file defines the Gibraltar implementation functions which use the CUDA 
 * driver API.  
 */

/* TODO:  
   - CPU failover.
     - Error check macro support.
   - Noncontiguous only occurs on CPU!
   - Document the environment variables used.
*/

/* If compute capability 1.3 or higher is available, this should be set. 
 * If it's set by the user at compile time, respect it.
 */
#ifndef GIB_USE_MMAP
#define GIB_USE_MMAP 1
#endif

/* Size of each GPU buffer; n+m will be allocated */
#if !GIB_USE_MMAP
int gib_buf_size = 1024*1024; 
#endif

#include "../inc/gibraltar.h"
#include "../inc/gib_galois.h"
#include "../inc/gib_cpu_funcs.h"
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <cuda.h>

int cudaInitialized = 0;

struct gpu_context_t {
  CUdevice dev;
  CUmodule module;
  CUcontext pCtx;
  CUfunction checksum;
  CUfunction recover_sparse;
  CUfunction recover;
  CUdeviceptr buffers;
};

typedef struct gpu_context_t * gpu_context;

/* This macro checks for an error in the command given.  If it fails, the entire
 * program is killed.
 * TODO:  Fail over to CPU code if an error occurs.
 */
#define ERROR_CHECK_FAIL(cmd) {						\
    CUresult rc = cmd;							\
    if (rc != CUDA_SUCCESS) {						\
      fprintf(stderr, "%s failed with %i at %i in %s\n", #cmd, rc,	\
	      __LINE__,  __FILE__);					\
      exit(EXIT_FAILURE);						\
    }									\
  }

/* Massive performance increases come from compiling the CUDA kernels 
   specifically for the coding process at hand.  This does so with the
   following command line:
      nvcc --ptx -DN=n -DM=m src/gib_cuda_checksum.cu -o gib_cuda_n+m.ptx
   This is called in a separate process fork'd from the original.  This
   function should never return, and the parent process should wait on the
   return code from the compiler before resuming operation.
*/
void gib_cuda_compile(int n, int m, char *filename) { /* never returns */
  char *executable = "nvcc";

  if (getenv("PATH") == NULL) {
    fprintf(stderr, "Your path is not set.  Please set it, and include the path"
	    " to nvcc.");
    exit(1);
  }

  char argv1[100], argv2[100];
  sprintf(argv1, "-DN=%i", n);
  sprintf(argv2, "-DM=%i", m);
  if (getenv("GIB_SRC_DIR") == NULL) {
    fprintf(stderr,
	    "Your environment is not completely set.  Please indicate a \n"
	    "directory where Gibraltar kernel sources can be found.  This\n"
	    "should not be a publicly accessible directory.\n");
    exit(1);
  }

  char src_filename[100];
  sprintf(src_filename, "%s/gib_cuda_checksum.cu", getenv("GIB_SRC_DIR"));
  char *const argv[] = {executable,
			"--ptx", 
			argv1, 
			argv2, 
			src_filename, 
			"-o", 
			filename, 
			NULL};

  execvp(argv[0], argv);
  perror("execve(nvcc)");
  fflush(0);
  exit(-1);
}

/* Initializes the CPU and GPU runtimes. */
int gib_init ( int n, int m, gib_context *c ) {
  static CUcontext pCtx;
  static CUdevice dev;
  if (m < 2 || n < 2) {
    fprintf(stderr, 
	    "It makes little sense to use Reed-Solomon coding when n or m is\n"
	    "less than two.  Use XOR or replication instead.\n");
    exit(1);
  }
  int rc_i = gib_cpu_init(n,m,c);
  if (rc_i != GIB_SUC) {
    fprintf(stderr, "gib_cpu_init returned %i\n", rc_i);
    exit(EXIT_FAILURE);
  }

  int gpu_id = 0;
  if (!cudaInitialized) {
    /* Initialize the CUDA runtime */
    int device_count;
    ERROR_CHECK_FAIL(cuInit(0));
    ERROR_CHECK_FAIL(cuDeviceGetCount(&device_count));
    if (getenv("GIB_GPU_ID") != NULL) {
        gpu_id = atoi(getenv("GIB_GPU_ID"));
        if (device_count <= gpu_id) {
            fprintf(stderr,
                    "GIB_GPU_ID is set to an invalid value (%i).  There are \n"
                    "only %i GPUs in the system.  Please specify another \n"
                    "value.\n", gpu_id, device_count);
            exit(-1);
        }
    }
    cudaInitialized = 1;
  }
  ERROR_CHECK_FAIL(cuDeviceGet(&dev, gpu_id));
#if GIB_USE_MMAP
    ERROR_CHECK_FAIL(cuCtxCreate(&pCtx, CU_CTX_MAP_HOST, dev));	
#else
    ERROR_CHECK_FAIL(cuCtxCreate(&pCtx, 0, dev));
#endif

  /* Initialize the Gibraltar context */
  gpu_context gpu_c = (gpu_context) malloc(sizeof(struct gpu_context_t));
  gpu_c->dev = dev;
  gpu_c->pCtx = pCtx;
  (*c)->acc_context = (void *)gpu_c;
	
  /* Determine whether the PTX has been generated or not by attempting to
   * open it read-only.
   */
  if (getenv("GIB_CACHE_DIR") == NULL) {
    fprintf(stderr, 
	   "Your environment is not completely set.  Please indicate a \n"
	   "directory where generated files may be placed with the \n"
	   "GIB_CACHE_DIR environment variable.  This directory should\n"
	   "not be publicly accessible and should exist.\n");
    exit(-1);
  }

  /* Try to open the appropriate ptx file.  If it doesn't exist, compile a
   * new one.
   */
  int filename_len = strlen(getenv("GIB_CACHE_DIR")) + 
    strlen("/gib_cuda_+.ptx") + log10(n)+1 + log10(m)+1 + 1;
  char *filename = (char *)malloc(filename_len);
  sprintf(filename, "%s/gib_cuda_%i+%i.ptx", getenv("GIB_CACHE_DIR"), n, m);

  FILE *fp = fopen(filename, "r");
  if (fp == NULL) {
    /* Compile the ptx and open it */
    int pid = fork();
    if (pid == -1) {
      perror("Forking for nvcc");
      exit(-1);
    }
    if (pid == 0) {
      gib_cuda_compile(n, m, filename); /* never returns */
    }
    int status;
    wait(&status);
    if (status != 0) {
      printf("Waiting for the compiler failed.\n");
      printf("The exit status was %i\n", WEXITSTATUS(status));
      printf("The child did%s exit normally.\n", 
	     (WIFEXITED(status)) ? "" : " NOT");
      
      exit(-1);
    }
    fp = fopen(filename, "r");
    if (fp == NULL) {
      perror(filename);
      exit(-1);
    }
  }
  fclose(fp);

  /* If we got here, the ptx file exists.  Use it. */
  ERROR_CHECK_FAIL(cuModuleLoad(&(gpu_c->module), filename));
  ERROR_CHECK_FAIL(cuModuleGetFunction(&(gpu_c->checksum), (gpu_c->module), 
	       "_Z14gib_checksum_dP11shmem_bytesi"));
  ERROR_CHECK_FAIL(cuModuleGetFunction(&(gpu_c->recover),
	       (gpu_c->module), 
	       "_Z13gib_recover_dP11shmem_bytesii"));
	
  /* Initialize the math libraries */
  gib_galois_init();
  unsigned char F[256*256];
  gib_galois_gen_F(F, m, n);

  /* Initialize/Allocate GPU-side structures */
  CUdeviceptr log_d, ilog_d, F_d;
  ERROR_CHECK_FAIL(cuModuleGetGlobal(&log_d, NULL, gpu_c->module, "gf_log_d"));
  ERROR_CHECK_FAIL(cuMemcpyHtoD(log_d, gib_gf_log, 256));
  ERROR_CHECK_FAIL(cuModuleGetGlobal(&ilog_d, NULL, gpu_c->module, 
				     "gf_ilog_d"));
  ERROR_CHECK_FAIL(cuMemcpyHtoD(ilog_d, gib_gf_ilog, 256));
  ERROR_CHECK_FAIL(cuModuleGetGlobal(&F_d, NULL, gpu_c->module, "F_d"));
  ERROR_CHECK_FAIL(cuMemcpyHtoD(F_d, F, m*n));
#if !GIB_USE_MMAP
  ERROR_CHECK_FAIL(cuMemAlloc(&(gpu_c->buffers), (n+m)*gib_buf_size));
#endif
  ERROR_CHECK_FAIL(cuCtxPopCurrent((&gpu_c->pCtx)));
  free(filename);
  return GIB_SUC;
}

int gib_destroy ( gib_context c ) {
  /* TODO:  Make sure everything created in gib_init is destroyed here. */
  ERROR_CHECK_FAIL(cuCtxPushCurrent(((gpu_context)(c->acc_context))->pCtx));
  int rc_i = gib_cpu_destroy(c);
  if (rc_i != GIB_SUC) {
    printf("gib_cpu_destroy returned %i\n", rc_i);
    exit(EXIT_FAILURE);
  }
  gpu_context gpu_c = (gpu_context) c->acc_context;
#if !GIB_USE_MMAP
  ERROR_CHECK_FAIL(cuMemFree(gpu_c->buffers));
#endif
  ERROR_CHECK_FAIL(cuModuleUnload(gpu_c->module));
  ERROR_CHECK_FAIL(cuCtxDestroy(gpu_c->pCtx));
  return GIB_SUC;
}

int gib_alloc ( void **buffers, int buf_size, int *ld, gib_context c ) {
  ERROR_CHECK_FAIL(cuCtxPushCurrent(((gpu_context)(c->acc_context))->pCtx));
#if GIB_USE_MMAP
  ERROR_CHECK_FAIL(cuMemHostAlloc(buffers, (c->n+c->m)*buf_size, 
				  CU_MEMHOSTALLOC_DEVICEMAP));
#else
  ERROR_CHECK_FAIL(cuMemAllocHost(buffers, (c->n+c->m)*buf_size));
#endif
  *ld = buf_size;
  ERROR_CHECK_FAIL(cuCtxPopCurrent(&((gpu_context)(c->acc_context))->pCtx));
  return GIB_SUC;
}

int gib_free ( void *buffers, gib_context c ) {
  ERROR_CHECK_FAIL(cuCtxPushCurrent(((gpu_context)(c->acc_context))->pCtx));
  ERROR_CHECK_FAIL(cuMemFreeHost(buffers));
  ERROR_CHECK_FAIL(cuCtxPopCurrent(&((gpu_context)(c->acc_context))->pCtx));
  return GIB_SUC;
}

int gib_generate ( void *buffers, int buf_size, gib_context c ) {
  ERROR_CHECK_FAIL(cuCtxPushCurrent(((gpu_context)(c->acc_context))->pCtx));
  /* Do it all at once if the buffers are small enough */
#if !GIB_USE_MMAP
  /* This is too large to do at once in the GPU memory we have allocated.
   * Split it into several noncontiguous jobs. 
   */
  if (buf_size > gib_buf_size) {
    int rc = gib_generate_nc(buffers, buf_size, buf_size, c);
    ERROR_CHECK_FAIL(cuCtxPopCurrent(&((gpu_context)(c->acc_context))->pCtx));
    return rc;
  }
#endif

  int nthreads_per_block = 128;
  int fetch_size = sizeof(int)*nthreads_per_block;
  int nblocks = (buf_size + fetch_size - 1)/fetch_size;
  gpu_context gpu_c = (gpu_context) c->acc_context;
  
  unsigned char F[256*256];
  gib_galois_gen_F(F, c->m, c->n);
  CUdeviceptr F_d;
  ERROR_CHECK_FAIL(cuModuleGetGlobal(&F_d, NULL, gpu_c->module, "F_d"));
  ERROR_CHECK_FAIL(cuMemcpyHtoD(F_d, F, (c->m)*(c->n)));
  
#if !GIB_USE_MMAP
  /* Copy the buffers to memory */
  ERROR_CHECK_FAIL(cuMemcpyHtoD(gpu_c->buffers, buffers, 
				(c->n)*buf_size));
#endif
  /* Configure and launch */
  ERROR_CHECK_FAIL(cuFuncSetBlockShape(gpu_c->checksum, nthreads_per_block,
				       1, 1));
  int offset = 0;
  void *ptr;
#if GIB_USE_MMAP
  CUdeviceptr cpu_buffers;
  ERROR_CHECK_FAIL(cuMemHostGetDevicePointer(&cpu_buffers, buffers, 0));
  ptr = (void *)cpu_buffers;
#else
  ptr = (void *)(gpu_c->buffers);
#endif
  ERROR_CHECK_FAIL(cuParamSetv(gpu_c->checksum, offset, &ptr, sizeof(ptr)));
  offset += sizeof(ptr);
  ERROR_CHECK_FAIL(cuParamSetv(gpu_c->checksum, offset, &buf_size,
			       sizeof(buf_size)));
  offset += sizeof(buf_size);
  ERROR_CHECK_FAIL(cuParamSetSize(gpu_c->checksum, offset));
  ERROR_CHECK_FAIL(cuLaunchGrid(gpu_c->checksum, nblocks, 1));

  /* Get the results back */
#if !GIB_USE_MMAP
  CUdeviceptr tmp_d = gpu_c->buffers + c->n*buf_size;
  void *tmp_h = (void *)((unsigned char *)(buffers) + c->n*buf_size);
  ERROR_CHECK_FAIL(cuMemcpyDtoH(tmp_h, tmp_d, (c->m)*buf_size));
#else
  ERROR_CHECK_FAIL(cuCtxSynchronize());
#endif
  ERROR_CHECK_FAIL(cuCtxPopCurrent(&((gpu_context)(c->acc_context))->pCtx));
  return GIB_SUC; 
}

int gib_recover ( void *buffers, int buf_size, int *buf_ids, int recover_last,
		  gib_context c ) {
  ERROR_CHECK_FAIL(cuCtxPushCurrent(((gpu_context)(c->acc_context))->pCtx));
#if !GIB_USE_MMAP
  if (buf_size > gib_buf_size) {
    int rc = gib_cpu_recover(buffers, buf_size, buf_ids, recover_last, c);
    ERROR_CHECK_FAIL(cuCtxPopCurrent(&((gpu_context)(c->acc_context))->pCtx));
    return rc;
  }
#endif

  int i, j;
  int n = c->n;
  int m = c->m;
  unsigned char A[128*128], inv[128*128], modA[128*128];
  for (i = n; i < n+recover_last; i++)
    if (buf_ids[i] >= n) {
      fprintf(stderr, "Attempting to recover a parity buffer, not allowed\n");
      return GIB_ERR;
    }

  gib_galois_gen_A(A, m+n, n);

  /* Modify the matrix to have the failed drives reflected */
  for (i = 0; i < n; i++) 
    for (j = 0; j < n; j++) 
      modA[i*n+j] = A[buf_ids[i]*n+j];

  gib_galois_gaussian_elim(modA, inv, n, n);

  /* Copy row buf_ids[i] into row i */
  for (i = n; i < n+recover_last; i++)
    for (j = 0; j < n; j++)
      modA[i*n+j] = inv[buf_ids[i]*n+j];

  int nthreads_per_block = 128;
  int fetch_size = sizeof(int)*nthreads_per_block;
  int nblocks = (buf_size + fetch_size - 1)/fetch_size;
  gpu_context gpu_c = (gpu_context) c->acc_context;

  CUdeviceptr F_d;
  ERROR_CHECK_FAIL(cuModuleGetGlobal(&F_d, NULL, gpu_c->module, "F_d"));
  ERROR_CHECK_FAIL(cuMemcpyHtoD(F_d, modA+n*n, (c->m)*(c->n)));

#if !GIB_USE_MMAP
  ERROR_CHECK_FAIL(cuMemcpyHtoD(gpu_c->buffers, buffers, (c->n)*buf_size));
#endif
  ERROR_CHECK_FAIL(cuFuncSetBlockShape(gpu_c->recover, nthreads_per_block, 
				       1, 1));
  int offset = 0;
  void *ptr;
#if GIB_USE_MMAP
  CUdeviceptr cpu_buffers;
  ERROR_CHECK_FAIL(cuMemHostGetDevicePointer(&cpu_buffers, buffers, 0));
  ptr = (void *)cpu_buffers;
#else
  ptr = (void *)gpu_c->buffers;
#endif
  ERROR_CHECK_FAIL(cuParamSetv(gpu_c->recover, offset, &ptr, sizeof(ptr)));
  offset += sizeof(ptr);
  ERROR_CHECK_FAIL(cuParamSetv(gpu_c->recover, offset, &buf_size, 
			       sizeof(buf_size)));
  offset += sizeof(buf_size);
  ERROR_CHECK_FAIL(cuParamSetv(gpu_c->recover, offset, &recover_last, 
			       sizeof(recover_last)));
  offset += sizeof(recover_last);
  ERROR_CHECK_FAIL(cuParamSetSize(gpu_c->recover, offset));
  ERROR_CHECK_FAIL(cuLaunchGrid(gpu_c->recover, nblocks, 1));
#if !GIB_USE_MMAP
  CUdeviceptr tmp_d = gpu_c->buffers + c->n*buf_size;
  void *tmp_h = (void *)((unsigned char *)(buffers) + c->n*buf_size);
  ERROR_CHECK_FAIL(cuMemcpyDtoH(tmp_h, tmp_d, recover_last*buf_size));
#else
  cuCtxSynchronize();
#endif
  ERROR_CHECK_FAIL(cuCtxPopCurrent(&((gpu_context)(c->acc_context))->pCtx));
  return GIB_SUC;
}

/* The inclusion of memory mapping has obviated the need for this before
   it was implemented.  It's done in the CPU to make it work, but there is
   no attempt to make it fast as there appears to be little need.  A GPU 
   upgrade fixes this. 
   
   TODO:  The MMapped version can benefit from this if the buffer isn't full.
   Bring this to life for that implementation only.
*/
int gib_generate_nc ( void *buffers, int buf_size, int work_size,  
		      gib_context c) {
  return gib_cpu_generate_nc(buffers, buf_size, work_size, c);
}
int gib_recover_nc ( void *buffers, int buf_size, int work_size, int *buf_ids, 
		     int recover_last, gib_context c ) {
  return gib_cpu_recover_nc(buffers, buf_size, work_size, buf_ids, 
			    recover_last, c);
}

