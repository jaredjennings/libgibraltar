/* Author:  Matthew Curry
 * Email:   mlcurry@sandia.gov
 *
 * Defines the internal CPU interface for Gibraltar.
 */
#include "gibraltar.h"
#ifdef __cplusplus
extern "C" {
#endif

int gib_cpu_init ( int n, int m, gib_context *c );
int gib_cpu_destroy ( gib_context c );
int gib_cpu_alloc ( void **buffers, int buf_size, int *ld, gib_context c );
int gib_cpu_free ( void *buffers );
int gib_cpu_generate ( void *buffers, int buf_size, gib_context c );
int gib_cpu_generate_nc ( void *buffers, int buf_size, int work_size, 
		gib_context c);
int gib_cpu_recover_sparse ( void *buffers, int buf_size, char *failed_bufs, 
		gib_context c );
int gib_cpu_recover_sparse_nc ( void *buffers, int buf_size, int work_size, 
		char *failed_bufs, gib_context c );
int gib_cpu_recover ( void *buffers, int buf_size, int *buf_ids, int recover_last,
		gib_context c );
int gib_cpu_recover_nc ( void *buffers, int buf_size, int work_size, int *buf_ids, int recover_last,
		gib_context c );

#ifdef __cplusplus
}
#endif

