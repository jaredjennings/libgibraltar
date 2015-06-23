/* Author:  Matthew Curry
 * Email:   mlcurry@sandia.gov
 *
 * A CPU implementation of the Galois arithmetic operations needed for both the
 * low-performance CPU version of Gibraltar and initializing the GPU version.
 */

#include "../inc/gib_galois.h"
#include "../inc/gibraltar.h" /* For error codes */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

gib_scalar gib_gf_log[GIB_GALOIS_DEGREE];
gib_scalar gib_gf_ilog[GIB_GALOIS_DEGREE];

gib_scalar gib_galois_mul(gib_scalar a, gib_scalar b) {
  int sum_log;
  if (a == 0 || b == 0) return 0;
  sum_log = gib_gf_log[a] + gib_gf_log[b];
  if (sum_log >= (GIB_GALOIS_DEGREE-1)) sum_log -= (GIB_GALOIS_DEGREE-1);
  return gib_gf_ilog[sum_log];
}

gib_scalar gib_galois_div(gib_scalar a, gib_scalar b) {
  int diff_log;
  if (a == 0) return 0;
  if (b == 0) return -1;
  diff_log = gib_gf_log[a] - gib_gf_log[b];
  if (diff_log < 0) diff_log += (GIB_GALOIS_DEGREE-1);
  return gib_gf_ilog[diff_log];
}

int gib_galois_init() {
  static int run = 0;
  int i, j, b, log;
  if (run) return 0;
  run = 1;
  
  /* This polynomial (and its use) was given as an example in James Plank's 
   * tutorial on Reed-Solomon coding for RAID.
   */
  int prim_poly = GIB_GENERATOR;
  memset(gib_gf_ilog, 0, GIB_GALOIS_DEGREE);
  memset(gib_gf_log, 0, GIB_GALOIS_DEGREE);
  
  b = 1;
  for (log = 0; log < GIB_GALOIS_DEGREE; log++) {
    gib_gf_log[b] = (gib_scalar) log;
    gib_gf_ilog[log] = (gib_scalar) b;
    b = b << 1;
    if (b & GIB_GALOIS_DEGREE) b = b ^ prim_poly;
  }
  
  return 0;
}

int gib_galois_gen_F(gib_scalar *mat, int rows, int cols) {
  /* F forms the lower portion (m x n) of A */
  int i, j, rc;
  gib_scalar *tmpA = NULL;
  tmpA = (gib_scalar *)malloc((rows+cols)*(cols)*sizeof(gib_scalar));
  if (tmpA == NULL)
    return GIB_OOM;
  if ((rc = gib_galois_gen_A(tmpA, rows+cols, cols)))
    return rc;
  if (tmpA == NULL) return GIB_OOM;
  for (i = cols; i < rows+cols; i++)
    for (j = 0; j < cols; j++)
      mat[(i-cols)*cols+j] = tmpA[i*cols+j];
  free(tmpA);
  return 0;
}

int gib_galois_gen_A(gib_scalar *mat, int rows, int cols) {
  int i, j, p;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      mat[i*cols+j] = 1;
      for (p = 0; p < j; p++)
	mat[i*cols+j] = gib_galois_mul(mat[i*cols+j], i);
    }
  }
  gib_galois_gaussian_elim(mat, NULL, rows, cols);
  return 0;
}

int gib_galois_gaussian_elim(gib_scalar *mat, gib_scalar *inv, int rows, 
			     int cols) {
  /* If the caller wants an inverse, inv will be not null. */
  if (inv != NULL && rows != cols) {
    /* The system is overqualified, and needs to be reduced in order to
     * compute the inverse requested.
     */
    return GIB_ERR;
  }
  
  int i, j, e;
  if (inv != NULL) {
    /* Initialize to identity */
    for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
	inv[i*cols+j] = (i == j) ? 1 : 0;
  }
  
  for (i = 0; i < cols; i++) {
    /* Make sure A[i][i] is nonzero by swapping */
    if (mat[i*cols+i] == 0) {
      for (j = i+1; mat[i*cols+j] == 0 && j < cols; j++);
      for (e = 0; e < rows; e++) {
	int tmp = mat[e*cols+i];
	mat[e*cols+i] = mat[e*cols+j];
	mat[e*cols+j] = tmp;
	if (inv != NULL) {
	  tmp = inv[e*cols+i];
	  inv[e*cols+i] = inv[e*cols+j];
	  inv[e*cols+j] = tmp;
	}
      }
    }
    
    int inverse = gib_galois_div(1, mat[i*cols+i]);
    /* Make mat[i,i] == 1 by dividing down the column by mat[i,i] */
    for (e = 0; e < rows; e++) {
      mat[e*cols+i] = gib_galois_mul(inverse, mat[e*cols+i]);
      if (inv != NULL)
	inv[e*cols+i] = gib_galois_mul(inverse, inv[e*cols+i]);
    }
    
    /* Subtract a multiple of this column from all columns so that all
     * mat[i,j]==0 where i != j
     */
    for (j = 0; j < cols; j++) {
      int fij = mat[i*cols+j];
      if (j == i) continue;
      for (e = 0; e < rows; e++) {
	mat[e*cols+j] = mat[e*cols+j] ^ 
	  gib_galois_mul(fij, mat[e*cols+i]);
	if (inv != NULL) {
	  inv[e*cols+j] = inv[e*cols+j] ^
	    gib_galois_mul(fij, inv[e*cols+i]);
	}
      }
    }			
  }
  return 0;
}
