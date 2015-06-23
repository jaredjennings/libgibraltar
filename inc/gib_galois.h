/* Author:  Matthew Curry
 * Email:   mlcurry@sandia.gov
 *
 * Internal CPU-based Galois field arithmetic interface.
 */
#ifndef GIB_GALOIS_H_
#define GIB_GALOIS_H_

#include "gib_context.h"

#ifdef __cplusplus
extern "C" {
#endif

extern gib_scalar gib_gf_log[GIB_GALOIS_DEGREE];
extern gib_scalar gib_gf_ilog[GIB_GALOIS_DEGREE];
int gib_galois_init();
int gib_galois_gen_F(gib_scalar *mat, int rows, int cols);
int gib_galois_gen_A(gib_scalar *mat, int rows, int cols);
int gib_galois_gaussian_elim(gib_scalar *mat, gib_scalar *inv, int rows, 
		int cols);
#ifdef __cplusplus
}
#endif

#endif /*GIB_GALOIS_H_*/
