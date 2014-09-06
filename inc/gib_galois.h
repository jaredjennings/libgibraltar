/* Author:  Matthew Curry
 * Email:   mlcurry@sandia.gov
 *
 * Internal CPU-based Galois field arithmetic interface.
 */
#ifndef GIB_GALOIS_H_
#define GIB_GALOIS_H_

#ifdef __cplusplus
extern "C" {
#endif

extern unsigned char gib_gf_log[256];
extern unsigned char gib_gf_ilog[256];
extern unsigned char gib_gf_table[256][256];
int gib_galois_init();
int gib_galois_gen_F(unsigned char *mat, int rows, int cols);
int gib_galois_gen_A(unsigned char *mat, int rows, int cols);
int gib_galois_gaussian_elim(unsigned char *mat, unsigned char *inv, int rows, 
		int cols);
#ifdef __cplusplus
}
#endif

#endif /*GIB_GALOIS_H_*/
