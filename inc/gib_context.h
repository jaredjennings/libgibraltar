/* Author:  Matthew Curry
 * Email:   mlcurry@sandia.gov
 *
 * Defines the Gibraltar context object, which contains the properties, matrix,
 * and accelerator context.  It is used for most Gibraltar calls.
 */
#ifndef GIB_CONTEXT_H_
#define GIB_CONTEXT_H_

/*
#define GIB_GALOIS_DEGREE 65536
static const int GIB_GENERATOR = 0x1100b;
typedef unsigned short gib_scalar;
*/
#define GIB_GALOIS_DEGREE 256
static const int GIB_GENERATOR = 0435;
typedef unsigned char gib_scalar;

struct gib_context_t {
	int n, m;
	gib_scalar *F;
	/* The stuff below is only used in the GPU case */
	void *acc_context;
};

typedef struct gib_context_t* gib_context;

#endif /*GIB_CONTEXT_H_*/
