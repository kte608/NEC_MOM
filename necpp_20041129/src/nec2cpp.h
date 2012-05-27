/*
 * nec2cpp.h - header file for nec2 c++
 */

#include "math_util.h"

using namespace std;

int 	main(int argc, char **argv);
void 	fblock(int nrow, int ncol, int imax, int ipsym);

void 	readmn(FILE* input_fp, FILE* output_fp, char *gm, int *i1, int *i2, int *i3, int *i4, nec_float *f1,
	nec_float *f2, nec_float *f3, nec_float *f4, nec_float *f5, nec_float *f6);

#include "misc.h"


