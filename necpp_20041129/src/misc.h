#ifndef __MISC__
#define __MISC__

#include <stdio.h>
#include "math_util.h"

/* carriage return and line feed */
#define	CR	0x0d
#define	LF	0x0a

/* max length of a line read from input file */
#define	LINE_LEN	132

/*  usage()
 *
 *  prints usage information
 */
void usage(void);


/*  abort_on_error()
 *
 *  prints an error message and exits
 */
void abort_on_error( int why );

/* Returns process time (user+system) BUT in _msec_ */
void secnds( nec_float *x);


/* Does the STOP function of fortran but with return value */
int stop( int flag );

/*  load_line()
 *
 *  loads a line from a file, aborts on failure. lines beginning
 *  with a '#' are ignored as comments. at the end of file EOF is
 *  returned.
 */
int load_line( char *buff, FILE *pfile );



#endif /* __MISC__ */
