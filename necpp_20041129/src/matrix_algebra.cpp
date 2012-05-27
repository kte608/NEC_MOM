/*
	Copyright (C) 2004  Timothy C.A. Molteno
	
	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include "math_util.h"
#include <iostream>
using namespace std;

#include "matrix_algebra.h"
#include "nec_output.h"

#ifndef ATLAS
/*
	Subroutine to factor a matrix into a unit lower triangular matrix 
	and an upper triangular matrix using the Gauss-Doolittle algorithm 
	presented on pages 411-416 of A. Ralston -- a first course in 
	numerical analysis.
	
	Comments below refer to comments in Ralstons text.
	
	(matrix transposed.)
*/
void lu_decompose(nec_output_file& s_output, int n, complex_array& a, int_array& ip, int ndim)
{
	DEBUG_TRACE("lu_decompose(" << n << "," << ndim << ")");
	
	/* Allocate scratch memory */
	complex_array scm;
	scm.resize(n);
	
	/* Un-transpose the matrix for Gauss elimination */
	for (int i = 1; i < n; i++ )
	{
		int i_offset = i * ndim;
		int j_offset = 0;
		for (int j = 0; j < i; j++ )
		{
			nec_complex aij = a[i+j_offset];
			a[i+j_offset] = a[j+i_offset];
			a[j+i_offset] = aij;
			
			j_offset += ndim;
		}
	}
	
	bool iflg=false;
	/* step 1 */
	for (int r = 0; r < n; r++ )
	{
		int r_offset = r*ndim;
		
		for (int k = 0; k < n; k++ )
			scm[k]= a[k+r_offset];
		
		/* steps 2 and 3 */
		int rm1 = r;
		for (int j = 0; j < rm1; j++ )
		{
			int pj= ip[j]-1;
			nec_complex arj = scm[pj];
			a[j+r_offset]= arj;
			scm[pj]= scm[j];
			int jp1 = j+1;
		
			int j_offset = j*ndim;
			for (int i = jp1; i < n; i++ )
				scm[i] -= a[i+j_offset]* arj;
		
		}
		
		/* step 4 */
		nec_float dmax = norm(scm[r]);
		
		int rp1 = r+1;
		ip[r]= rp1;
		for (int i = rp1; i < n; i++ )
		{
			nec_float elmag = norm(scm[i]);
			if ( elmag >= dmax)
			{
				dmax = elmag;
				ip[r] = i+1;	// set the permute array element
			}
		}
		
		if ( dmax < 1.e-10)
			iflg=true;
		
		int pr = ip[r]-1;
		a[r+r_offset] = scm[pr];
		scm[pr] = scm[r];
		
		/* step 5 */
		if ( rp1 < n)
		{
			nec_complex arr = cplx_10() / a[r+r_offset];
			
			for (int i = rp1; i < n; i++ )
				a[i+r_offset]= scm[i]* arr;
		}
		
		if ( true == iflg )
		{
			s_output.string("\n  PIVOT(");
			s_output.integer(r);
			s_output.string(")= ");
			s_output.real(dmax);
			iflg=false;
		}	
	} /* for( r=0; r < n; r++ ) */
	

#ifdef _NEC_ERROR_CHECK

	// Debug output to try and figure out the LAPACK stuff
	for (int j = 0; j < n; j++ )
	{
		cout << "ip[" << j <<"] = " << ip[j] << endl;
		for (int i = 0; i < n; i++ )
		{
			cout << a[i+j*ndim] << ", ";
		}
		cout << endl;
	}
#endif
	
}

#else
/* 
	Subroutine to factor a matrix into a unit lower triangular matrix 
	and an upper triangular matrix using ATLAS (matrix transposed.)

	
	int clapack_zgetrf(const enum CBLAS_ORDER Order, const int M, const int N,
                   void *A, const int lda, int *ipiv);
		   
		CblasRowMajor 	assumes we are storing A in C style
		CblasColMajor	assume we are storing A in FORTRAN style
		   
	ipiv	(output) INTEGER array, dimension (min(M,N))   
		The pivot indices; for 1 <= i <= min(M,N), row i of the   
		matrix was interchanged with row IPIV(i).
	
	input parameters
		ndim		actual full matrix dimension
		n		dimension of submatrix to be factored
		a_in[ndim,ndim]	The full input matrix
		ip[ndim]	The pivot points
		

 */
extern "C"
{
#include <atlas_enum.h>
#include <clapack.h>
}

void lu_decompose(nec_output_file& s_output,  int n, complex_array& a_in, int_array& ip, int ndim)
{
	DEBUG_TRACE("factor_lapack(" << n << "," << ndim << ")");
	ASSERT(n <= ndim);
	
	// copy the input matrix a_in into a temporary array.
	// transposing as we go... Should use cblas_zgemm...
	complex_array A(n,n);
	int_array piv(n);
	
	for (int i = 0; i < n; i++ )
	{
		int row_start = i*ndim;
		for (int j = 0; j < n; j++ )
		{
			nec_complex& x = a_in[row_start + j];
			A.set(i,j,x);
		}
	}

	int lead_dim = std::max(1, n);
		
	// Now call the LAPACK LU-Decomposition
	int info = clapack_zgetrf (CblasRowMajor, n, n, 
		(void*) A.get_ptr(), lead_dim, piv.get_ptr());
	
	if (0 != info)
	{
		/*
			The factorization has been completed, but the factor U is exactly singular,
			and division by zero will occur if it is used to solve a system of equations. 
		*/
		throw new nec_exception("nec++: LU Decomposition Failed: %d\n", info);
		exit(1);
	}
		
	
	/*
	IPIV	(output) INTEGER array, dimension (min(M,N))   
		The pivot indices; for 1 <= i <= min(M,N), row i of the   
		matrix was interchanged with row IPIV(i).
	*/  
	for (int j = 0; j < n; j++ )
	{
		ip[j] = piv[j] + 1;
	}
		
	// copy the output back into the a_in array.
	for (int i = 0; i < n; i++ )
	{
		int row_start = i*ndim;
		
		for (int j = 0; j < n; j++ )
		{
			a_in[row_start + j] = A.get(i,j);
		}
	}
	

#ifdef _NEC_ERROR_CHECK

	for (int i = 0; i < n; i++ )
	{
		cout << "piv[" << i << "]" << ip[i] << ",";
	}
	
	for (int i = 0; i < n; i++ )
	{
		int row_start = i*ndim;
		
		for (int j = 0; j < n; j++ )
		{
			cout << a_in[row_start + j] << ", ";
		}
		cout << endl;
	}
#endif
} 
#endif


/*-----------------------------------------------------------------------*/

/*	factrs

	For symmetric structure, transforms submatricies to form
	matricies of the symmetric modes and calls routine to LU decompose
	matricies.
	
	If no symmetry [nrow = np], the routine is called to LU decompose the
	complete matrix.
*/
void factrs(nec_output_file& s_output,  int np, int nrow, complex_array& a, int_array& ip )
{
	DEBUG_TRACE("factrs(" << np << "," << nrow << ")");
	if (nrow == np) // no symmetry
	{
		lu_decompose(s_output,  np, a, ip, nrow );
		return;
	}
	
	int num_symmetric_modes = nrow / np;
	DEBUG_TRACE("\tnum_symmetric_modes = " << num_symmetric_modes);
	
	for (int mode = 0; mode < num_symmetric_modes; mode++ )
	{
		int mode_offset = mode * np;
		
		complex_array a_temp = a.sub_array(mode_offset);
		int_array ip_temp = ip.sub_array(mode_offset);
		
		lu_decompose(s_output,  np, a_temp, ip_temp, nrow );
	}
}

/*-----------------------------------------------------------------------*/

/* 	
	Subroutine to solve the matrix equation lu*x=b where l is a unit
	lower triangular matrix and u is an upper triangular matrix both
	of which are stored in a.  the rhs vector b is input and the
	solution is returned through vector b.   (matrix transposed)
      
	  COMPLEX*16 A,B,Y,SUM
      INTEGER PI
      COMMON /SCRATM/ Y(2*MAXSEG)
      DIMENSION A(NDIM,NDIM), IP(NDIM), B(NDIM)
C
C     FORWARD SUBSTITUTION
C
      DO 3 I=1,N
      PI=IP(I)
      Y(I)=B(PI)
      B(PI)=B(I)
      IP1=I+1
      IF (IP1.GT.N) GO TO 2
      DO 1 J=IP1,N
      B(J)=B(J)-A(J,I)*Y(I)
1     CONTINUE
2     CONTINUE
3     CONTINUE
C
C     BACKWARD SUBSTITUTION
C
      DO 6 K=1,N
      I=N-K+1
      SUM=(0.,0.)
      IP1=I+1
      IF (IP1.GT.N) GO TO 5
      DO 4 J=IP1,N
      SUM=SUM+A(I,J)*B(J)
4     CONTINUE
5     CONTINUE
      B(I)=(Y(I)-SUM)/A(I,I)
6     CONTINUE
      RETURN
      END
	
*/
void solve( int n, complex_array& a, int_array& ip,
    complex_array& b, int ndim )
{
	DEBUG_TRACE("solve(" << n << "," << ndim << ")");
/*	DEBUG_TRACE("a.size=(" << a.size() << ")");
	DEBUG_TRACE("ip.size=(" << ip.size() << ")");
	DEBUG_TRACE("b.size=(" << b.size() << ")");
	ASSERT(a.size() >= n*ndim);
	ASSERT(ip.size() == ndim);
	ASSERT(b.size() >= ndim); // b can be greater since it was assigned a dimension of n+m somewhere (in netwk())
*/	
	complex_array y;
	y.resize(n);
	
	/* forward substitution */
	for (int i = 0; i < n; i++ )
	{
		int pia= ip[i]-1;
		y[i]= b[pia];
		b[pia]= b[i];
		int ip1= i+1;
		
		int i_offset = i*ndim;
		for (int j = ip1; j < n; j++ )
			b[j] -= a[j+i_offset] * y[i];
	}
	
	/* backward substitution */
	for (int k = 0; k < n; k++ )
	{
		int i= n-k-1;
		nec_complex sum(cplx_00());
		int ip1= i+1;
		
		for (int j = ip1; j < n; j++ )
			sum += a[i+j*ndim]* b[j];
		
		b[i] = (y[i]- sum) / a[i+i*ndim];
	}
}



/* subroutine solves, for symmetric structures, handles the */
/* transformation of the right hand side vector and solution */
/* of the matrix eq. */
void solves(complex_array& a, int_array& ip, complex_array& b, int neq,
	int nrh, int np, int n, int mp, int m, int nop, 
	complex_array& symmetry_array
)
{
	DEBUG_TRACE("solves(" << neq << "," << nrh << "," << np << "," << n << ")");
	int  ic, i, kk, ia, ib, j, k;

	nec_complex  sum;
	
	/* Allocate to scratch memory */
	complex_array scm;
	scm.resize(n + 2*m);
	
	int npeq= np+ 2*mp;
	nec_float fnop = nop;
	nec_float fnorm = 1.0/ fnop;
	int nrow= neq;
	
	if ( nop != 1)
	{
		for( ic = 0; ic < nrh; ic++ )
		{
			if ( (n != 0) && (m != 0) )
			{
				for( i = 0; i < neq; i++ )
					scm[i]= b[i+ic*neq];
			
				kk=2* mp;
				ia= np-1;
				ib= n-1;
				j= np-1;
			
				for( k = 0; k < nop; k++ )
				{
					if ( k != 0 )
					{
						for( i = 0; i < np; i++ )
						{
							ia++;
							j++;
							b[j+ic*neq]= scm[ia];
						}
				
						if ( k == (nop-1) )
							continue;
					} /* if ( k != 0 ) */
				
					for( i = 0; i < kk; i++ )
					{
						ib++;
						j++;
						b[j+ic*neq]= scm[ib];
					}
				} /* for( k = 0; k < nop; k++ ) */
			
			} /* if ( (n != 0) && (m != 0) ) */
		
			/* transform matrix eq. rhs vector according to symmetry modes */
			for( i = 0; i < npeq; i++ )
			{
				for( k = 0; k < nop; k++ )
				{
					ia= i+ k* npeq;
					scm[k]= b[ia+ic*neq];
				}
			
				sum= scm[0];
				for( k = 1; k < nop; k++ )
					sum += scm[k];
			
				b[i+ic*neq]= sum* fnorm;
			
				for( k = 1; k < nop; k++ )
				{
					ia= i+ k* npeq;
					sum= scm[0];
				
					for( j = 1; j < nop; j++ )
						sum += scm[j]* conj( symmetry_array[k+j*nop]);
				
					b[ia+ic*neq]= sum* fnorm;
				}
			} /* for( i = 0; i < npeq; i++ ) */
		
		} /* for( ic = 0; ic < nrh; ic++ ) */
	
	} /* if ( nop != 1) */
	
	/* solve each mode equation */
	for( kk = 0; kk < nop; kk++ )
	{
		ia= kk* npeq;
		ib= ia;
	
		for( ic = 0; ic < nrh; ic++ )
		{
			complex_array a_sub = a.sub_array(ib);
			complex_array b_sub = b.sub_array(ia+ic*neq);
			int_array ip_sub = ip.sub_array(ia);
			solve( npeq, a_sub, ip_sub, b_sub, nrow );
		}
	
	} /* for( kk = 0; kk < nop; kk++ ) */
	
	if ( nop == 1)
	{
		return;
	}
	
	/* inverse transform the mode solutions */
	for( ic = 0; ic < nrh; ic++ )
	{
		for( i = 0; i < npeq; i++ )
		{
			for( k = 0; k < nop; k++ )
			{
				ia= i+ k* npeq;
				scm[k]= b[ia+ic*neq];
			}
		
			sum= scm[0];
			for( k = 1; k < nop; k++ )
				sum += scm[k];
		
			b[i+ic*neq]= sum;
			for( k = 1; k < nop; k++ )
			{
				ia= i+ k* npeq;
				sum= scm[0];
			
				for( j = 1; j < nop; j++ )
					sum += scm[j]* symmetry_array[k+j*nop];
			
				b[ia+ic*neq]= sum;
			}
		
		} /* for( i = 0; i < npeq; i++ ) */
	
		if ( (n == 0) || (m == 0) )
			continue;
	
		for( i = 0; i < neq; i++ )
			scm[i]= b[i+ic*neq];
	
		kk=2* mp;
		ia= np-1;
		ib= n-1;
		j= np-1;
	
		for( k = 0; k < nop; k++ )
		{
			if ( k != 0 )
			{
				for( i = 0; i < np; i++ )
				{
					ia++;
					j++;
					b[ia+ic*neq]= scm[j];
				}
			
				if ( k == nop)
					continue;
			
			} /* if ( k != 0 ) */
		
			for( i = 0; i < kk; i++ )
			{
				ib++;
				j++;
				b[ib+ic*neq]= scm[j];
			}
		} /* for( k = 0; k < nop; k++ ) */
	
	} /* for( ic = 0; ic < nrh; ic++ ) */
}

/*-----------------------------------------------------------------------*/
/* test for convergence in numerical integration */
void test(
	nec_float f1r, nec_float f2r, nec_float *tr,
 	nec_float f1i, nec_float f2i, nec_float *ti,
	nec_float dmin )
{
	static nec_float _min_val =  1.0e-37;
	
/*
{
  double den;

  den= fabs( f2r);
  *tr= fabs( f2i);

  if( den < *tr)
    den= *tr;
  if( den < dmin)
    den= dmin;

  if( den < 1.0e-37)
  {
    *tr=0.;
    *ti=0.;
    return;
  }

  *tr= fabs(( f1r- f2r)/ den);
  *ti= fabs(( f1i- f2i)/ den);
}

*/	
	nec_float den = fabs( f2r);
	nec_float temp_tr = fabs( f2i);
	
	if( den < temp_tr)
		den = temp_tr;
	if( den < dmin)
		den = dmin;
	
	if( den < _min_val)
	{
		*tr = 0.0;
		*ti = 0.0;
		return;
	}
	
	*tr= fabs((f1r - f2r)/ den);
	*ti= fabs((f1i - f2i)/ den); 
}

/*
	Simpler test for convergence in numerical integration.
	This tests only one number. It is a special case of the
	test() function above.
*/
nec_float test_simple( nec_float f1r, nec_float f2r, nec_float dmin )
{
	static nec_float _min_val =  1.0e-37;
	
	nec_float den = fabs(f2r);
	
	if( den < dmin)
		den = dmin;
	if (den < _min_val)
	{
		return 0.0;
	}
	
	return fabs((f1r - f2r) / den);	
}
