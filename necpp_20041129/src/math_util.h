#ifndef __math_util__
#define __math_util__

/*
	Various Useful Math Utilities for nec2cpp
	
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
 
//#include <cmath>
#include <complex>

#include "safe_array.h"

/**
	This typedef allows us to use nec2++ with
	a different complex number precision. For example
	float, or long double.
*/
typedef double nec_float;
typedef std::complex<nec_float> nec_complex;

typedef safe_array<int> int_array;
typedef safe_array<nec_float>  real_array;
typedef safe_array<nec_complex>  complex_array;


inline nec_complex cplx_00()
{
	static nec_complex _cplx00(0.0,0.0); return _cplx00;
}

inline nec_complex cplx_01()
{
	static nec_complex _cplx01(0.0,1.0); return _cplx01;
}

inline nec_complex cplx_10()
{
	static nec_complex _cplx10(1.0,0.0); return _cplx10;
}

inline nec_complex cplx_11()
{
	static nec_complex _cplx11(1.0,1.0); return _cplx11;
}


inline nec_complex cplx_exp(nec_float x)
{
	return nec_complex(cos(x),sin(x));
}


inline nec_float pi()
{
	static nec_float _pi = 3.1415926536; return _pi;
}

inline nec_float two_pi()
{
	static nec_float _tmp = 2.0 * pi(); return _tmp;
}

inline nec_float pi_two()
{
	static nec_float _tmp = pi() / 2.0; return _tmp;
}

inline nec_complex two_pi_j()
{
	static nec_complex _tmp(0.0,two_pi()); return _tmp;
}




inline nec_float rad_to_degrees(nec_float in_radians)
{
	static nec_float _rad_to_deg = 360.0 / (2 * pi()); // 57.29577951
	
	return in_radians * _rad_to_deg;
}

inline nec_float degrees_to_rad(nec_float in_degrees)
{
	static nec_float _deg_to_rad = (2 * pi()) / 360.0;
	
	return in_degrees * _deg_to_rad;
}

/**
	Create a complex number from a magnitude and
	an angle in degrees.
*/
inline nec_complex deg_polar(nec_float r, nec_float theta)
{
	return std::polar(r, degrees_to_rad(theta));
}


/**
	Get the angle of a complex number in degrees.
*/
inline nec_float arg_degrees(nec_complex z)
{
	return rad_to_degrees(arg(z));
}


/**
	atgn2 is arctangent function modified to return 0 when x=y=0.
*/
inline nec_float atgn2( nec_float x, nec_float y)
{
	if ((0.0 == y) && (0.0 == x))
		return 0.0;
		
	return( std::atan2(y, x) );
}


/** function db10 returns db for magnitude (field) */
inline nec_float db10( nec_float x )
{
	if ( x < 1.0e-20 )
		return( -999.99 );
	
	return( 10.0 * log10(x) );
}

/*-----------------------------------------------------------------------*/

/** function db20 returns db for mag**2 (power) i */
inline nec_float db20( nec_float x )
{
	if ( x < 1.0e-20 )
		return( -999.99 );
	
	return( 20.0 * log10(x) );
}



#endif /* __math_util__ */
