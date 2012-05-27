/***************************************************************************
 *   Copyright (C) 2004 by Tim Molteno                                     *
 *   tim@physics.otago.ac.nz                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/*! \file nec.h
    \brief nec++ Library Functions.
    
\verbatim
How to use libNEC. 
	
	Enter the following file into test_nec.c, and compile with
	gcc -o test_nec test_nec.c -L . -lnecpp -lm -lstdc++
	
	#include "nec.h"
	#include <stdio.h>
	
	int main(int argc, char **argv)
	{
		nec_context* nec;
		double gain;
			
		nec = nec_create();
		nec_wire(nec, 0, 36, 0, 0, 0, -0.042, 0.008, 0.017, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 21, -0.042, 0.008, 0.017, -0.048, 0.021, -0.005, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 70, -0.048, 0.021, -0.005, 0.039, 0.032, -0.017, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 70, -0.048, 0.021, -0.005, 0.035, 0.043, 0.014, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 50, -0.042, 0.008, 0.017, 0.017, -0.015, 0.014, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 66, 0.017, -0.015, 0.014, -0.027, 0.04, -0.031, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 85, -0.027, 0.04, -0.031, 0.046, -0.01, 0.028, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 47, 0.046, -0.01, 0.028, -0.013, -0.005, 0.031, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 70, 0.017, -0.015, 0.014, -0.048, -0.038, -0.04, 0.001, 1.0, 1.0);
		nec_wire(nec, 0, 77, -0.048, -0.038, -0.04, 0.049, -0.045, -0.04, 0.001, 1.0, 1.0);
		nec_geometry_complete(nec, 0, 0);
		
		nec_gn_card(nec, -1,0,0.0, 0.0, 0.0,0.0, 0.0, 0.0);
		nec_ld_card(nec, 5,0,0,0,3.72e7,0.0,0.0);
		nec_pt_card(nec, -1, 0, 0, 0);
		nec_ex_card(nec, 1, 1, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
		nec_fr_card(nec, 0, 2, 2400.0, 100.0);
		nec_rp_card(nec, 0, 1, 1, 500, 90.0, 90.0, 0.0, 0.0, 0.0, 0.0);
		gain = nec_get_maximum_gain(nec);
		
		nec_delete(nec);
		
		printf("Gain is %f",gain);
		
		return 0;
	}
\endverbatim
*/

#ifndef __nec__
#define __nec__

/*! A Struct to represent the nec_context class.*/
struct nec_context;
typedef struct nec_context nec_context;

#ifdef __cplusplus
extern "C" {
#endif

/*! Construct and initialize an nec_context */
nec_context* nec_create();
/*! Delete an nec_context object. */
void nec_delete(nec_context* in_context);


/*
	Geometry description functions.
*/
void nec_wire(nec_context* in_context, int tag_id, int segment_count,
		double xw1, double yw1, double zw1,
		double xw2, double yw2, double zw2, 
		double rad, double rdel, double rrad);
void nec_geometry_complete(nec_context* in_context, int card_int_1, int card_int_2);


/*
	NEC card functions.
*/

/*!
 * FR crd
 *	@param in_context The nec_context created with nec_create()
 *	@param in_ifrq 0 is a linear range of frequencies, 1 is a log range.
 *	@param in_nfrq The number of frequencies
 *	@param in_freq_mhz The starting frequency in MHz.
 *	@param in_del_freq The frequency step (in MHz for ifrq = 0)
 */
void nec_fr_card(nec_context* in_context, int in_ifrq, int in_nfrq, double in_freq_mhz, double in_del_freq);
void nec_ld_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3);
void nec_gn_card(nec_context* in_context, int itmp1, int itmp2, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
void nec_ex_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
void nec_tl_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
void nec_nt_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
void nec_xq_card(nec_context* in_context, int itmp1);
void nec_gd_card(nec_context* in_context, double tmp1, double tmp2, double tmp3, double tmp4);
void nec_rp_card(nec_context* in_context, int itmp1, int n_theta, int n_phi, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
void nec_pt_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4);
void nec_pq_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4);
void nec_kh_card(nec_context* in_context, double tmp1);
void nec_ne_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
void nec_nh_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4, double tmp1, double tmp2, double tmp3, double tmp4, double tmp5, double tmp6);
void nec_ek_card(nec_context* in_context, int itmp1);
void nec_cp_card(nec_context* in_context, int itmp1, int itmp2, int itmp3, int itmp4);
void nec_pl_card(nec_context* in_context, char* ploutput_filename, int itmp1, int itmp2, int itmp3, int itmp4);


double nec_get_maximum_gain(nec_context* in_context);



#ifdef __cplusplus
}
#endif


#endif
