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
#ifndef __nec_ground__
#define __nec_ground__

#include "math_util.h"

/*
	\section{NEC Ground Specification}
	
	\subsection{the GN Card}
	
	The GN card specifies the relative dielectric constant and conductivity of ground in the vicinity of the antenna. 
	In addition, a second set of ground parameters for a second medium can be specified,
	or a radial wire ground screen can be modeled using a reflection coefficient approximation.
	
	The parameters of the second medium can also be specified on another data card whose mnemonic is GD.
	With the GD card, the parameters of the second medium can be varied and only the radiated fields
	need to be recalculated.
	Furthermore, if a radial wire ground screen has been specified on the GN card,
	the GD card is the only way to include a second medium.
	See Section~\ref{sec_card_gd} for details.
	
	\subsection{"gd" card: ground representation}
	\label{sec_card_gd}
	
	To specify the ground parameters of a second medium which is not in the immediate vicinity of the antenna.
	This card may only be used if a GN card has also been used.
	It does not affect the fields of surface patches.
	
	<ground type=finite>
		<dielectric_constant = 80/>
		<conductivity = 4/>
		<radial-wire>
			<count=12>
			<length=15 unit=cm>
			<
		</radial-wire>
	</ground>
	
*/
class nec_ground
{
public:

	nec_ground()
	{
		default_values();
	}
	
	void default_values()
	{
		ksymp=1;
		radial_wire_count=0;
		iperf=0;
	}
	/*
		Parse a GN card. The input parameters here are the fields of the
		GN card.	
	*/
	void parse_gn(int itmp1, int itmp2,
		nec_float tmp1, nec_float tmp2,
		nec_float tmp3, nec_float tmp4,
		nec_float tmp5, nec_float tmp6
		);
	
	/*
		Setup a cliff (two medium ground)
	*/
	void setup_cliff(nec_float in_eprs2,
		nec_float in_sig2,
		nec_float clt, nec_float cht);
	
	nec_complex get_zrati2(nec_float wavelength);
	
	nec_float get_cl(nec_float wavelength) // cliff edge in wavelengths, 
	{
		return cliff_edge_distance / wavelength;
	}
	
	nec_float get_ch(nec_float wavelength) // cliff Height in wavelengths.
	{
		return cliff_height / wavelength;
	}
	 
	
	// accessors for the ground type
	inline bool type_finite_reflection()	{	return (0 == iperf); }
	inline bool type_perfect()				{	return (1 == iperf); }
	inline bool type_sommerfeld_norton()	{	return (2 == iperf); }
	
	
	bool is_valid()
	{
		if (iperf < 0) return false;
		if (iperf > 2) return false;
		
		return true;
	}
	
	int ksymp;
	
	
	nec_float epsr;	// relative dielectric constant
	nec_float sig;	// Conductivity
	
	// radial wire ground
	int radial_wire_count;
	nec_float radial_wire_length;
	nec_float radial_wire_radius;
	
	// second medium parameters
	nec_float cliff_edge_distance;
	nec_float cliff_height;
	nec_float epsr2;	// Relative dielectric constant
	nec_float sig2;		// Conductivity in mhos/meter

private:
	/* iperf: Ground-type flag. The options are: 
		-1 - nullifies ground parameters previously used and sets free- space condition. The remainder of the card is left blank in this case. 
		O - finite ground, reflection-coefficient approximation. 
		1 - perfectly conducting ground. 
		2 - finite ground, Sommerfeld/Norton method.
	*/
	int iperf;
};

#endif /* __nec_ground__ */

