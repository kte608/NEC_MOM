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
#ifndef __nec_radiation_pattern__
#define __nec_radiation_pattern__

#include "nec_results.h"
#include "math_util.h"

class nec_context;

class nec_radiation_pattern : public nec_base_result
{
public:
	// Radiation Pattern
	nec_radiation_pattern(nec_context* in_context, int in_n_theta, int in_n_phi, int in_ifar, nec_float in_wavelength, nec_float pinr, nec_float pnlr)
	{
		n_theta = in_n_theta;
		n_phi = in_n_phi;
		
		int n_angles = n_theta * n_phi;
		
		_gain.resize(n_angles);
		_power_gain_vert.resize(n_angles);
		_power_gain_horiz.resize(n_angles);
		_power_gain_tot.resize(n_angles);
		_polarization_axial_ratio.resize(n_angles);
		_polarization_tilt.resize(n_angles);
		_polarization_sense_index.resize(n_angles);
				
		_e_theta.resize(n_angles);
		_e_phi.resize(n_angles);
		_e_r.resize(n_angles);
		
		_ifar = in_ifar;
		_wavelength = in_wavelength;
		_pinr = pinr;
		_pnlr = pnlr;
		
		m_context = in_context;
		m_analysis_done = false;
		_maximum_gain = -999.0;
	}
		
	virtual void write_to_file(ostream& os)
	{
		analyze();
		write_to_file2(os);
	}

	void analyze();
	
	void write_gain_normalization()
	{
		if (_ifar != 1)
		{
			nec_float norm = get_maximum_gain_db();
			printf("Max Gain: %f\n",norm);
		}
	}
	
	nec_float get_maximum_gain_db()
	{
		return get_gain_normalization_factor(0);
	}
	
	/*! \brief Get an array representing restuls of the radiation pattern
	*/
	real_array get_radiation_pattern()
	{
		return _gain;
	}

private:
	nec_context* m_context;
	bool m_analysis_done;
	
	real_array _gain;
	int n_theta, n_phi;
	
	int _ifar;
	
	nec_float _wavelength;
	nec_float _pinr;
	nec_float _pnlr;
	
	// quantitles...
	nec_float _average_power_gain;
	nec_float _average_power_solid_angle;
	nec_float _maximum_gain;

	real_array	_power_gain_vert;
	real_array	_power_gain_horiz;
	real_array	_power_gain_tot;
	real_array	_polarization_axial_ratio;
	real_array	_polarization_tilt;
	int_array	_polarization_sense_index;
				
	complex_array	_e_theta;
	complex_array	_e_phi;
	complex_array	_e_r;
	void write_to_file2(ostream& os);
	
	nec_float get_gain_normalization_factor(nec_float gnor);
	
	void write_normalized_gain(ostream& os, nec_float gnor, int inor);
};

// some auxiliary functions to be made private once
// the radiation pattern calculation is done entirely
// inside this class...
void gfld(nec_context* m_context, nec_float rho, nec_float phi, nec_float rz,
    nec_complex *eth, nec_complex *epi,
    nec_complex *erd, nec_complex ux, int ksymp, nec_float _wavelength  );

void ffld(nec_context* m_context, nec_float thet, nec_float phi,
    nec_complex *eth, nec_complex *eph, nec_float _wavelength  );

void fflds(nec_context* m_context, nec_float rox, nec_float roy, nec_float roz,
    complex_array& scur, nec_complex *ex,
    nec_complex *ey, nec_complex *ez );

#endif /* __nec_radiation_pattern__ */

