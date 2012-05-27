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
#ifndef __nec_results__
#define __nec_results__
/*
	This class contains the results of the NEC analysis. Set methods
	will store the results in this class, but will NOT print the results
	to a file
*/
#include <vector>
#include <ostream>
#include <iostream>
#include <iomanip>
#include <string>

#include "math_util.h"
using namespace std;


class nec_base_result
{
private:
	bool _write_file;

public:
	virtual void write_to_file(ostream& os) = 0;
	
	nec_base_result()
		: _write_file(true)
	{
	}
	
	virtual ~nec_base_result()
	{
	}

	inline void section_start(ostream& os)
	{
		os << endl << endl << endl;
	}
	
	inline void int_out(ostream& os, int w, int i)
	{
		os << setw(w) << i << " ";
	}
	
	inline void string_out(ostream& os, int w, char* s)
	{
		os << right <<setw(w) << s << " ";
	}
	
	inline void real_out(ostream& os, int w, int p, nec_float f, bool sci = true)
	{
		ios::fmtflags flags = ios::showpoint | ios::uppercase | ios::right;
		if (sci)
			flags |= ios::scientific;
		else
			flags |= ios::fixed;
		
		os.unsetf(ios::adjustfield | ios::basefield | ios::floatfield);
		os.setf(flags);
		os.precision(p);
		os.width(w);
		os << f << " ";
	}
	
	inline void complex_out(ostream& os, int w, int p, nec_complex c, bool sci = true)
	{
		real_out(os,w,p,real(c),sci);
		real_out(os,w,p,imag(c),sci);
	}

	inline void polar_out(ostream& os, int w, int p, nec_complex c, bool sci = true)
	{
		real_out(os,w,p,abs(c),sci);
		real_out(os,w,p,arg_degrees(c),sci);
	}
	
	inline bool write_file() const
	{
		return _write_file;
	}
	
	inline void set_write_file(bool f)
	{
		_write_file = f;
	}
};



class nec_rx_pattern : public nec_base_result
{
	// Receiving Pattern
	nec_float _norm_factor;
	nec_float _eta, _axial_ratio;
	int _segment_number;
	string _type;
	
	long n_theta;
	long n_phi;
	nec_float _theta0, _theta_step;
	nec_float _phi0, _phi_step;
	
	real_array _mag;
	
public:
	nec_rx_pattern(
		int in_n_theta, int in_n_phi,
		real_array& in_mag,
		nec_float theta0, nec_float theta_step,
		nec_float phi0, nec_float phi_step,
		nec_float in_eta, 
		nec_float in_axial_ratio, 
		int in_segment_number, 
		string in_type)
	{
		n_theta = in_n_theta;
		n_phi = in_n_phi;
		
		_mag.copy(in_mag);
		_mag.resize(n_theta, n_phi);
		
		_theta0 = theta0;
		_theta_step = theta_step;
		
		_phi0 = phi0;
		_phi_step = phi_step;
		
		_eta = in_eta;
		_axial_ratio = in_axial_ratio;
		_segment_number = in_segment_number;
		_type = in_type;
		
		_mag.resize(n_theta, n_phi);
	}

	virtual ~nec_rx_pattern()
	{
	}
	
	void set_input(int theta_index, int phi_index, nec_float mag)
	{
		_mag.set(theta_index,phi_index,mag);
	}
	
	nec_float get_norm_factor()
	{
		return _mag.max();
	}
	
	virtual void write_to_file(ostream& os)
	{
		if (n_theta == 0)
			return;
		if (n_phi == 0)
			return;
			
		nec_float norm_factor = get_norm_factor();
		
		section_start(os);
		os << "                      ---- NORMALIZED RECEIVING PATTERN ----" << endl;
		os << "                      NORMALIZATION FACTOR: ";real_out(os,11,4,norm_factor);os << endl;
		os << "                      ETA: ";real_out(os,7,2,_eta,false); os << "DEGREES" << endl;
		os << "                      TYPE: " << _type << endl;
		os << "                      AXIAL RATIO: "; real_out(os,6,3,_axial_ratio,false); os << endl;
		os << "                      SEGMENT No: "; int_out(os, 5, _segment_number); os << endl << endl;
		os << "                      THETA     PHI       ---- PATTERN ----" << endl;
		os << "                      (DEG)    (DEG)       DB     MAGNITUDE" << endl;
		
		nec_float theta = _theta0;
		
		for (int t=0; t<n_theta; t++)
		{
			nec_float phi = _phi0;
			
			for (int p=0; p<n_phi;p++)
			{
				nec_float magnitude = _mag.get(t,p) / norm_factor;
				nec_float gain = db20(magnitude);
				
				os << "                    ";
				real_out(os,7,2, theta, false);
				real_out(os,7,2, phi, false);
				os << "  ";
				real_out(os,7,2, gain, false);
				os << " ";
				real_out(os,11,4, magnitude);
				os << endl;
				
				phi += _phi_step;
			}
			theta += _theta_step;
		}
	}
};

class nec_antenna_input : public nec_base_result
{
	// Antenna Input Parameters
	vector<int> _tag, _segment;
	vector<nec_float> _power;
	vector<nec_complex> _voltage, _current, _impedance, _admittance;
	long n_items;
	
public:
	nec_antenna_input()
	{
		n_items = 0;
	}

	virtual ~nec_antenna_input()
	{
	}
	
	void set_input(int tag, int segment, nec_complex voltage, nec_complex current, nec_complex impedance, nec_complex admittance, nec_float power)
	{
		n_items++;
		_tag.push_back(tag);
		_segment.push_back(segment);
		_voltage.push_back(voltage);
		_current.push_back(current);
		_impedance.push_back(impedance);
		_admittance.push_back(admittance);
		_power.push_back(power);
	}
	
	virtual void write_to_file(ostream& os)
	{
		if (n_items == 0)
			return;
			
		section_start(os);
		os << "                        --------- ANTENNA INPUT PARAMETERS ---------" << endl;	
		os << "  TAG   SEG       VOLTAGE (VOLTS)         CURRENT (AMPS)         IMPEDANCE (OHMS)        ADMITTANCE (MHOS)     POWER" << endl;
		os << "  NO.   NO.     REAL      IMAGINARY     REAL      IMAGINARY     REAL      IMAGINARY    REAL       IMAGINARY   (WATTS)" << endl;
		for (int i=0; i<n_items; i++)
		{
			int_out(os,4, _tag[i]);
			int_out(os,5, _segment[i]);
			complex_out(os,11,4, _voltage[i]);
			complex_out(os,11,4, _current[i]);
			complex_out(os,11,4, _impedance[i]);
			complex_out(os,11,4, _admittance[i]);
			real_out(os,11,4, _power[i]);
			os << endl;
		}
	}
};

/*
	Stores a nec_result 

	Usage
	
		nec_antenna_input* ai = new nec_antenna_input();
		s_results.add(ai);
		ai->set_intput(tag, segment, voltage, current, impedance, admittance, power);
		ai->set_intput(tag, segment, voltage, current, impedance, admittance, power);
		ai->set_intput(tag, segment, voltage, current, impedance, admittance, power);
	
*/
class nec_results
{
	vector<nec_base_result*> _results;
	int _n;
	bool _file_out;
	
public:
	
	nec_results()
	{
		_n = 0;
		_file_out = false;
	}

	// On destruction we write to a file.
	~nec_results()
	{
		// write_to_file();
		for (int i=0;i<_n;i++)
		{
			delete _results[i];
			_results[i] = NULL;
		}
	}
	
	void add(nec_base_result* br)
	{
		_results.push_back(br);
		_n++;
	}
	
	void write_to_file()
	{
		if (false == _file_out)
			return;
			
		for (int i=0;i<_n;i++)
		{
			if (_results[i]->write_file())
			{
				_results[i]->write_to_file(cout);
				_results[i]->set_write_file(false);
			}
		}
	}
	
	void set_stdout(bool f)
	{
		_file_out = f;
	}
};

#endif /* __nec_results__ */

