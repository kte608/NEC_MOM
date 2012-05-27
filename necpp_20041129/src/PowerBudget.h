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
#ifndef __PowerBudget__
#define __PowerBudget__

#include <vector>
#include "BaseInput.h"

/*
			---------- POWER BUDGET ---------
			INPUT POWER   =  3.3203E-05 Watts
			RADIATED POWER=  3.3203E-05 Watts
			STRUCTURE LOSS=  0.0000E+00 Watts
			NETWORK LOSS  =  0.0000E+00 Watts
			EFFICIENCY    =  100.00 Percent
*/
class PowerBudget : public BaseInput
{
public:	
	vector<double> input_power, radiated_power, structure_loss, network_loss, efficiency;
	long n_items;
	
	double get_power_line(istream& if1)
	{
		string line = readline(if1);
		// delete up to '=' character (+1 includes the = character)
		line.erase(0, line.find("=",0) + 1);
		stringstream ss(line);
		double ret = read_fixed(ss);
		return ret;
	}
	
	PowerBudget(istream& if1)
	{
		n_items = 0;
		string searchString("POWER BUDGET");
		while (if1.good())
		{
		string line = readline(if1);

		if (line.find(searchString,0) != string::npos)
		{
			{
				double p = get_power_line(if1);
				input_power.push_back(p);
			}
			{
				double p = get_power_line(if1);
				radiated_power.push_back(p);
			}
			{
				double p = get_power_line(if1);
				structure_loss.push_back(p);
			}
			{
				double p = get_power_line(if1);
				network_loss.push_back(p);
			}
			{
				double p = get_power_line(if1);
				efficiency.push_back(p);
			}
			n_items++;
		}	
		}
	}

	bool equalto(const PowerBudget& pb)
	{
		if (difference(pb) > 1e-4)
			return false;

		return true;
	}

	double difference(const PowerBudget& pb)
	{
		double ret = 0;

		if (n_items != pb.n_items)
			return 1.0;

		try
		{
		for (long i=0; i<n_items; i++)
		{
			ret += diff(pb.input_power[i],input_power[i]);
			ret += diff(pb.radiated_power[i],radiated_power[i]);
			ret += diff(pb.structure_loss[i],structure_loss[i]);
			ret += diff(pb.network_loss[i],network_loss[i]);
			ret += diff(pb.efficiency[i],efficiency[i]);
		}
		}
		catch(string message)
		{
			cout << "diff : " << message << endl;
		}
		return ret;
	};
};

#endif /* __PowerBudget__ */
