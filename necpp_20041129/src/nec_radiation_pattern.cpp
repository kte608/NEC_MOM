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
#include "nec_radiation_pattern.h"
#include "nec_context.h"

#define	CONST3	nec_complex(0.0,-29.97922085)


/*! \brief Write the analyzed data to a file
*/
void nec_radiation_pattern::write_to_file2(ostream& os)
{
	analyze();
	static char  *hpol[4] = { "LINEAR", "RIGHT ", "LEFT  ", " " };
	static char  *igtp[2] = { "----- POWER GAINS ----- ", "--- DIRECTIVE GAINS ---" };
	static char  *igax[4] = { " MAJOR", " MINOR", " VERTC", " HORIZ" };
	
	int i, itmp1, itmp2;
	nec_float exrm=0., exra=0., pint, tmp1, tmp2;
	nec_float phi, pha, thet, tha ;
	
	nec_complex  eth, eph, erd;
		
	if ( _ifar >= 2) 
	{
		section_start(os);
		os << "                                 ------ FAR FIELD GROUND PARAMETERS ------" << endl << endl;
		
		if ( _ifar > 3)
		{
			os << endl;
			os << "                                        RADIAL WIRE GROUND SCREEN" << endl;
			os << "                                        "; int_out(os,5, m_context->ground.radial_wire_count); os << "WIRES" << endl;
			os << "                                        WIRE LENGTH= "; real_out(os,8,2, m_context->ground.radial_wire_length,false); os << "METERS" << endl;
			os << "                                        WIRE RADIUS= "; real_out(os,10,3, m_context->ground.radial_wire_radius); os << "METERS" << endl;
		} /* if ( _ifar > 3) */
		
		if ( _ifar != 4 ) 
		{
			std::string hclif;
			if ( (_ifar == 2) || (_ifar == 5) )
				hclif = "LINEAR";
			if ( (_ifar == 3) || (_ifar == 6) )
				hclif= "CIRCLE";
		
		
			os << endl;
			os << "                                        " << hclif << " CLIFF" << endl;
			os << "                                        EDGE DISTANCE= "; real_out(os,9,2,m_context->ground.cliff_edge_distance,false); os << "METERS" << endl;
			os << "                                        HEIGHT= "; real_out(os,8,2,m_context->ground.cliff_height,false); os << "METERS" << endl;
			os << "                                        SECOND MEDIUM -" << endl;
			os << "                                        RELATIVE DIELECTRIC CONST.= "; real_out(os,7,3,m_context->ground.epsr2, false); os << endl;
			os << "                                        CONDUCTIVITY= "; real_out(os,10,3,m_context->ground.sig2,false); os << "MHOS" << endl;
		} /* if ( m_context->ifar != 4 ) */
		
	} /* if ( _ifar >= 2) */

	if ( _ifar == 1)
	{
		section_start(os);
		os << "                             ------- RADIATED FIELDS NEAR GROUND --------" << endl << endl;
		os << "    ------- LOCATION -------     --- E(THETA) ---     ---- E(PHI) ----    --- E(RADIAL) ---" << endl;
		os << "      RHO    PHI        Z           MAG    PHASE         MAG    PHASE        MAG     PHASE" << endl;
		os << "    METERS DEGREES    METERS      VOLTS/M DEGREES      VOLTS/M DEGREES     VOLTS/M  DEGREES" << endl;
	}
	else
	{
		itmp1 = 2 * m_context->iax;
		itmp2 = itmp1+1;
		
		section_start(os);
		os << "                             ---------- RADIATION PATTERNS -----------" << endl << endl;
		
		if ( m_context->rfld >= 1.0e-20)
		{
			exrm=1./ m_context->rfld;
			exra= m_context->rfld/ _wavelength;
			exra=-360.*( exra- floor( exra));
		
			os << "                             RANGE: "; real_out(os,13,6,m_context->rfld); os << "METERS" << endl;
			os << "                             EXP(-JKR)/R: "; real_out(os,12,5,exrm); os << "AT PHASE: "; real_out(os,7,2,exra,false); os << "DEGREES" << endl;
		}
		
		os << " ---- ANGLES -----     "; string_out(os,23,igtp[m_context->ipd]); os << "      ---- POLARIZATION ----   ---- E(THETA) ----    ----- E(PHI) ------" << endl;
		os << "  THETA      PHI      "; string_out(os,6,igax[itmp1]); os << "  "; string_out(os,6,igax[itmp2]); os << "   TOTAL       AXIAL      TILT  SENSE   MAGNITUDE    PHASE    MAGNITUDE     PHASE" << endl;
		os << " DEGREES   DEGREES        DB       DB       DB       RATIO   DEGREES            VOLTS/M   DEGREES     VOLTS/M   DEGREES" << endl;
		
	} /* if ( _ifar == 1) */


	i=0;
	pint=0.;
	tmp1= m_context->dph* TA;
	tmp2=.5* m_context->dth* TA;
	phi= m_context->phis- m_context->dph;

	for(int kph = 1; kph <= m_context->nph; kph++ )
	{
		phi += m_context->dph;
		pha= phi* TA;
		thet= m_context->thets- m_context->dth;
		
		for(int kth = 1; kth <= m_context->nth; kth++ )
		{
			thet += m_context->dth;
			if ( (m_context->ground.ksymp == 2) && (thet > 90.01) && (_ifar != 1) )
				continue;
		
			tha= thet* TA;
		
			/* elliptical polarization calc. */
			if ( _ifar == 1)
			{
				nec_complex e_theta = _e_theta[i];
				nec_complex e_phi = _e_phi[i];
				nec_complex e_r = _e_r[i];
				
				os << " ";
				real_out(os,9,2,m_context->rfld,false);
				real_out(os,7,2,phi,false); real_out(os,9,2,thet,false);
				real_out(os,11,4,abs(e_theta)); real_out(os,7,2,arg_degrees(e_theta),false); 
				real_out(os,11,4,abs(e_phi)); real_out(os,7,2,arg_degrees(e_phi),false);
				real_out(os,11,4,abs(e_r)); real_out(os,7,2,arg_degrees(e_r),false);
				os << endl;
			}
			else
			{
				nec_complex e_theta = _e_theta[i];
				nec_complex e_phi = _e_phi[i];
				
				char* pol_sense = hpol[_polarization_sense_index[i]];
				
				os << " "; real_out(os,7,2,thet,false); real_out(os,9,2,phi,false);
				os << " "; real_out(os,8,2,_power_gain_vert[i],false); real_out(os,8,2,_power_gain_horiz[i],false); real_out(os,8,2,_power_gain_tot[i],false); real_out(os,11,4,_polarization_axial_ratio[i],false);
				real_out(os,9,2,_polarization_tilt[i],false); string_out(os,6,pol_sense);
				real_out(os,11,4,abs(e_theta)); real_out(os,9,2,arg_degrees(e_theta),false); real_out(os,11,4,abs(e_phi)); real_out(os,9,2,arg_degrees(e_phi),false);
				os << endl;
				
				m_context->plot_card.plot_patterns(thet, phi,
					e_theta, e_phi,
					_power_gain_vert[i], _power_gain_horiz[i], _power_gain_tot[i]);
				
			} /* if ( _ifar != 1) */
		
			i++;
		} /* for( kth = 1; kth <= m_context->nth; kth++ ) */
	
	} /* for( kph = 1; kph <= m_context->nph; kph++ ) */

	if ( m_context->iavp != 0)
	{		
		section_start(os);
		os << "  AVERAGE POWER GAIN: "; real_out(os,11,4,_average_power_gain); os << " - SOLID ANGLE USED IN AVERAGING: ("; real_out(os,7,4,_average_power_solid_angle,false); os << ")*PI STERADIANS" << endl;
	}

	if ( m_context->inor != 0)
		write_normalized_gain(os, m_context->gnor, m_context->inor);
}


/*! \brief Generate the data for the radiation pattern
*/
void nec_radiation_pattern::analyze()
{
	if (m_analysis_done)
		return;
		
	int pol_sense_index;
	int i;
	nec_float exrm=0., exra=0., prad, gcon, gcop, pint, tmp1, tmp2;
	nec_float phi, pha, thet, tha, erdm=0., erda=0., ethm2, ethm;
	nec_float etha, ephm2, ephm, epha, tilta, emajr2, eminr2, pol_axial_ratio;
	nec_float dfaz, dfaz2, cdfaz, tstor1=0., tstor2, stilta, gnmj;
	nec_float gnmn, gnv, gnh, tmp3, tmp4, tmp5, tmp6;
	nec_complex  eth, eph, erd;
		

	if ( (m_context->excitation_type == 0) || (m_context->excitation_type == 5) )
	{
		gcop= _wavelength* _wavelength*2.* pi()/(376.73* _pinr);
		prad= _pinr- m_context->structure_power_loss- _pnlr;
		gcon= gcop;
		if ( m_context->ipd != 0)
			gcon= gcon* _pinr/ prad;
	}
	else 
	if ( m_context->excitation_type == 4)
	{
		_pinr=394.51* m_context->xpr6* m_context->xpr6* _wavelength* _wavelength;
		gcop= _wavelength* _wavelength*2.* pi()/(376.73* _pinr);
		prad= _pinr- m_context->structure_power_loss- _pnlr;
		gcon= gcop;
		if ( m_context->ipd != 0)
			gcon= gcon* _pinr/ prad;
	}
	else
	{
		prad=0.;
		gcon=4.* pi()/(1.0+ m_context->xpr6* m_context->xpr6);
		gcop= gcon;
	}

	i=0;
	pint=0.;
	tmp1= m_context->dph* TA;
	tmp2=.5* m_context->dth* TA;
	phi= m_context->phis- m_context->dph;

	for(int kph = 1; kph <= m_context->nph; kph++ )
	{
		phi += m_context->dph;
		pha= phi* TA;
		thet= m_context->thets- m_context->dth;
		
		for(int kth = 1; kth <= m_context->nth; kth++ )
		{
			thet += m_context->dth;
			if ( (m_context->ground.ksymp == 2) && (thet > 90.01) && (_ifar != 1) )
				continue;
		
			tha= thet* TA;
			if ( _ifar != 1)
				ffld(m_context, tha, pha, &eth, &eph, _wavelength);
			else
			{
				gfld(m_context, m_context->rfld/_wavelength, pha, thet/_wavelength,
				&eth, &eph, &erd, m_context->zrati, m_context->ground.ksymp, _wavelength );
				erdm= abs( erd);
				erda= arg_degrees( erd);
			}
		
			ethm2= norm(eth);
			ethm= sqrt(ethm2);
			etha= arg_degrees(eth);
			ephm2= norm(eph);
			ephm= sqrt( ephm2);
			epha= arg_degrees( eph);
		
			if ( _ifar == 1)
			{
				nec_complex e_theta = deg_polar(ethm, etha);
				nec_complex e_phi = deg_polar(ephm, epha);
				nec_complex e_r = deg_polar(erdm, erda);
				
				
				_e_theta[i] = e_theta;
				_e_phi[i] = e_phi;
				_e_r[i] = e_r;
			}
			else // _ifar != 1
			{ 
				/* elliptical polarization calc. */
				if ( (ethm2 <= 1.0e-20) && (ephm2 <= 1.0e-20) )
				{
					tilta=0.;
					emajr2=0.;
					eminr2=0.;
					pol_axial_ratio=0.;
					pol_sense_index = 3;
				}
				else
				{
					dfaz= epha- etha;
					if ( epha >= 0.)
						dfaz2= dfaz-360.;
					else
						dfaz2= dfaz+360.;
				
					if ( fabs(dfaz) > fabs(dfaz2) )
						dfaz= dfaz2;
				
					cdfaz= cos( dfaz* TA);
					tstor1= ethm2- ephm2;
					tstor2=2.* ephm* ethm* cdfaz;
					tilta=.5* atan2( tstor2, tstor1);
					stilta= sin( tilta);
					tstor1= tstor1* stilta* stilta;
					tstor2= tstor2* stilta* cos( tilta);
					emajr2=- tstor1+ tstor2+ ethm2;
					eminr2= tstor1- tstor2+ ephm2;
					if ( eminr2 < 0.)
						eminr2=0.;
				
					pol_axial_ratio= sqrt( eminr2/ emajr2);
					tilta= rad_to_degrees(tilta);
					if ( pol_axial_ratio <= 1.0e-5)
						pol_sense_index = 0;
					else
					if ( dfaz <= 0.)
						pol_sense_index = 1;
					else
						pol_sense_index = 2;
				
				} /* if ( (ethm2 <= 1.0e-20) && (ephm2 <= 1.0e-20) ) */
			
				gnmj= db10( gcon* emajr2);
				gnmn= db10( gcon* eminr2);
				gnv = db10( gcon* ethm2);
				gnh = db10( gcon* ephm2);
				nec_float gtot= db10( gcon*(ethm2+ ephm2) );
			
				if ( m_context->inor > 0)
				{
					switch( m_context->inor )
					{
					case 1:
						tstor1= gnmj;
						break;
				
					case 2:
						tstor1= gnmn;
						break;
				
					case 3:
						tstor1= gnv;
						break;
				
					case 4:
						tstor1= gnh;
						break;
				
					case 5:
						tstor1= gtot;
					}
				
					_gain[i]= tstor1;
				
				} /* if ( m_context->inor > 0) */
			
				if ( m_context->iavp != 0)
				{
					tstor1= gcop*( ethm2+ ephm2);
					tmp3= tha- tmp2;
					tmp4= tha+ tmp2;
				
					if ( kth == 1)
						tmp3= tha;
					else
					if ( kth == m_context->nth)
						tmp4= tha;
				
					nec_float da = fabs( tmp1*( cos( tmp3)- cos( tmp4)));
					if ( (kph == 1) || (kph == m_context->nph) )
						da *=.5;
					pint += tstor1 * da;
				
					if ( m_context->iavp == 2)
						continue;
				}
			
				if ( m_context->iax != 1)
				{
					tmp5= gnmj;
					tmp6= gnmn;
				}
				else
				{
					tmp5= gnv;
					tmp6= gnh;
				}
			
				ethm= ethm* _wavelength;
				ephm= ephm* _wavelength;
			
				if ( m_context->rfld >= 1.0e-20 )
				{
					ethm= ethm* exrm;
					etha= etha+ exra;
					ephm= ephm* exrm;
					epha= epha+ exra;
				}
			
				nec_complex e_theta = deg_polar(ethm, etha);
				nec_complex e_phi = deg_polar(ephm, epha);
				
				_power_gain_vert[i] = tmp5;
				_power_gain_horiz[i] = tmp6;
				_power_gain_tot[i] = gtot;
				_polarization_axial_ratio[i] = pol_axial_ratio;
				_polarization_tilt[i] = tilta;
				_polarization_sense_index[i] = pol_sense_index;
				
				_e_theta[i] = e_theta;
				_e_phi[i] = e_phi;
			} /* if ( _ifar == 1) */
			i++;
		} /* for( kth = 1; kth <= m_context->nth; kth++ ) */
	
	} /* for( kph = 1; kph <= m_context->nph; kph++ ) */

	if ( m_context->iavp != 0)
	{
		tmp3= m_context->thets* TA;
		tmp4= tmp3+ m_context->dth* TA* (nec_float)( m_context->nth-1);
		tmp3= fabs( m_context->dph* TA* (nec_float)( m_context->nph-1)*( cos( tmp3)- cos( tmp4)));
		pint /= tmp3;
		tmp3 /= pi();
		
		_average_power_gain = pint;
		_average_power_solid_angle = tmp3;
	}

	_maximum_gain = _gain.max();
	m_analysis_done = true;
}


nec_float nec_radiation_pattern::get_gain_normalization_factor(nec_float gnor)
{
	if ( fabs(gnor) > 1.0e-20)
		return gnor;
		
	analyze();
	return _maximum_gain;
}

void nec_radiation_pattern::write_normalized_gain(ostream& os, nec_float gnor, int inor)
{
	// if  is non-zero then use it as a normaliztion factor
	
	nec_float normalization_factor = get_gain_normalization_factor(gnor);
		
	string norm_type;
	switch (inor)
	{
		case 1:
			norm_type = "  MAJOR AXIS"; break;
		case 2:
			norm_type = "  MINOR AXIS"; break;
		case 3:
			norm_type = "    VERTICAL"; break;
		case 4:
			norm_type = "  HORIZONTAL"; break;
		case 5:
			norm_type = "      TOTAL "; break;
		
		default: throw 100; // eek
	}
	
	section_start(os);
	os << "                              ---------- NORMALIZED GAIN ----------" << endl;
	os << "                                      " << norm_type << " GAIN" << endl;
	os << "                                   NORMALIZATION FACTOR: "; real_out(os,7,2,normalization_factor,false); os << "db" << endl << endl;
	os << "    ---- ANGLES ----                ---- ANGLES ----                ---- ANGLES ----" << endl;
	os << "    THETA      PHI        GAIN      THETA      PHI        GAIN      THETA      PHI       GAIN" << endl;
	os << "   DEGREES   DEGREES        DB     DEGREES   DEGREES        DB     DEGREES   DEGREES       DB" << endl;

	
	int row_count = 0;
	int n_cols = 3;
	
	int item_count = 0;
	for (int p=0;p<m_context->nph;p++)
	{
		nec_float phi = p*m_context->dph ;
		for (int t=0;t<m_context->nth;t++)
		{
			nec_float theta = t * m_context->dth;
			nec_float norm_gain = _gain[item_count++] - normalization_factor;
			
			os << " "; real_out(os,9,2,theta,false); real_out(os,9,2,phi,false); real_out(os,9,2,norm_gain,false);
		
			if (item_count % n_cols == 0)
			{
				row_count++;
				os << endl;
			}
		}
	}
	os << endl;
}

/*! \brief gfld computes the radiated field including ground wave.
*/
void gfld(nec_context* m_context, nec_float rho, nec_float phi, nec_float rz,
    nec_complex *eth, nec_complex *epi,
    nec_complex *erd, nec_complex ux, int ksymp, nec_float _wavelength )
{
	int i, k;
	nec_float b, r, thet, arg, phx, phy, rx, ry, dx, dy, dz, rix, riy, rhs, rhp;
	nec_float rhx, rhy, calp, cbet, sbet, cph, sph, el, rfl, riz, thx, thy, thz;
	nec_float rxyz, rnx, rny, rnz, omega, sill, top, bot, a, too, boo, c, rr, ri;
	nec_complex cix, ciy, ciz, exa, erv;
	nec_complex ezv, erh, eph, ezh, ex, ey;
	
	r= sqrt( rho*rho+ rz*rz );
	if ( (ksymp == 1) || (abs(ux) > .5) || (r > 1.e5) )
	{
		/* computation of space wave only */
		if ( rz >= 1.0e-20)
		thet= atan( rho/ rz);
		else
		thet= pi()*.5;
	
		ffld(m_context, thet, phi, eth, epi, _wavelength);
		arg= -two_pi() * r;
		exa= nec_complex( cos( arg), sin( arg))/ r;
		*eth= *eth* exa;
		*epi= *epi* exa;
		*erd=cplx_00();
		return;
	} /* if ( (ksymp == 1) && (abs(ux) > .5) && (r > 1.e5) ) */
	
	/* computation of space and m_context->ground waves. */
	m_context->ground_wave.set_u(ux);
	phx=- sin( phi);
	phy= cos( phi);
	rx= rho* phy;
	ry=- rho* phx;
	cix=cplx_00();
	ciy=cplx_00();
	ciz=cplx_00();
	
	/* summation of field from individual segments */
	for( i = 0; i < m_context->geometry.n; i++ )
	{
		dx= m_context->geometry.cab[i];
		dy= m_context->geometry.sab[i];
		dz= m_context->geometry.salp[i];
		rix= rx- m_context->geometry.x[i];
		riy= ry- m_context->geometry.y[i];
		rhs= rix* rix+ riy* riy;
		rhp= sqrt( rhs);
	
		if ( rhp >= 1.0e-6)
		{
			rhx= rix/ rhp;
			rhy= riy/ rhp;
		}
		else
		{
			rhx=1.;
			rhy=0.;
		}
	
		calp=1.- dz* dz;
		if ( calp >= 1.0e-6)
		{
			calp= sqrt( calp);
			cbet= dx/ calp;
			sbet= dy/ calp;
			cph= rhx* cbet+ rhy* sbet;
			sph= rhy* cbet- rhx* sbet;
		}
		else
		{
			cph= rhx;
			sph= rhy;
		}
	
		el= pi()* m_context->geometry.segment_length[i];
		rfl=-1.;
	
		/* Integration of (current)*(phase factor) over segment and image for 
		   constant, sine, and cosine current distributions */
		for( k = 0; k < 2; k++ )
		{
			rfl=- rfl;
			riz= rz- m_context->geometry.z[i]* rfl;
			rxyz= sqrt( rix* rix+ riy* riy+ riz* riz);
			rnx= rix/ rxyz;
			rny= riy/ rxyz;
			rnz= riz/ rxyz;
			omega=-( rnx* dx+ rny* dy+ rnz* dz* rfl);
			sill= omega* el;
			top= el+ sill;
			bot= el- sill;
		
			if ( fabs( omega) >= 1.0e-7)
				a=2.* sin( sill)/ omega;
			else
				a=(2.- omega* omega* el* el/3.)* el;
		
			if ( fabs( top) >= 1.0e-7)
				too= sin( top)/ top;
			else
				too=1.- top* top/6.;
		
			if ( fabs( bot) >= 1.0e-7)
				boo= sin( bot)/ bot;
			else
				boo=1.- bot* bot/6.;
		
			b= el*( boo- too);
			c= el*( boo+ too);
			rr= a* m_context->air[i]+ b* m_context->bii[i]+ c* m_context->cir[i];
			ri= a* m_context->aii[i]- b* m_context->bir[i]+ c* m_context->cii[i];
			arg= two_pi()*( m_context->geometry.x[i]* rnx+ m_context->geometry.y[i]* rny+ m_context->geometry.z[i]* rnz* rfl);
			exa= nec_complex( cos( arg), sin( arg))* nec_complex( rr, ri)/two_pi();
		
			if ( k != 1 )
			{
				m_context->ground_wave.xx1= exa;
				m_context->ground_wave.r1= rxyz;
				m_context->ground_wave.zmh= riz;
				continue;
			}
		
			m_context->ground_wave.xx2 = exa;
			m_context->ground_wave.r2= rxyz;
			m_context->ground_wave.zph= riz;
		
		} /* for( k = 0; k < 2; k++ ) */
	
		/* call subroutine to compute the field */
		/* of segment including m_context->ground wave. */
		gwave( &erv, &ezv, &erh, &ezh, &eph, m_context->ground_wave);
		erh= erh* cph* calp+ erv* dz;
		eph= eph* sph* calp;
		ezh= ezh* cph* calp+ ezv* dz;
		ex= erh* rhx- eph* rhy;
		ey= erh* rhy+ eph* rhx;
		cix= cix+ ex;
		ciy= ciy+ ey;
		ciz= ciz+ ezh;
	
	} /* for( i = 0; i < n; i++ ) */
	
	arg= -two_pi() * r;
	exa= nec_complex( cos( arg), sin( arg));
	cix= cix* exa;
	ciy= ciy* exa;
	ciz= ciz* exa;
	rnx= rx/ r;
	rny= ry/ r;
	rnz= rz/ r;
	thx= rnz* phy;
	thy=- rnz* phx;
	thz=- rho/ r;
	*eth= cix* thx+ ciy* thy+ ciz* thz;
	*epi= cix* phx+ ciy* phy;
	*erd= cix* rnx+ ciy* rny+ ciz* rnz;
}

/* ffld calculates the far zone radiated electric fields, */
/* the factor exp(j*k*r)/(r/lamda) not included */
void ffld(nec_context* m_context, nec_float thet, nec_float phi,
    nec_complex *eth, nec_complex *eph, nec_float _wavelength )
{
  int k, i, ip;
  bool jump;
  nec_float phx, phy, roz, rozs, thx, thy, thz, rox, roy;
  nec_float tthet=0., darg=0., omega, el, sill, top, bot, a;
  nec_float too, boo, b, c, d, rr, ri, arg, dr, rfl, rrz;
  nec_complex cix, ciy, ciz, exa, ccx, ccy, ccz, cdp;
  nec_complex zrsin, rrv, rrh, rrv1, rrh1, rrv2, rrh2;
  nec_complex tix, tiy, tiz, zscrn, ex, ey, ez, gx, gy, gz;

  phx=- sin( phi);
  phy= cos( phi);
  roz= cos( thet);
  rozs= roz;
  thx= roz* phy;
  thy=- roz* phx;
  thz=- sin( thet);
  rox=- thz* phy;
  roy= thz* phx;

  jump = false;
  if ( m_context->geometry.n != 0)
  {
    /* loop for structure image if any */
    /* calculation of reflection coeffecients */
    for( k = 0; k < m_context->ground.ksymp; k++ )
    {
      if ( k != 0 )
      {
	/* for perfect m_context->ground */
	if (m_context->ground.type_perfect()) // (  m_context->ground.iperf == 1)
	{
	  rrv=-cplx_10();
	  rrh=-cplx_10();
	}
	else
	{
	  /* for infinite planar m_context->ground */
	  zrsin= sqrt(1.- m_context->zrati* m_context->zrati* thz* thz);
	  rrv=-( roz- m_context->zrati* zrsin)/( roz+ m_context->zrati* zrsin);
	  rrh=( m_context->zrati* roz- zrsin)/( m_context->zrati* roz+ zrsin);

	} /* if (  m_context->ground.iperf == 1) */

	/* for the cliff problem, two reflction coefficients calculated */
	if ( m_context->ifar > 1)
	{
	  rrv1= rrv;
	  rrh1= rrh;
	  tthet= tan( thet);

	  if ( m_context->ifar != 4)
	  {
	    nec_complex zrati2 = m_context->ground.get_zrati2(_wavelength);
		
		zrsin = sqrt(1.-  zrati2 *  zrati2 * thz* thz);
	    rrv2 =-( roz-  zrati2* zrsin)/( roz+  zrati2* zrsin);
	    rrh2 =(  zrati2* roz- zrsin)/(  zrati2* roz+ zrsin);
	    darg = -two_pi() * 2.0 * m_context->ground.get_ch(_wavelength) * roz;
	  }
	} /* if ( m_context->ifar > 1) */

	roz=- roz;
	ccx= cix;
	ccy= ciy;
	ccz= ciz;

      } /* if ( k != 0 ) */

      cix=cplx_00();
      ciy=cplx_00();
      ciz=cplx_00();

      /* loop over structure segments */
      for( i = 0; i < m_context->geometry.n; i++ )
      {
	omega=-( rox* m_context->geometry.cab[i]+ roy* m_context->geometry.sab[i]+ roz* m_context->geometry.salp[i]);
	el= pi()* m_context->geometry.segment_length[i];
	sill= omega* el;
	top= el+ sill;
	bot= el- sill;

	if ( fabs( omega) >= 1.0e-7)
	  a=2.* sin( sill)/ omega;
	else
	  a=(2.- omega* omega* el* el/3.)* el;

	if ( fabs( top) >= 1.0e-7)
	  too= sin( top)/ top;
	else
	  too=1.- top* top/6.;

	if ( fabs( bot) >= 1.0e-7)
	  boo= sin( bot)/ bot;
	else
	  boo=1.- bot* bot/6.;

	b= el*( boo- too);
	c= el*( boo+ too);
	rr= a* m_context->air[i]+ b* m_context->bii[i]+ c* m_context->cir[i];
	ri= a* m_context->aii[i]- b* m_context->bir[i]+ c* m_context->cii[i];
	arg= two_pi()*( m_context->geometry.x[i]* rox+ m_context->geometry.y[i]* roy+ m_context->geometry.z[i]* roz);

	if ( (k != 1) || (m_context->ifar < 2) )
	{
	  /* summation for far field integral */
	  exa= nec_complex( cos( arg), sin( arg))* nec_complex( rr, ri);
	  cix= cix+ exa* m_context->geometry.cab[i];
	  ciy= ciy+ exa* m_context->geometry.sab[i];
	  ciz= ciz+ exa* m_context->geometry.salp[i];
	  continue;
	}

	/* calculation of image contribution */
	/* in cliff and m_context->ground screen problems */

	/* specular point distance */
	dr= m_context->geometry.z[i]* tthet;

	d= dr* phy+ m_context->geometry.x[i];
	if ( m_context->ifar == 2)
	{
	  if (( m_context->ground.get_cl(_wavelength) - d) > 0.0)
	  {
	    rrv= rrv1;
	    rrh= rrh1;
	  }
	  else
	  {
	    rrv= rrv2;
	    rrh= rrh2;
	    arg= arg+ darg;
	  }
	} /* if ( m_context->ifar == 2) */
	else
	{
	  d= sqrt( d*d + (m_context->geometry.y[i]-dr*phx)*(m_context->geometry.y[i]-dr*phx) );
	  if ( m_context->ifar == 3)
	  {
	    if (( m_context->ground.get_cl(_wavelength) - d) > 0.0)
	    {
	      rrv= rrv1;
	      rrh= rrh1;
	    }
	    else
	    {
	      rrv= rrv2;
	      rrh= rrh2;
	      arg= arg+ darg;
	    }
	  } /* if ( m_context->ifar == 3) */
	  else
	  {
	    if (( m_context->scrwl- d) >= 0.)
	    {
	      /* radial wire m_context->ground screen reflection coefficient */
	      d= d+ m_context->t2;
	      zscrn= m_context->t1* d* log( d/ m_context->t2);
	      zscrn=( zscrn* m_context->zrati)/( ETA* m_context->zrati+ zscrn);
	      zrsin= sqrt(1.- zscrn* zscrn* thz* thz);
	      rrv=( roz+ zscrn* zrsin)/(- roz+ zscrn* zrsin);
	      rrh=( zscrn* roz+ zrsin)/( zscrn* roz- zrsin);
	    } /* if (( m_context->scrwl- d) < 0.) */
	    else
	    {
	      if ( m_context->ifar == 4)
	      {
		rrv= rrv1;
		rrh= rrh1;
	      } /* if ( m_context->ifar == 4) */
	      else
	      {
		if ( m_context->ifar == 5)
		  d= dr* phy+ m_context->geometry.x[i];

		if (( m_context->ground.get_cl(_wavelength) - d) > 0.)
		{
		  rrv= rrv1;
		  rrh= rrh1;
		}
		else
		{
		  rrv= rrv2;
		  rrh= rrh2;
		  arg= arg+ darg;
		} /* if (( cl- d) > 0.) */

	      } /* if ( m_context->ifar == 4) */

	    } /* if (( m_context->scrwl- d) < 0.) */

	  } /* if ( m_context->ifar == 3) */

	} /* if ( m_context->ifar == 2) */

	/* contribution of each image segment modified by */
	/* reflection coef, for cliff and m_context->ground screen problems */
	exa= nec_complex( cos( arg), sin( arg))* nec_complex( rr, ri);
	tix= exa* m_context->geometry.cab[i];
	tiy= exa* m_context->geometry.sab[i];
	tiz= exa* m_context->geometry.salp[i];
	cdp=( tix* phx+ tiy* phy)*( rrh- rrv);
	cix= cix+ tix* rrv+ cdp* phx;
	ciy= ciy+ tiy* rrv+ cdp* phy;
	ciz= ciz- tiz* rrv;

      } /* for( i = 0; i < n; i++ ) */

      if ( k == 0 )
	continue;

      /* calculation of contribution of structure image for infinite m_context->ground */
      if ( m_context->ifar < 2)
      {
	cdp=( cix* phx+ ciy* phy)*( rrh- rrv);
	cix= ccx+ cix* rrv+ cdp* phx;
	ciy= ccy+ ciy* rrv+ cdp* phy;
	ciz= ccz- ciz* rrv;
      }
      else
      {
	cix= cix+ ccx;
	ciy= ciy+ ccy;
	ciz= ciz+ ccz;
      }

    } /* for( k=0; k < m_context->ground.ksymp; k++ ) */

    if ( m_context->geometry.m > 0)
      jump = true;
    else
    {
      *eth=( cix* thx+ ciy* thy+ ciz* thz)* CONST3;
      *eph=( cix* phx+ ciy* phy)* CONST3;
      return;
    }

  } /* if ( n != 0) */

  if ( ! jump )
  {
    cix=cplx_00();
    ciy=cplx_00();
    ciz=cplx_00();
  }

  /* electric field components */
  roz= rozs;
  rfl=-1.;
  for( ip = 0; ip < m_context->ground.ksymp; ip++ )
  {
    rfl=- rfl;
    rrz= roz* rfl;
    complex_array temp = m_context->current_vector.sub_array(m_context->geometry.n);
    m_context->geometry.fflds(rox, roy, rrz, temp, &gx, &gy, &gz);

    if ( ip != 1 )
    {
      ex= gx;
      ey= gy;
      ez= gz;
      continue;
    }

    if (m_context->ground.type_perfect()) // (  m_context->ground.iperf == 1)
    {
      gx=- gx;
      gy=- gy;
      gz=- gz;
    }
    else
    {
      rrv= sqrt(1.- m_context->zrati* m_context->zrati* thz* thz);
      rrh= m_context->zrati* roz;
      rrh=( rrh- rrv)/( rrh+ rrv);
      rrv= m_context->zrati* rrv;
      rrv=-( roz- rrv)/( roz+ rrv);
      *eth=( gx* phx+ gy* phy)*( rrh- rrv);
      gx= gx* rrv+ *eth* phx;
      gy= gy* rrv+ *eth* phy;
      gz= gz* rrv;

    } /* if (  m_context->ground.iperf == 1) */

    ex= ex+ gx;
    ey= ey+ gy;
    ez= ez- gz;

  } /* for( ip = 0; ip < m_context->ground.ksymp; ip++ ) */

  ex= ex+ cix* CONST3;
  ey= ey+ ciy* CONST3;
  ez= ez+ ciz* CONST3;
  *eth= ex* thx+ ey* thy+ ez* thz;
  *eph= ex* phx+ ey* phy;
}

