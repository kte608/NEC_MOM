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

#include "common.h"
#include "c_ggrid.h"
#include "c_geometry.h"
#include "math_util.h"
#include "matrix_algebra.h"
#include "electromag.h"
#include "nec_radiation_pattern.h"
#include "nec_results.h"
#include "nec_output.h"
#include "nec_ground.h"
#include "c_plot_card.h"


enum excitation_return
{
	FREQ_PRINT_NORMALIZATION = 0,
	FREQ_LOOP_CONTINUE = 1,
	FREQ_LOOP_CARD_CONTINUE = 2
};

/*!
 *  A nec_context class. A more elaborate class description.
 */
class nec_context
{
public:
	nec_context();
	
	// Called after construction...
	void initialize();
	

	void calc_prepare();
	
	c_geometry& get_geometry()
	{
		return geometry;
	}
	
	/*! \brief Get the maximum gain in dB.
	
	This function requires a previous rp_card() method to have been called.
	
	\return The maximum gain in dB or -999.0 if no radiation pattern had been previously requested. This only works
	for a SINGLE FRQUENCY rp_card request.
	*/
	double get_maximum_gain()
	{
		return m_maximum_gain_db;
	}
	
	/*! \brief Get the Radiation Pattern
	
	This function requires a previous rp_card() method to have been called. It is designed for 
	efficient plotting of NEC results.
	
	\return The requested radiation pattern data. This only works
	for a SINGLE FRQUENCY rp_card() request.
	
	*/
	double get_radiation_pattern(int theta_index, int phi_index)
	{
		return m_radiation_pattern[theta_index*nth + phi_index];
	}
	
	void set_output(nec_output_file in_output, nec_output_flags in_output_flags)
	{
		m_output = in_output;
		m_output_flags = in_output_flags;
		
		m_output_fp = m_output.get_fp();
	}
	
	void set_results_stdout(bool flag)
	{
		m_results.set_stdout(flag);
	}
	
	void set_gain_only(bool flag)
	{
		m_output_flags.set_gain_only(flag);
	}
	
	/*! \brief Signal the end of a geometry description.
	
	This function prepares for a calculation by calling calc_prepare().
	*/
	void geometry_complete(int card_int_1, int card_int_2)
	{
		DEBUG_TRACE("geometry_complete()");
		geometry.geometry_complete(this, card_int_1, card_int_2);
		calc_prepare();
	}

	
	
	
	/*! wire card,
	
	All co-ordinates are in meters.
	
		\param tag_id The tag ID.
		\param segment_count The number of segments.
		\param xw1 The x coordinate of the wire starting point.
		\param yw1 The y coordinate of the wire starting point.
		\param zw1 The z coordinate of the wire starting point.
		\param xw2 The x coordinate of the wire ending point.
		\param yw2 The y coordinate of the wire ending point.
		\param zw2 The z coordinate of the wire ending point.
		\param rad The wire radius (meters)
		\param rdel For tapered wires, the. Otherwise set to 1.0
		\param rrad For tapered wires, the. Otherwise set to 1.0
	*/
	void wire(int tag_id, int segment_count,
		nec_float xw1, nec_float yw1, nec_float zw1,
		nec_float xw2, nec_float yw2, nec_float zw2,
		nec_float rad, nec_float rdel, nec_float rrad)
	{
		geometry.wire(tag_id, segment_count, xw1, yw1, zw1, xw2, yw2, zw2, rad, rdel, rrad);
	}
	
	/*! "fr" card, frequency parameters
	\verbatim	
	FREQUENCY
	I1- O= LINEAR STEP, 1=MULTIPLICATIVE
	I2- NO. STEPS, BLANK=1
	F1- FREQUENCY OR START FREQUENCY
	F2- FREQ INCREMENT, ADD OR MULTIPLY
	\endverbatim
	*/
	void fr_card(int in_ifrq, int in_nfrq, nec_float in_freq_mhz, nec_float in_del_freq);
	
	/*! 1: "ld" card, loading parameters
	\verbatim
	LD	LOADING
		itmp1- 	-1 CANCEL LOADS,
			0=SERIES RLC LUMP,
			1=PARALLEL RLC LUMP,
			2=SERIES DIST.,
			3=PARALLEL DIST. (A),
			4=Z (B),
			5=WIRE COND. (C)
		itmp2- TAG# TO BE LOADED, BLANK/0= USE ABSOLUTE #s	
		itmp3- SEG# OF TAG # TO START LOADS, OR ABSOLUTE SEG#
		itmp4- SEG# OF TAG# TO END LOADS, OR OR ABSOLUTE SEG#
		F1- RES., OHMS, OR (A) OHMS/UNIT LENGTH, OR (B) RES. OR (C) OHMS/METER
		F2- IND., HENRY, OR (A) HY/LENGTH OR (B) REACT. OR (C) BLANK
		F3- CAP,. FARAD, OR (A,B) BLANK
	\endverbatim
	*/
	void ld_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3);
	
	/*! \brief Ground parameters under the antenna
	
	\remark Specifies the relative dielectric constant and conductivity of ground in the vicinity of the antenna. In addition, a second set of ground parameters for a second medium can be specified, or a radial wire ground screen can be modeled using a reflection coefficient approximation.
	
	\param ground_type (was IPERF) Ground-type flag. The options are:
		\arg \c -1 - Nullifies ground parameters previously used and sets free-space condition. The remainder of the parameters should be zero in this case. 
		\arg \c O - Finite ground, reflection-coefficient approximation.
		\arg \c 1 - Perfectly conducting ground.
		\arg \c 2 - Finite ground, Sommerfeld/Norton method. 
	
	\param rad_wire_count (was NRADL) - Number of radial wires in the ground screen approximation; Set to zero implies no ground screen.
	
	\param EPSE (F1) - Relative dielectric constant for ground in the vicinity of the antenna. Set to zero in case of a perfect ground.
	\param SIG (F2) - Conductivity in mhos/meter of the ground in the vicinity of the antenna. Set to zero in the case of a perfect ground. If SIG is input as a negative number, the complex dielectric constant Ec = Er -j*sigma/(omega*epsilonzero) is set to EPSR - |SIG|.
	
	\remark
	 
	Options for Remaining Floating Point Fields (F3-F6):
		\li a. For an infinite ground plane, F3 through F6 are blank. 
		\li b. Radial wire ground screen approximation (NRADL nonzero). The ground screen is always centered at the origin, i.e., at (0,0,0), and lies in the XY plane. (F3) - The radius of the screen in meters. (F4) - Radius of the wires used in the screen, in meters. (F5) & (F6) - Blank.
		\li c. Second medium parameters (NRADL = O) for medium outside the region of the first medium (cliff problem). These parameters alter the far field patterns but do not affect the antenna impedance or current distribution. (F3) - Relative dielectric constant of medium 2. (F4) - Conductivity of medium 2 in mhos/meter. (F5) - Distance in meters from the origin of the coordinate system to the join between medium 1 and 2. This distance is either the radius of the circle where the two media join or the distance out the positive X axis to where the two media join in a line parallel to the Y axis. Specification of the circular or linear option is on the RP card. See Figure 16. (F6) - Distance in meters (positive or zero) by which the surface of medium 2 is below medium 1.
	*/
	void gn_card(int ground_type, int rad_wire_count, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);
	
	
	/*! "ex" card, excitation parameters
	\verbatim
			EX	EXCITE STRUCTURE, LAST ENCOUNTERED=USED
				I1- 0=E VOLTAGE (A), 1=LINEAR WAVE (B), 2= R CIRC WAVE (B)
				3=L CIRC WAVE (B), 4= CURRENT (C), 5= VOLTAGE DISC. (A)
				I2- (A) SOURCE TAG#, (B) # TH ANGLS, (C) BLANK
				I3- (A) SOURCE SEG#, (B) # PH ANGLS, (C) BLANK
				I4- (A) XX= ADMIT.,IMPED. PRINT, X=0 NO/1 DO, (BC), 1= ADM. PRINT
				F1- (A) EREAL, (B) TH ANGL, (C) X OF SOURCE
				F2- (A) EIMAG, (B) PH ANGL, (C) Y OF SOURCE
				F3- (A) NORM FOR I4, (B) ET ANGL, Z OF SOURCE
				F4- (A) BLANK, (B) TH INC, (C) ALPHA ANGLE FROM XY
				F5- (A) BLANK, (B) PH INC, (C) BETA ANGLE FROM X
				F6- (A) BLANK, (B) MIN/MAJ AXIS, PRODUCT AMPS X LENGTH
				
				// NOT YET DONE... F7- (A) BLANK, (B) INCIDENT AMPLITUDE (Volts/m)
	\endverbatim
	*/
	void ex_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);
	
	
	
	/*! 5: "tl" card, transmission line parameters
	
	\remark
	To generate a transmission line between any two points on the structure. Characteristic impedance, length, and shunt admittance are the defining parameters.
	
	\verbatim
	TL TRANSMISSION LINE 
		I1- PORT 1 TAG #, BLANK/0, USE I2 AS ABSOLUTE
		I2- SEGMENT#, OR ABSOLUTE END 1 SEGMENT, -1=CANCEL NETS/LINES
		I3- AS I1 FOR PORT 2
		I4- AS I2 FOR PORT 2
		F1- LINE Zo, -=CROSSED LINE
		F2- LINE LENGTH METERS, BLANK=STRAIGHT LINE P1 TO P2
		F3- REAL SHUNT ADM., END 1 MHOS
		F4- IMAG SHUNT ADM., END 1
		F5- REAL SHUNT ADM., END 2
		F6- IMAG SHUNT ADM., END 2
	\endverbatim
	*/
	void tl_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);
	
	/*! 4: "nt" card, network parameters
	\verbatim
		NT	NETWORKS
			I1- PORT 1 TAG #, BLANK/0, USE I2 AS ABSOLUTE
			I2- SEGMENT#, OR ABSOLUTE END 1 SEGMENT, -1=CANCEL NETS/LINES
			I3- AS I1 FOR PORT 2
			I4- AS I2 FOR PORT 2
			F1- REAL OF Y(11), MHOS
			F2- IMAG OF Y(11)
			F3- REAL OF Y(12)
			F4- IMAG OF Y(12)
			F5- REAL OF Y(22)
			F6- IMAG OF Y(22)
	\endverbatim
	*/
	void nt_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);
	
	/*! "xq" execute card - calc. including radiated fields
	
	\verbatim
	XQ	EXECUTE ACCUMULATED CARD DECK
		itmp1-
			0=NO PATTERN,
			1=XY PATTERN,
			2= YZ PATTERN,
			3=BOTH
		(DO NOT USE FOR RADIAL GND SCREEN OR 2ND GND MEDIUM)
	
		NOTES: FOR A SINGLE FREQUENCY, XQ, NE, NH, RP CAUSE IMMEDIATE EXECUTION
		FOR MULTIPLE FREQS, ONLY XQ, RP CAUSE EXECUTION
	\endverbatim
	*/
	void xq_card(int itmp1);
	
	/*! "gd" card, ground representation */
	void gd_card(nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4);
	
	/*! \brief Standard radiation pattern parameters 
	
	\param calc_mode This integer selects the mode of calculation for the radiated field. Some values of (calc_mode) will affect the meaning of the remaining parameters on the card. Options available for calc_mode are:
		\arg \c O - normal mode. Space-wave fields are computed. An infinite ground plane is included if it has been specified previously on a GN card; otherwise, the antenna is in free space.
		\arg \c 1 - surface wave propagating along ground is added to the normal space wave. This option changes the meaning of some of the other parameters on the RP card as explained below, and the results appear in a special output format. Ground parameters must have been input on a GN card. The following options cause calculation of only the space wave but with special ground conditions. Ground conditions include a two-medium ground (cliff where the media join in a circle or a line), and a radial wire ground screen. Ground parameters and dimensions must be input on a GN or GD card before the RP card is read. The RP card only selects the option for inclusion in the field calculation. (Refer to the GN and GD cards for further explanation.)
		\arg \c  2 - linear cliff with antenna above upper level. Lower medium parameters are as specified for the second medium on the GN card or on the GD card.
		\arg \c 3 - circular cliff centered at origin of coordinate system: with antenna above upper level. Lower medium parameters are as specified for the second medium on the GN card or on the GD card.
		\arg \c 4 - radial wire ground screen centered at origin.
		\arg \c 5 - both radial wire ground screen and linear cliff.
		\arg \c 6 - both radial wire ground screen ant circular cliff.
	
	\param n_theta The number of theta angles.
	\param n_phi The number of phi angles. 
		 
	The field point is specified in spherical coordinates (R, sigma, theta), except when the surface wave is computed. For computing the surface wave field (calc_mode = l), cylindrical coordinates (phi, theta, z) are used to accurately define points near the ground plane at large radial distances.
		 
	The rp_card() function allows automatic stepping of the field point to compute the field over a region about the antenna at uniformly spaced points.
	The integers n_theta and n_phi and floating point numbers theta0, phi0, delta_theta, delta_phi, radial_distance, and gain_norm control the field-point stepping.
		 
	\param theta0 - Initial theta angle in degrees (initial z coordinate in meters if calc_mode = 1).
	
	\param phi0 - Initial phi angle in degrees.
	
	\param delta_theta - Increment for theta in degrees (increment for z in meters if calc_mode = 1).
	
	\param delta_phi - Increment for phi in degrees.
	
	\param radial_distance - Radial distance (R) of field point from the origin in meters. radial_distance is optional. If it is zero, the radiated electric field will have the factor exp(-jkR)/R omitted. If a value of R is specified, it should represent a point in the far-field region since near components of the field cannot be obtained with an RP card. (If calc_mode = 1, then radial_distance represents the cylindrical coordinate phi in meters and is not optional. It must be greater than about one wavelength.)
	
	\param gain_norm - Determines the gain normalization factor if normalization has been requested in the I4 field. If gain_norm is zero, the gain will be normalized to its maximum value. If gain_norm is not zero, the gain wi11 be normalized to the value of gain_norm.
	
	\remark
	\li The rp_card() function will call simulate(), causing the interaction matrix to be computed and factored and the structure currents to be computed if these operations have not already been performed. Hence, all required input parameters must be set before the fp_card() function is called. 
	\li At a single frequency, any number of rp_card() calls may occur in sequence so that different field-point spacings may be used over different regions of space. If automatic frequency stepping is being used (i.e., in_nfrq on the fr_card() function is greater than one), only one rp_card() function will act as data inside the loop. Subsequent calls to rp_card() will calculate patterns at the final frequency. 
	\li When both n_theta and n_phi are greater than one, the angle theta (or Z) will be stepped faster than phi. 
	\li When a ground plane has been specified, field points should not be requested below the ground (theta greater than 90 degrees or Z less than zero.)
	
	*/
	void rp_card(int calc_mode, int n_theta, int n_phi, int itmp4, nec_float theta0, nec_float phi0, nec_float delta_theta, nec_float delta_phi, nec_float radial_distance, nec_float gain_norm);
	
	 /*! "pt" card, print control for current */
	void pt_card(int itmp1, int itmp2, int itmp3, int itmp4);
	
	
	 /*! "pq" card, print control for charge */
	void pq_card(int itmp1, int itmp2, int itmp3, int itmp4);
	
	
	
	/*! "kh" card, matrix integration limit */
	void kh_card(nec_float tmp1);
	
	
	/*! Near field calculation parameters 
	
	\remark
	\li If the number of frequencies is not equal to one (as specified by the fr_card() function, then the ne_card() function will call simulate(), causing the interaction matrix to be computed and factored and the structure currents to be computed.
	*/
	void ne_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);
	
	/*! Near field calculation parameters 
	
	\remark
	\li If the number of frequencies is not equal to one (as specified by the fr_card() function, then the ne_card() function will call simulate(), causing the interaction matrix to be computed and factored and the structure currents to be computed.
	*/
	void nh_card(int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);
	
	/*! "ek" card,  extended thin wire kernel option */
	void ek_card(int itmp1);
	
	
	/*! "cp" card, maximum coupling between antennas */
	void cp_card(int itmp1, int itmp2, int itmp3, int itmp4);
	
	
	/*! "pl" card, plot flags 
		\exception int Throws int on error.
	*/
	void pl_card(const char* ploutput_filename, int itmp1, int itmp2, int itmp3, int itmp4);
	
	
	/*!****************************************************
	*** normal exit of nec2++ when all jobs complete ok ***
	******************************************************/
	void all_jobs_completed()
	{
		// put in here for the moment...
		m_results.write_to_file();
	}
	
	
	
	/*!	\brief Start a simulation
	
		This function will trigger a calculation. In the traditional NEC
		world, This signals the end of the main input section and the
		beginning of the frequency do loop.
		
		\param far_field_flag is true if last card was XQ or RP
		\warning far_field_flag is should never be specified as true
		because both the xq_card() and rp_card() functions will call
		this function automatically.
	*/
	void simulate(bool far_field_flag = false);

	//! an object to pipe output through...
	nec_output_file m_output;
	nec_ground ground;
	c_geometry geometry;
	c_plot_card plot_card;
	
	c_ggrid ggrid;
	c_ground_wave ground_wave;

//private:
	//! the maximum gain of the last radiation pattern
	nec_float m_maximum_gain_db;
	real_array m_radiation_pattern;
		
	//! some geometry temp flatgs
	int ifrtmw;
	int ifrtmp;

	
	//! pq card flags
	int iptflq;
	int iptaq, iptaqf, iptaqt;
	
	
	//! pt card flags...
	int iptflg;
	int iptag, iptagf, iptagt;
	
	
	int iflow;
	int ifrq, nfrq;
	nec_float delfrq;
	int_array ldtyp, ldtag, ldtagf, ldtagt;
	real_array zlr, zli, zlc, fnorm;
	
	int nthi, nphi;
	nec_float thetis, phiss;
	
	
	nec_results m_results;
	
	
	//! an object to pipe output through...
	nec_output_flags m_output_flags;
	
	
	nec_float wavelength;
	
	/* common  /cmb/ */
	complex_array cm;	// primary interaction matrix
	
	/* common  /matpar/ */
	int icase, npblk, nlast;
	int imat, nbbx, npbx, nlbx, nbbl, npbl, nlbl;
	
	/* common  /save/ */
	int_array ip;
	nec_float freq_mhz;
	
	/* common  /crnt/ */
	real_array air, aii; //! coefficients of the constant terms in the current interpolation functions for the current vector 
	real_array bir, bii; //! coefficients of the sine terms in the current interpolation functions
	real_array cir, cii; //! coefficients of the cosine terms in the current interpolation functions
	complex_array current_vector;	//! the current vector
	
	int ifar;
	nec_float t2, scrwl, scrwr;
	nec_complex zrati, t1, frati;
	
	/* common  /zload/ */
	int nload;
	complex_array zarray;
	
	/* common  /yparm/ */
	int ncoup, icoup;
	int_array nctag, ncseg;
	complex_array y11a, y12a;
	
	/* common  /vsorc/ */
	int_array ivqd, source_segment_array, iqds;
	int nvqd, voltage_source_count, nqds;
	complex_array vqd, vqds, source_voltage_array;
	
	/* common  /netcx/ */
	int masym, neq, npeq, neq2, network_count, ntsol, nprint;
	int_array iseg1, iseg2, ntyp;
	real_array x11r, x11i, x12r;
	real_array x12i, x22r, x22i;
	nec_float input_power, network_power_loss;
	nec_complex zped;
	
	/* common  /fpat/ */
	int near, nfeh, nrx, nry, nrz, nth, nph, ipd, iavp, inor, iax, excitation_type;
	nec_float thets, phis, dth, dph, rfld, gnor; 
	nec_float xpr6, structure_power_loss, xnr, ynr, znr, dxnr, dynr, dznr;
	
	
	/* common  /dataj/ */
	int iexk, ind1, indd1, ind2, indd2, ipgnd;
	nec_float s, b, xj, yj, zj, cabj, sabj, salpj;
	nec_float rkh; /* matrix integration limit */
	nec_float t1xj, t1yj, t1zj, t2xj, t2yj, t2zj;
	nec_complex  exk, eyk, ezk, exs, eys, ezs, exc, eyc, ezc;
	
	/* common  /smat/ */
	int nop; /* My addition */
	complex_array symmetry_array;
	
	/* common  /incom/ */
	int isnor;
	nec_float xo, yo, zo, sn, xsn, ysn;
	
	/* common  /tmi/ */
	int ija; /* changed to ija to avoid conflict */
	nec_float zpk, rkb2;
	
	/*common  /tmh/ */
	nec_float zpka, rhks;
	
private:

	/*! \brief A private convenience function called by ne_card() and nh_card()
	*/
	void ne_nh_card(int in_nfeh, int itmp1, int itmp2, int itmp3, int itmp4, nec_float tmp1, nec_float tmp2, nec_float tmp3, nec_float tmp4, nec_float tmp5, nec_float tmp6);
	
	
	void print_freq_int_krnl(
		nec_float f, 
		nec_float lambda, 
		nec_float int_dist, 
		bool using_extended_kernel);
		
	void 	antenna_env(void);
	void	print_structure_currents(char *pattype, int iptflg, int iptflq,
		real_array& fnorm, int iptag, int iptagf, int iptagt, int iptaq, 
		int iptaqf, int iptaqt);
	void	print_network_data(void);
	void 	print_norm_rx_pattern(int iptflg, int nthi, int nphi,
		real_array& fnorm, nec_float thetis, nec_float phiss);
	void	print_input_impedance(int iped, int ifrq, int nfrq, nec_float delfrq,
		real_array& fnorm);
	void	print_power_budget(void);
	void	structure_segment_loading(int_array& ldtyp, int_array& ldtag, int_array& ldtagf,
		int_array& ldtagt, real_array& zlr, real_array& zli, real_array& zlc);
		
	
	enum excitation_return
		excitation_loop(int in_freq_loop_state, int mhz, real_array& fnorm, 
			int iptflg, int iptflq, int iptag, int iptagf, int iptagt, 
			int iptaq, int iptaqf, int iptaqt, nec_float thetis, 
			int nfrq, int iflow, int nthi, int nphi, int iped, 
			int ib11, int ic11, int id11,
			int internal_inc);
			
	void	setup_excitation(int iptflg);
	
	
private:
	/* pointers to output files */
	FILE *m_output_fp;
	
	int inc, processing_state, isave;
	int nthic, nphic;
	int iped;
	
	nec_float impedance_norm_factor; // was zpnorm
	
	nec_float xpr1, xpr2, xpr3, xpr4, xpr5, xpr7;
	
	
	void load(int_array& ldtyp, int_array& ldtag,
			int_array& ldtagf, int_array& ldtagt,
			real_array& zlr, real_array& zli, real_array& zlc);

	void cmset(int nrow, complex_array& cm, nec_float rkhx, int iexkx);
	void compute_matrix_ss(int j1, int j2, int im1, int im2,
			complex_array& cm, int nrow, int itrp);
	void cmsw(int j1, int j2, int i1, int i2, complex_array& cm,
			complex_array& cw, int ncw, int nrow, int itrp);
	void cmws(int j, int i1, int i2, complex_array& cm, int nr,
			complex_array& cw, int nw, int itrp);
	void cmww(int j, int i1, int i2, complex_array& cm, int nr,
			complex_array& cw, int nw, int itrp);
	void couple(complex_array& cur, nec_float wlam);

	void efld(nec_float xi, nec_float yi, nec_float zi, nec_float ai, int ij);
	void eksc(nec_float s, nec_float z, nec_float rh, nec_float xk, int ij,
			nec_complex *ezs, nec_complex *ers, nec_complex *ezc,
			nec_complex *erc, nec_complex *ezk, nec_complex *erk);
	void ekscx(nec_float bx, nec_float s, nec_float z, nec_float rhx, nec_float xk,
			int ij, int inx1, int inx2, nec_complex *ezs,
			nec_complex *ers, nec_complex *ezc, nec_complex *erc,
			nec_complex *ezk, nec_complex *erk);
	void etmns(nec_float p1, nec_float p2, nec_float p3, nec_float p4, nec_float p5,
			nec_float p6, nec_float incident_amplitude, int excitation_type, complex_array& e);

	void fblock( int nrow, int ncol, int imax, int ipsym );

	void gf(nec_float zk, nec_float *co, nec_float *si);
	void gh(nec_float zk, nec_float *hr, nec_float *hi);
	void gx(nec_float zz, nec_float rh, nec_float xk,
			nec_complex *gz, nec_complex *gzp);
	void gxx(nec_float zz, nec_float rh, nec_float a, nec_float a2, nec_float xk,
			int ira, nec_complex *g1, nec_complex *g1p, nec_complex *g2,
			nec_complex *g2p, nec_complex *g3, nec_complex *gzp);
	void hfk(nec_float el1, nec_float el2, nec_float rhk,
			nec_float zpkx, nec_float *sgr, nec_float *sgi);
	void hintg(nec_float xi, nec_float yi, nec_float zi);
	void hsfld(nec_float xi, nec_float yi, nec_float zi, nec_float ai);
	void hsflx(nec_float s, nec_float rh, nec_float zpx, nec_complex *hpk,
			nec_complex *hps, nec_complex *hpc);

	void intx(nec_float el1, nec_float el2, nec_float b, int ij,
			nec_float *sgr, nec_float *sgi);

	void nefld(nec_float xob, nec_float yob, nec_float zob, nec_complex *ex,
			nec_complex *ey, nec_complex *ez);
	void netwk(complex_array& cm, nec_complex *cmb, nec_complex *cmc,
			nec_complex *cmd, int_array& ip, complex_array& einc);
	void nfpat(void);
	void nhfld(nec_float xob, nec_float yob, nec_float zob, nec_complex *hx,
			nec_complex *hy, nec_complex *hz);
	void pcint(nec_float xi, nec_float yi, nec_float zi, nec_float cabi,
			nec_float sabi, nec_float salpi, complex_array& e);
	void impedance_print(int in1, int in2, int in3, nec_float fl1, nec_float fl2,
			nec_float fl3, nec_float fl4, nec_float fl5, nec_float fl6, 
			char *ia, int ichar);
	void qdsrc(int is, nec_complex v, complex_array& e);
	void print_radiation_pattern(nec_float pinr, nec_float pnlr);

	
	void rom2(nec_float a, nec_float b, complex_array& sum, nec_float dmin);
	void sflds(nec_float t, complex_array& e);
	void solgf(nec_complex *a, nec_complex *b, nec_complex *c,
			nec_complex *d, nec_complex *xy, int *ip, int np, int n1,
			int n, int mp, int m1, int m, int n1c, int n2c, int n2cz);
	void unere(nec_float xob, nec_float yob, nec_float zob);
	nec_complex zint(nec_float sigl, nec_float rolam);

}; /* nec_context */
