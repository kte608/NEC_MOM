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
#include "c_geometry.h"

#include "nec_context.h"
#include "nec_exception.h"

c_geometry::c_geometry()
{
	n = 0;
	np = 0; 	// n is the number of segments
	
	m = 0;
	mp = 0;		// m is the number of patches?
	
	m_ipsym = 0;
	
	n_plus_m = 0;
	n_plus_2m = 0;
	n_plus_3m = 0;
	
	jsno = 0;
	nscon = 0;
	maxcon = 0;
	
	m_context = NULL;
	m_output = NULL;
	m_output_fp = NULL;
}

void c_geometry::set_context(nec_context* in_context)
{
	m_context = in_context;
	m_output = &m_context->m_output;
	m_output_fp = m_output->get_fp();
}

/*
	get_segment_number returns the segment number of the mth segment having the
	tag number in_tag.  if in_tag=0 segment number m is returned.
*/
int c_geometry::get_segment_number( int in_tag, int m)
{
	int segment_count = n;
	
	ASSERT(in_tag >= 0);
	ASSERT(m >= 0);
	ASSERT(segment_count >= 0);
	
	
	if (m <= 0)
	{
		throw new nec_exception("CHECK DATA, PARAMETER SPECIFYING SEGMENT POSITION IN A GROUP OF EQUAL TAGS MUST NOT BE ZERO" );
	}
	
	if ( 0 == in_tag)
	{
		return( m );
	}
	
	int tag_seg_count=0;
	for (int i = 0; i < segment_count; i++ )
	{
		if ( segment_tags[i] == in_tag )
		{	
			tag_seg_count++;
			if ( tag_seg_count == m)
			{
				return( i+1 );
			}
		}		
	}
	
	throw new nec_exception("NO SEGMENT HAS AN ITAG OF ", in_tag);
	
	return 0;
}




/*-----------------------------------------------------------------------*/
/*
	parse_geometry is the main routine for input of geometry data.

	ANTLR Grammar for geometry input
	
		geometry
			:	(geometry_element {transform}*)+
				geometry_end
			;
		
		geometry_element
			:	wire
			|	arc
			|	helix
			|	patch
			|	multi_patch
			;
		
		transform
			:	scale
			|	reflect
			|	rotate
			|	move
			;
		
		wire = "GW" int int real real real real real real {taper}
		taper = "GC" xxx
		arc = "GA"
		helix = "GH"
		
		patch = "SP"
		multi_patch = "SM" {patch_data}+
		patch_data = "SC"
		
		scale = "GS"
		reflect = "GX"
		rotate = "GR"
		move = "GM"
		
		geometry_end = "GE" .+
		
		int:
			[0-9]+;
			
		real:
			int:{.int};
		
*/

#include "c_plot_card.h"
	

void c_geometry::parse_geometry(nec_context* m_context, FILE* input_fp )
{
	char gm[3];
	char ifx[2] = {'*', 'X'}, ify[2]={'*','Y'}, ifz[2]={'*','Z'};
	char ipt[4] = { 'P', 'R', 'T', 'Q' };
	
	/* input card mnemonic list */
	/* "XT" stands for "exit", added for testing */
	#define GM_NUM  12
	char *atst[GM_NUM] =
	{
		"GW", "GX", "GR", "GS", "GE", "GM", \
		"SP", "SM", "GA", "SC", "GH", "GF"
	};
	
	bool print_structure_spec = true;
	
	int nwire, isct, i1, i2, iy, iz;
	int ix; 
	int card_int_1, card_int_2; /* The two integer parameters from the geometry card */
	nec_float rad, xs1, xs2, ys1, ys2, zs1, zs2, x4=0, y4=0, z4=0;
	nec_float x3=0, y3=0, z3=0, xw1, xw2, yw1, yw2, zw1, zw2;
	nec_float dummy;
	
	m_ipsym=0;
	nwire=0;
	n=0;
	np=0;
	m=0;
	mp=0;
	isct=0;
	
	/* read geometry data card and branch to */
	/* section for operation requested */
	do
	{
		read_geometry_card(input_fp, gm, &card_int_1, &card_int_2, &xw1, &yw1, &zw1, &xw2, &yw2, &zw2, &rad);
	
		/* identify card id mnemonic */
		int gm_num;
		for( gm_num = 0; gm_num < GM_NUM; gm_num++ )
			if ( strncmp( gm, atst[gm_num], 2) == 0 )
				break;

		if ( print_structure_spec )
		{
			m_output->end_section();
			m_output->set_indent(32);
			m_output->line("-------- STRUCTURE SPECIFICATION --------");
			m_output->line("COORDINATES MUST BE INPUT IN" );
			m_output->line("METERS OR BE SCALED TO METERS" );
			m_output->line("BEFORE STRUCTURE INPUT IS ENDED" );
			m_output->set_indent(0);
		
			m_output->line("  WIRE                                                                                 SEG FIRST  LAST  TAG");
			m_output->line("   No:        X1         Y1         Z1         X2         Y2         Z2       RADIUS   No:   SEG   SEG  No:" );
		
			print_structure_spec = false;
		}

		if ( gm_num != 10 )
			isct=0;

		switch( gm_num )
		{
		case 0:
		/* "gw" card, generate segment data for straight wire.
			GW		STRAIGHT WIRE, ENDS 1,2
				card_int_1- TAG NO.
				card_int_2- NO. SEGMENTS
				xw1- X1
				F2- Y1
				F3- Z1
				F4- X2
				F5- Y2
				F6- Z2
				F7- WIRE RAD., 0=USE GC FOR TAPERED WIRE
		*/
		{
			int wire_segment_count = card_int_2;
			int wire_tag = card_int_1;
			
			nwire++;
		
			// output some wire diagnostics.
			m_output->nec_printf( "\n"
				" %5d  %10.5f %10.5f %10.5f %10.5f"
				" %10.5f %10.5f %10.5f %5d %5d %5d %4d",
				nwire, xw1, yw1, zw1, xw2, yw2, zw2, rad, wire_segment_count, n+1, n + wire_segment_count, wire_tag );
		
			if ( rad != 0)	// rad == 0 implies a tapered wire
			{
				xs1 = 1.0;
				ys1 = 1.0;
			}
			else
			{
				read_geometry_card(input_fp, gm, &ix, &iy, &xs1, &ys1, &zs1,
					&dummy, &dummy, &dummy, &dummy);
			
				if ( strcmp(gm, "GC" ) != 0 )
				{
					throw new nec_exception("GEOMETRY DATA CARD ERROR" );
				}
			
				m_output->nec_printf(
					"\n  ABOVE WIRE IS TAPERED.  SEGMENT LENGTH RATIO: %9.5f\n"
					"                                 "
					"RADIUS FROM: %9.5f TO: %9.5f", xs1, ys1, zs1 );
			
				if ( (ys1 == 0) || (zs1 == 0) )
				{
					throw new nec_exception("GEOMETRY DATA CARD ERROR" );
				}
			
				rad= ys1;
				ys1= pow( (zs1/ys1), (1./(wire_segment_count-1.)) );
			}
			//printf("radius:%g\n",rad);
			wire(wire_tag, wire_segment_count, xw1, yw1, zw1, xw2, yw2, zw2, rad, xs1, ys1);
		}
		continue;

		/* reflect structure along x,y, or z */
		/* axes or rotate to form cylinder.  */
		case 1: /* "gx" card */
			iy= card_int_2/10;
			iz= card_int_2- iy*10;
			ix= iy/10;
			iy= iy- ix*10;
		
			if ( ix != 0)
			ix=1;
			if ( iy != 0)
			iy=1;
			if ( iz != 0)
			iz=1;
		
			m_output->nec_printf(
				"\n  STRUCTURE REFLECTED ALONG THE AXES %c %c %c"
				" - TAGS INCREMENTED BY %d\n",
				ifx[ix], ify[iy], ifz[iz], card_int_1 );
		
			reflect( ix, iy, iz, card_int_1, card_int_2);
		
			continue;
	
		case 2: /* "gr" card */
		
			m_output->nec_printf(
				"\n  STRUCTURE ROTATED ABOUT Z-AXIS %d TIMES"
				" - LABELS INCREMENTED BY %d\n", card_int_2, card_int_1 );
		
			ix=-1;
			iz = 0;
			reflect( ix, iy, iz, card_int_1, card_int_2);
		
			continue;
	
		case 3: /* "gs" card, scale structure dimensions by factor xw1. */
		
			if ( n > 0)
			{
				for(int i = 0; i < n; i++ )
				{
					x[i]= x[i]* xw1;
					y[i]= y[i]* xw1;
					z[i]= z[i]* xw1;
					x2[i]= x2[i]* xw1;
					y2[i]= y2[i]* xw1;
					z2[i]= z2[i]* xw1;
					bi[i]= bi[i]* xw1;
				}
			} /* if ( n >= n2) */
		
			if ( m > 0)
			{
				yw1= xw1* xw1;
				for(int i = 0; i < m; i++ )
				{
					px[i]= px[i]* xw1;
					py[i]= py[i]* xw1;
					pz[i]= pz[i]* xw1;
					pbi[i]= pbi[i]* yw1;
				}
			} /* if ( m >= m2) */
		
			m_output->nec_printf(
				"\n     STRUCTURE SCALED BY FACTOR: %10.5f", xw1 );
		
			continue;

		case 4: /* "ge" card, terminate structure geometry input. */

			geometry_complete(m_context, card_int_1, card_int_2);
			
			return;

		/* "gm" card, move structure or reproduce */
		/* original structure in new positions.   */
		case 5:
		
			m_output->nec_printf(
				"\n     THE STRUCTURE HAS BEEN MOVED, MOVE DATA CARD IS:\n"
				"   %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f",
				card_int_1, card_int_2, xw1, yw1, zw1, xw2, yw2, zw2, rad );
		
			xw1= xw1* TA;
			yw1= yw1* TA;
			zw1= zw1* TA;
		
			move( xw1, yw1, zw1, xw2, yw2, zw2, (int)( rad+.5), card_int_2, card_int_1);
			continue;
		
		case 6: /* "sp" card, generate single new patch */
		
			i1= m+1;
			card_int_2++;
		
			if ( card_int_1 != 0)
			{
				throw new nec_exception("PATCH DATA ERROR" );
			}
		
			m_output->nec_printf( "\n"
				" %5d%c %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f",
				i1, ipt[card_int_2-1], xw1, yw1, zw1, xw2, yw2, zw2 );
		
			if ( (card_int_2 == 2) || (card_int_2 == 4) )
			isct=1;
		
			if ( card_int_2 > 1)
			{
				read_geometry_card(input_fp, gm, &ix, &iy, &x3, &y3, &z3, &x4, &y4, &z4, &dummy);
			
				if ( (card_int_2 == 2) || (card_int_1 > 0) )
				{
					x4= xw1+ x3- xw2;
					y4= yw1+ y3- yw2;
					z4= zw1+ z3- zw2;
				}
			
				m_output->nec_printf( "\n"
					"      %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f",
					x3, y3, z3, x4, y4, z4 );
			
				if ( strcmp(gm, "SC") != 0 )
				{
					throw new nec_exception("PATCH DATA ERROR" );
				}
			} /* if ( card_int_2 > 1) */
			else
			{
				xw2= xw2* TA;
				yw2= yw2* TA;
			}
		
			patch( card_int_1, card_int_2, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);
			continue;
	
		case 7: /* "sm" card, generate multiple-patch surface */
		
			i1= m+1;
			m_output->nec_printf( "\n"
				" %5d%c %10.5f %11.5f %11.5f %11.5f %11.5f %11.5f"
				"     SURFACE - %d BY %d PATCHES",
				i1, ipt[1], xw1, yw1, zw1, xw2, yw2, zw2, card_int_1, card_int_2 );
		
			if ( (card_int_1 < 1) || (card_int_2 < 1) )
			{
				throw new nec_exception("PATCH DATA ERROR" );
			}
		
			read_geometry_card(input_fp, gm, &ix, &iy, &x3, &y3, &z3, &x4, &y4, &z4, &dummy);
		
			if ( (card_int_2 == 2) || (card_int_1 > 0) )
			{
				x4= xw1+ x3- xw2;
				y4= yw1+ y3- yw2;
				z4= zw1+ z3- zw2;
			}
		
			m_output->nec_printf( "\n"
				"      %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f",
				x3, y3, z3, x4, y4, z4 );
		
			if ( strcmp(gm, "SC" ) != 0 )
			{
				throw new nec_exception("PATCH DATA ERROR" );
			}
		
			patch( card_int_1, card_int_2, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);
			continue;
		
		case 8: /* "ga" card, generate segment data for wire arc */
		
			nwire++;
			i1= n+1;
			i2= n+ card_int_2;
		
			m_output->nec_printf( "\n"
				" %5d  ARC RADIUS: %9.5f  FROM: %8.3f TO: %8.3f DEGREES"
				"       %11.5f %5d %5d %5d %4d",
				nwire, xw1, yw1, zw1, xw2, card_int_2, i1, i2, card_int_1 );
		
			arc( card_int_1, card_int_2, xw1, yw1, zw1, xw2);
			continue;
	
		case 9: /* "sc" card */
		
			if ( isct == 0)
			{
				throw new nec_exception("PATCH DATA ERROR" );
			}
		
			i1= m+1;
			card_int_2++;
		
			if ( (card_int_1 != 0) || ((card_int_2 != 2) && (card_int_2 != 4)) )
			{
				throw new nec_exception("PATCH DATA ERROR" );
			}
		
			xs1= x4;
			ys1= y4;
			zs1= z4;
			xs2= x3;
			ys2= y3;
			zs2= z3;
			x3= xw1;
			y3= yw1;
			z3= zw1;
		
			if ( card_int_2 == 4)
			{
				x4= xw2;
				y4= yw2;
				z4= zw2;
			}
		
			xw1= xs1;
			yw1= ys1;
			zw1= zs1;
			xw2= xs2;
			yw2= ys2;
			zw2= zs2;
		
			if ( card_int_2 != 4)
			{
				x4= xw1+ x3- xw2;
				y4= yw1+ y3- yw2;
				z4= zw1+ z3- zw2;
			}
		
			m_output->nec_printf( "\n"
				" %5d%c %10.5f %11.5f %11.5f %11.5f %11.5f %11.5f",
				i1, ipt[card_int_2-1], xw1, yw1, zw1, xw2, yw2, zw2 );
		
			m_output->nec_printf( "\n"
				"      %11.5f %11.5f %11.5f  %11.5f %11.5f %11.5f",
				x3, y3, z3, x4, y4, z4 );
		
			patch( card_int_1, card_int_2, xw1, yw1, zw1, xw2, yw2, zw2, x3, y3, z3, x4, y4, z4);
		
			continue;
	
		case 10: /* "gh" card, generate helix */
		
			nwire++;
			i1= n+1;
			i2= n+ card_int_2;
		
			m_output->nec_printf( "\n"
				" %5d HELIX STRUCTURE - SPACING OF TURNS: %8.3f AXIAL"
				" LENGTH: %8.3f  %8.3f %5d %5d %5d %4d\n      "
				" RADIUS X1:%8.3f Y1:%8.3f X2:%8.3f Y2:%8.3f ",
				nwire, xw1, yw1, rad, card_int_2, i1, i2, card_int_1, zw1, xw2, yw2, zw2 );
		
			helix( xw1, yw1, zw1, xw2, yw2, zw2, rad, card_int_2, card_int_1);
		
			continue;

		case 11: /* "gf" card, not supported */
			abort_on_error(-5);
	
		default: /* error message */
		
			m_output->nec_printf( "\n  GEOMETRY DATA CARD ERROR" );
			m_output->nec_printf( "\n"
				" %2s %3d %5d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f",
				gm, card_int_1, card_int_2, xw1, yw1, zw1, xw2, yw2, zw2, rad );
		
			stop(-1);
    } /* switch( gm_num ) */

	} /* do */
	while( true );
}

/**
	We have finished with the geometry description, now connect 
	things up.
*/
void c_geometry::geometry_complete(nec_context* m_context, int card_int_1, int card_int_2)
{
	if (0 == np + mp)
		throw new nec_exception("Geometry has no wires or patches.");
		
	/* TCAM: The following does not make sense for the semantics of the plot
	card. I have left this in, I hope someone will tell me why it is required.
	*/
	if (card_int_2 != 0)
		m_context->plot_card.set_plot_real_imag_currents();
		
	// now proceed and complete the geometry setup...
	
	connect_segments( card_int_1);

	if ( n != 0)
	{
		/* Allocate wire buffers */
		segment_length.resize(n);
		sab.resize(n);
		cab.resize(n);
		salp.resize(n);
	
		m_output->nec_printf( "\n\n\n"
			"                              "
			" ---------- SEGMENTATION DATA ----------\n"
			"                                       "
			" COORDINATES IN METERS\n"
			"                           "
			" I+ AND I- INDICATE THE SEGMENTS BEFORE AND AFTER I\n" );
	
		m_output->nec_printf( "\n"
			"   SEG    COORDINATES OF SEGM CENTER     SEGM    ORIENTATION"
			" ANGLES    WIRE    CONNECTION DATA   TAG\n"
			"   No:       X         Y         Z      LENGTH     ALPHA     "
			" BETA    RADIUS    I-     I    I+   NO:" );
	
		for(int i = 0; i < n; i++ )
		{
			nec_float xw1= x2[i]- x[i];
			nec_float yw1= y2[i]- y[i];
			nec_float zw1= z2[i]- z[i];
			x[i]=( x[i]+ x2[i])*.5;
			y[i]=( y[i]+ y2[i])*.5;
			z[i]=( z[i]+ z2[i])*.5;
			
			nec_float xw2= xw1* xw1+ yw1* yw1+ zw1* zw1;
			nec_float yw2= sqrt( xw2);
			yw2=( xw2/ yw2+ yw2)*.5;
			segment_length[i]= yw2;
			cab[i]= xw1/ yw2;
			sab[i]= yw1/ yw2;
			xw2= zw1/ yw2;
	
			if ( xw2 > 1.)
				xw2=1.;
			if ( xw2 < -1.)
				xw2=-1.;
	
			salp[i]= xw2;
			xw2= rad_to_degrees(asin( xw2));
			yw2= rad_to_degrees(atan2( yw1, xw1));
	
			m_output->nec_printf( "\n"
				" %5d %9.5f %9.5f %9.5f %9.5f"
				" %9.5f %9.5f %9.5f %5d %5d %5d %5d",
				i+1, x[i], y[i], z[i], segment_length[i], xw2, yw2,
				bi[i], icon1[i], i+1, icon2[i], segment_tags[i] );
			//printf("I am here! %g\n",bi[i]);
			m_context->plot_card.plot_segments(i,x,y,z,segment_length,xw2,yw2,bi,icon1,icon2);
	
			if ( (segment_length[i] <= 1.e-20) || (bi[i] <= 0.) )
			{
				throw new nec_exception("SEGMENT DATA ERROR" );
			}
		} /* for( i = 0; i < n; i++ ) */
	} /* if ( n != 0) */

	if ( m != 0)
	{
		m_output->nec_printf( "\n\n\n"
			"                                   "
			" --------- SURFACE PATCH DATA ---------\n"
			"                                            "
			" COORDINATES IN METERS\n\n"
			" PATCH      COORD. OF PATCH CENTER           UNIT NORMAL VECTOR      "
			" PATCH           COMPONENTS OF UNIT TANGENT VECTORS\n"
			"  NO:       X          Y          Z          X        Y        Z      "
			" AREA         X1       Y1       Z1        X2       Y2      Z2" );
	
		for(int i = 0; i < m; i++ )
		{
			nec_float xw1=( t1y[i]* t2z[i]- t1z[i]* t2y[i])* psalp[i];
			nec_float yw1=( t1z[i]* t2x[i]- t1x[i]* t2z[i])* psalp[i];
			nec_float zw1=( t1x[i]* t2y[i]- t1y[i]* t2x[i])* psalp[i];
	
			m_output->nec_printf( "\n"
				" %4d %10.5f %10.5f %10.5f  %8.5f %8.5f %8.5f"
				" %10.5f  %8.5f %8.5f %8.5f  %8.5f %8.5f %8.5f",
				i+1, px[i], py[i], pz[i], xw1, yw1, zw1, pbi[i],
				t1x[i], t1y[i], t1z[i], t2x[i], t2y[i], t2z[i] );
		} /* for( i = 0; i < m; i++ ) */
	} /* if ( m == 0) */

	n_plus_m  = n+m;
	n_plus_2m = n+2*m;
	n_plus_3m = n+3*m;
}


/* subroutine wire generates segment geometry */
/* data for a straight wire of segment_count segments. */
void c_geometry::wire( int tag_id, int segment_count, nec_float xw1, nec_float yw1, nec_float zw1,
    nec_float xw2, nec_float yw2, nec_float zw2, nec_float rad,
    nec_float rdel, nec_float rrad )
{
	nec_float xd, yd, zd, delz, rd, fns, radz;
	nec_float xs1, ys1, zs1, xs2, ys2, zs2;
	
	int istart = n;
	n = n + segment_count;
	np= n;
	mp= m;
	m_ipsym=0;
	
	if ( segment_count < 1)
		return;
	
	/* Reallocate tags buffer */
	segment_tags.resize(n + m);
	// segment_tags.resize((n+m) );/*????*/
	
	/* Reallocate wire buffers */
	x.resize(n);
	y.resize(n);
	z.resize(n);
	x2.resize(n);
	y2.resize(n);
	z2.resize(n);
	bi.resize(n);
	
	xd= xw2- xw1;
	yd= yw2- yw1;
	zd= zw2- zw1;
	
	if ( fabs( rdel-1.) >= 1.0e-6)
	{
		delz= sqrt( xd* xd+ yd* yd+ zd* zd);
		xd= xd/ delz;
		yd= yd/ delz;
		zd= zd/ delz;
		delz= delz*(1.- rdel)/(1.- pow(rdel, segment_count) );
		rd= rdel;
	}
	else
	{
		fns= segment_count;
		xd= xd/ fns;
		yd= yd/ fns;
		zd= zd/ fns;
		delz=1.0;
		rd=1.0;
	}
	
	radz= rad;
	xs1= xw1;
	ys1= yw1;
	zs1= zw1;
	
	for(int i = istart; i < n; i++ )
	{
		segment_tags[i]= tag_id;
		xs2= xs1+ xd* delz;
		ys2= ys1+ yd* delz;
		zs2= zs1+ zd* delz;
		x[i]= xs1;
		y[i]= ys1;
		z[i]= zs1;
		x2[i]= xs2;
		y2[i]= ys2;
		z2[i]= zs2;
		bi[i]= radz;
		delz= delz* rd;
		radz= radz* rrad;
		xs1= xs2;
		ys1= ys2;
		zs1= zs2;
	}
	
	x2[n-1]= xw2;
	y2[n-1]= yw2;
	z2[n-1]= zw2;
}

/*-----------------------------------------------------------------------*/

/* subroutine helix generates segment geometry */
/* data for a helix of segment_count segments */
void c_geometry::helix( nec_float s, nec_float hl, nec_float a1, nec_float b1,
    nec_float a2, nec_float b2, nec_float rad, int segment_count, int tag_id )
{
	int ist;
	nec_float turns, zinc, copy, sangle, hdia, turn, pitch, hmaj, hmin;
	
	ist= n;
	n += segment_count;
	np= n;
	mp= m;
	m_ipsym=0;
	
	if ( segment_count < 1)
		return;
	
	turns= fabs( hl/ s);
	zinc= fabs( hl/ segment_count);
	
	segment_tags.resize(n+m); /*????*/
	
	/* Reallocate wire buffers */
	x.resize(n);
	y.resize(n);
	z.resize(n);
	x2.resize(n);
	y2.resize(n);
	z2.resize(n);
	bi.resize(n);
	
	z[ist]=0.;
	for(int i = ist; i < n; i++ )
	{
		bi[i]= rad;
		segment_tags[i]= tag_id;
	
		if ( i != ist )
		z[i]= z[i-1]+ zinc;
	
		z2[i]= z[i]+ zinc;
	
		if ( a2 == a1)
		{
			if ( b1 == 0.)
				b1= a1;
		
			x[i]= a1* cos(2.* pi()* z[i]/ s);
			y[i]= b1* sin(2.* pi()* z[i]/ s);
			x2[i]= a1* cos(2.* pi()* z2[i]/ s);
			y2[i]= b1* sin(2.* pi()* z2[i]/ s);
		}
		else
		{
			if ( b2 == 0.)
				b2= a2;
		
			x[i]=( a1+( a2- a1)* z[i]/ fabs( hl))* cos(2.* pi()* z[i]/ s);
			y[i]=( b1+( b2- b1)* z[i]/ fabs( hl))* sin(2.* pi()* z[i]/ s);
			x2[i]=( a1+( a2- a1)* z2[i]/ fabs( hl))* cos(2.* pi()* z2[i]/ s);
			y2[i]=( b1+( b2- b1)* z2[i]/ fabs( hl))* sin(2.* pi()* z2[i]/ s);	
		} /* if ( a2 == a1) */
	
		if ( hl > 0.)
			continue;
	
		copy= x[i];
		x[i]= y[i];
		y[i]= copy;
		copy= x2[i];
		x2[i]= y2[i];
		y2[i]= copy;	
	} /* for( i = ist; i < n; i++ ) */
	
	if ( a2 != a1)
	{
		sangle= atan( a2/( fabs( hl)+( fabs( hl)* a1)/( a2- a1)));
		m_output->nec_printf(
			"\n       THE CONE ANGLE OF THE SPIRAL IS %10.5f", sangle );
		return;
	}
	
	if ( a1 == b1)
	{
		hdia=2.* a1;
		turn= hdia* pi();
		pitch= atan( s/( pi()* hdia));
		turn= turn/ cos( pitch);
		pitch=180.* pitch/ pi();
	}
	else
	{
		if ( a1 >= b1)
		{
			hmaj=2.* a1;
			hmin=2.* b1;
		}
		else
		{
			hmaj=2.* b1;
			hmin=2.* a1;
		}
	
		hdia= sqrt(( hmaj*hmaj+ hmin*hmin)/2* hmaj);
		turn=2.* pi()* hdia;
		pitch=(180./ pi())* atan( s/( pi()* hdia));
	
	} /* if ( a1 == b1) */
	
	m_output->nec_printf( "\n"
		"       THE PITCH ANGLE IS: %.5f    THE LENGTH OF WIRE/TURN IS: %.5f",
		pitch, turn );
}

/*-----------------------------------------------------------------------*/


/* subroutine move moves the structure with respect to its */
/* coordinate system or reproduces structure in new positions. */
/* structure is rotated about x,y,z axes by rox,roy,roz */
/* respectively, then shifted by xs,ys,zs */
void c_geometry::move( nec_float rox, nec_float roy, nec_float roz, nec_float xs,
    nec_float ys, nec_float zs, int its, int nrpt, int itgi )
{
  int nrp, ix, i1, k, ir, i, ii;
  nec_float sps, cps, sth, cth, sph, cph, xx, xy;
  nec_float xz, yx, yy, yz, zx, zy, zz, xi, yi, zi;

  if ( fabs( rox)+ fabs( roy) > 1.0e-10)
    m_ipsym= m_ipsym*3;

  sps= sin( rox);
  cps= cos( rox);
  sth= sin( roy);
  cth= cos( roy);
  sph= sin( roz);
  cph= cos( roz);
  xx= cph* cth;
  xy= cph* sth* sps- sph* cps;
  xz= cph* sth* cps+ sph* sps;
  yx= sph* cth;
  yy= sph* sth* sps+ cph* cps;
  yz= sph* sth* cps- cph* sps;
  zx=- sth;
  zy= cth* sps;
  zz= cth* cps;

  if ( nrpt == 0)
    nrp=1;
  else
    nrp= nrpt;

  ix=1;
  if ( n > 0)
  {
    i1= get_segment_number( its, 1);
    if ( i1 < 1)
      i1= 1;

    ix= i1;
    if ( nrpt == 0)
      k= i1-1;
    else
    {
      k= n;
      /* Reallocate tags buffer */
      segment_tags.resize(n+m + (n+1-i1)*nrpt);
      // mreq = n+m + (n+1-i1)*nrpt;
      // segment_tags.resize(mreq);

      /* Reallocate wire buffers */
      int new_size = (n+(n+1-i1)*nrpt);
      x.resize(new_size);
      y.resize(new_size);
      z.resize(new_size);
      x2.resize(new_size);
      y2.resize(new_size);
      z2.resize(new_size);
      bi.resize(new_size);
    }

    for( ir = 0; ir < nrp; ir++ )
    {
      for( i = i1-1; i < n; i++ )
      {
	xi= x[i];
	yi= y[i];
	zi= z[i];
	x[k]= xi* xx+ yi* xy+ zi* xz+ xs;
	y[k]= xi* yx+ yi* yy+ zi* yz+ ys;
	z[k]= xi* zx+ yi* zy+ zi* zz+ zs;
	xi= x2[i];
	yi= y2[i];
	zi= z2[i];
	x2[k]= xi* xx+ yi* xy+ zi* xz+ xs;
	y2[k]= xi* yx+ yi* yy+ zi* yz+ ys;
	z2[k]= xi* zx+ yi* zy+ zi* zz+ zs;
	bi[k]= bi[i];
	segment_tags[k]= segment_tags[i];
	if ( segment_tags[i] != 0)
	  segment_tags[k]= segment_tags[i]+ itgi;

	k++;

      } /* for( i = i1; i < n; i++ ) */

      i1= n+1;
      n= k;

    } /* for( ir = 0; ir < nrp; ir++ ) */

  } /* if ( n >= n2) */

  if ( m > 0)
  {
    i1 = 0;
    if ( nrpt == 0)
      k= 0;
    else
      k = m;

    /* Reallocate patch buffers */
    int new_size = m * (1+nrpt);
    px.resize(new_size);
    py.resize(new_size);
    pz.resize(new_size);
    t1x.resize(new_size);
    t1y.resize(new_size);
    t1z.resize(new_size);
    t2x.resize(new_size);
    t2y.resize(new_size);
    t2z.resize(new_size);
    pbi.resize(new_size);
    psalp.resize(new_size);

    for( ii = 0; ii < nrp; ii++ )
    {
      for( i = i1; i < m; i++ )
      {
	xi= px[i];
	yi= py[i];
	zi= pz[i];
	px[k]= xi* xx+ yi* xy+ zi* xz+ xs;
	py[k]= xi* yx+ yi* yy+ zi* yz+ ys;
	pz[k]= xi* zx+ yi* zy+ zi* zz+ zs;
	xi= t1x[i];
	yi= t1y[i];
	zi= t1z[i];
	t1x[k]= xi* xx+ yi* xy+ zi* xz;
	t1y[k]= xi* yx+ yi* yy+ zi* yz;
	t1z[k]= xi* zx+ yi* zy+ zi* zz;
	xi= t2x[i];
	yi= t2y[i];
	zi= t2z[i];
	t2x[k]= xi* xx+ yi* xy+ zi* xz;
	t2y[k]= xi* yx+ yi* yy+ zi* yz;
	t2z[k]= xi* zx+ yi* zy+ zi* zz;
	psalp[k]= psalp[i];
	pbi[k]= pbi[i];
	k++;

      } /* for( i = i1; i < m; i++ ) */

      i1= m;
      m = k;

    } /* for( ii = 0; ii < nrp; ii++ ) */

  } /* if ( m >= m2) */

  if ( (nrpt == 0) && (ix == 1) )
    return;

  np= n;
  mp= m;
  m_ipsym=0;

  return;
}

/*-----------------------------------------------------------------------*/

/* reflects partial structure along x,y, or z axes or rotates */
/* structure to complete a symmetric structure. */
void c_geometry::reflect( int ix, int iy, int iz, int itx, int nop )
{
	int iti, i, nx, itagi, k;
	nec_float e1, e2, fnop, sam, cs, ss, xk, yk;
	
	np= n;
	mp= m;
	m_ipsym=0;
	iti= itx;
	
	if ( ix >= 0)
	{
		if ( nop == 0)
		return;
	
		m_ipsym=1;
	
		/* reflect along z axis */
		if ( iz != 0)
		{
		m_ipsym=2;
	
		if ( n > 0 )
		{
		/* Reallocate tags buffer */
		segment_tags.resize(2*n + m);
		// segment_tags.resize((2*n+m));
	
		/* Reallocate wire buffers */
		int new_size = 2*n;
		x.resize(new_size);
		y.resize(new_size);
		z.resize(new_size);
		x2.resize(new_size);
		y2.resize(new_size);
		z2.resize(new_size);
		bi.resize(new_size);
	
		for( i = 0; i < n; i++ )
		{
		nx= i+ n;
		e1= z[i];
		e2= z2[i];
	
		if ( (fabs(e1)+fabs(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
		{
			nec_exception* nex = new nec_exception("GEOMETRY DATA ERROR--SEGMENT ");
			nex->append(i+1);
			nex->append("LIES IN PLANE OF SYMMETRY");
			throw nex;
/*			m_output->nec_printf(
				"\n  GEOMETRY DATA ERROR--SEGMENT %d"
				" LIES IN PLANE OF SYMMETRY", i+1 );
				stop(-1); */
		}
	
		x[nx]= x[i];
		y[nx]= y[i];
		z[nx]=- e1;
		x2[nx]= x2[i];
		y2[nx]= y2[i];
		z2[nx]=- e2;
		itagi= segment_tags[i];
	
		if ( itagi == 0)
			segment_tags[nx]=0;
		if ( itagi != 0)
			segment_tags[nx]= itagi+ iti;
	
		bi[nx]= bi[i];
	
		} /* for( i = 0; i < n; i++ ) */
	
		n= n*2;
		iti= iti*2;
	
		} /* if ( n > 0) */
	
		if ( m > 0 )
		{
		/* Reallocate patch buffers */
		int new_size = 2*m;
		px.resize(new_size);
		py.resize(new_size);
		pz.resize(new_size);
		t1x.resize(new_size);
		t1y.resize(new_size);
		t1z.resize(new_size);
		t2x.resize(new_size);
		t2y.resize(new_size);
		t2z.resize(new_size);
		pbi.resize(new_size);
		psalp.resize(new_size);
	
		for( i = 0; i < m; i++ )
		{
		nx = i+m;
		if ( fabs(pz[i]) <= 1.0e-10)
		{
			m_output->nec_printf(
			"\n  GEOMETRY DATA ERROR--PATCH %d"
			" LIES IN PLANE OF SYMMETRY", i+1 );
			stop(-1);
		}
	
		px[nx]= px[i];
		py[nx]= py[i];
		pz[nx]=- pz[i];
		t1x[nx]= t1x[i];
		t1y[nx]= t1y[i];
		t1z[nx]=- t1z[i];
		t2x[nx]= t2x[i];
		t2y[nx]= t2y[i];
		t2z[nx]=- t2z[i];
		psalp[nx]=- psalp[i];
		pbi[nx]= pbi[i];
		}
	
		m= m*2;
	
		} /* if ( m >= m2) */
	
		} /* if ( iz != 0) */
	
		/* reflect along y axis */
		if ( iy != 0)
		{
		if ( n > 0)
		{
		/* Reallocate tags buffer */
		segment_tags.resize(2*n + m);
		// segment_tags.resize((2*n+m));/*????*/
	
		/* Reallocate wire buffers */
		int new_size = 2*n;
		x.resize(new_size);
		y.resize(new_size);
		z.resize(new_size);
		x2.resize(new_size);
		y2.resize(new_size);
		z2.resize(new_size);
		bi.resize(new_size);
	
		for( i = 0; i < n; i++ )
		{
			nx= i+ n;
			e1= y[i];
			e2= y2[i];
		
			if ( (fabs(e1)+fabs(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
			{
				m_output->nec_printf(
				"\n  GEOMETRY DATA ERROR--SEGMENT %d"
				" LIES IN PLANE OF SYMMETRY", i+1 );
				stop(-1);
			}
		
			x[nx]= x[i];
			y[nx]=- e1;
			z[nx]= z[i];
			x2[nx]= x2[i];
			y2[nx]=- e2;
			z2[nx]= z2[i];
			itagi= segment_tags[i];
		
			if ( itagi == 0)
				segment_tags[nx]=0;
			if ( itagi != 0)
				segment_tags[nx]= itagi+ iti;
		
			bi[nx]= bi[i];
		
		} /* for( i = n2-1; i < n; i++ ) */
	
		n= n*2;
		iti= iti*2;
	
		} /* if ( n >= n2) */
	
		if ( m > 0 )
		{
			/* Reallocate patch buffers */
			int new_size = 2*m;
			px.resize(new_size);
			py.resize(new_size);
			pz.resize(new_size);
			t1x.resize(new_size);
			t1y.resize(new_size);
			t1z.resize(new_size);
			t2x.resize(new_size);
			t2y.resize(new_size);
			t2z.resize(new_size);
			pbi.resize(new_size);
			psalp.resize(new_size);
		
			for( i = 0; i < m; i++ )
			{
				nx= i+m;
				if ( fabs( py[i]) <= 1.0e-10)
				{
					m_output->nec_printf(
					"\n  GEOMETRY DATA ERROR--PATCH %d"
					" LIES IN PLANE OF SYMMETRY", i+1 );
					stop(-1);
				}
			
				px[nx]= px[i];
				py[nx]=- py[i];
				pz[nx]= pz[i];
				t1x[nx]= t1x[i];
				t1y[nx]=- t1y[i];
				t1z[nx]= t1z[i];
				t2x[nx]= t2x[i];
				t2y[nx]=- t2y[i];
				t2z[nx]= t2z[i];
				psalp[nx]=- psalp[i];
				pbi[nx]= pbi[i];
			
			} /* for( i = m2; i <= m; i++ ) */
		
			m= m*2;
		} /* if ( m >= m2) */
	
		} /* if ( iy != 0) */
	
		/* reflect along x axis */
		if ( ix == 0 )
			return;
	
		if ( n > 0 )
		{
		/* Reallocate tags buffer */
		segment_tags.resize(2*n + m);
		// segment_tags.resize((2*n+m));/*????*/
	
		/* Reallocate wire buffers */
		int new_size = 2*n;
		x.resize(new_size);
		y.resize(new_size);
		z.resize(new_size);
		x2.resize(new_size);
		y2.resize(new_size);
		z2.resize(new_size);
		bi.resize(new_size);
	
		for( i = 0; i < n; i++ )
		{
			nx= i+ n;
			e1= x[i];
			e2= x2[i];
		
			if ( (fabs(e1)+fabs(e2) <= 1.0e-5) || (e1*e2 < -1.0e-6) )
			{
				m_output->nec_printf(
					"\n  GEOMETRY DATA ERROR--SEGMENT %d"
					" LIES IN PLANE OF SYMMETRY", i+1 );
				stop(-1);
			}
		
			x[nx]=- e1;
			y[nx]= y[i];
			z[nx]= z[i];
			x2[nx]=- e2;
			y2[nx]= y2[i];
			z2[nx]= z2[i];
			itagi= segment_tags[i];
		
			if ( itagi == 0)
			segment_tags[nx]=0;
			if ( itagi != 0)
			segment_tags[nx]= itagi+ iti;
		
			bi[nx]= bi[i];
		}
	
		n= n*2;
	
		} /* if ( n > 0) */
	
		if ( m == 0 )
			return;
	
		/* Reallocate patch buffers */
		int new_size = 2*m;
		px.resize(new_size);
		py.resize(new_size);
		pz.resize(new_size);
		t1x.resize(new_size);
		t1y.resize(new_size);
		t1z.resize(new_size);
		t2x.resize(new_size);
		t2y.resize(new_size);
		t2z.resize(new_size);
		pbi.resize(new_size);
		psalp.resize(new_size);
	
		for( i = 0; i < m; i++ )
		{
			nx= i+m;
			if ( fabs( px[i]) <= 1.0e-10)
			{
				m_output->nec_printf(
					"\n  GEOMETRY DATA ERROR--PATCH %d"
					" LIES IN PLANE OF SYMMETRY", i+1 );
				stop(-1);
			}
		
			px[nx]=- px[i];
			py[nx]= py[i];
			pz[nx]= pz[i];
			t1x[nx]=- t1x[i];
			t1y[nx]= t1y[i];
			t1z[nx]= t1z[i];
			t2x[nx]=- t2x[i];
			t2y[nx]= t2y[i];
			t2z[nx]= t2z[i];
			psalp[nx]=- psalp[i];
			pbi[nx]= pbi[i];
		}
	
		m= m*2;
		return;
	
	} /* if ( ix >= 0) */
	
	/* reproduce structure with rotation to form cylindrical structure */
	fnop= (nec_float)nop;
	m_ipsym=-1;
	sam=two_pi() / fnop;
	cs= cos( sam);
	ss= sin( sam);
	
	if ( n > 0)
	{
		n *= nop;
		nx= np;
	
		/* Reallocate tags buffer */
		segment_tags.resize(n + m);
		//segment_tags.resize((n+m));/*????*/
	
		/* Reallocate wire buffers */
		x.resize(n);
		y.resize(n);
		z.resize(n);
		x2.resize(n);
		y2.resize(n);
		z2.resize(n);
		bi.resize(n);
	
		for( i = nx; i < n; i++ )
		{
			k= i- np;
			xk= x[k];
			yk= y[k];
			x[i]= xk* cs- yk* ss;
			y[i]= xk* ss+ yk* cs;
			z[i]= z[k];
			xk= x2[k];
			yk= y2[k];
			x2[i]= xk* cs- yk* ss;
			y2[i]= xk* ss+ yk* cs;
			z2[i]= z2[k];
			bi[i]= bi[k];
			itagi= segment_tags[k];
		
			if ( itagi == 0)
				segment_tags[i]=0;
			if ( itagi != 0)
				segment_tags[i]= itagi+ iti;
		}
	
	} /* if ( n >= n2) */
	
	if ( m == 0 )
		return;
	
	m *= nop;
	nx= mp;
	
	/* Reallocate patch buffers */
	px.resize(m);
	py.resize(m);
	pz.resize(m);
	t1x.resize(m);
	t1y.resize(m);
	t1z.resize(m);
	t2x.resize(m);
	t2y.resize(m);
	t2z.resize(m);
	pbi.resize(m);
	psalp.resize(m);
	
	for( i = nx; i < m; i++ )
	{
		k = i-mp;
		xk= px[k];
		yk= py[k];
		px[i]= xk* cs- yk* ss;
		py[i]= xk* ss+ yk* cs;
		pz[i]= pz[k];
		xk= t1x[k];
		yk= t1y[k];
		t1x[i]= xk* cs- yk* ss;
		t1y[i]= xk* ss+ yk* cs;
		t1z[i]= t1z[k];
		xk= t2x[k];
		yk= t2y[k];
		t2x[i]= xk* cs- yk* ss;
		t2y[i]= xk* ss+ yk* cs;
		t2z[i]= t2z[k];
		psalp[i]= psalp[k];
		pbi[i]= pbi[k];
	
	} /* for( i = nx; i < m; i++ ) */
}
	
/*-----------------------------------------------------------------------*/

/* connect sets up segment connection data in arrays icon1 and */
/* icon2 by searching for segment ends that are in contact. */
void c_geometry::connect_segments( int ignd )
{
	bool jump, ipf;
	int i, iz, ic, j, jx, ix, ixx, iseg, iend, jend, nsflg;
	nec_float sep=0., xi1, yi1, zi1, xi2, yi2, zi2;
	nec_float slen, xa, ya, za, xs, ys, zs;
	
	nscon= -1;
	maxcon = 1;
	
	if ( ignd != 0)
	{
		m_output->nec_printf( "\n\n     GROUND PLANE SPECIFIED." );
	
		if ( ignd > 0)
			m_output->nec_printf(
				"\n     WHERE WIRE ENDS TOUCH GROUND, CURRENT WILL"
				" BE INTERPOLATED TO IMAGE IN GROUND PLANE.\n" );
	
		if ( m_ipsym == 2)
		{
			np=2* np;
			mp=2* mp;
		}
	
		if ( abs( m_ipsym) > 2 )
		{
			np= n;
			mp= m;
		}
	
		/*** possibly should be error condition?? **/
		if ( np > n)
		{
			throw new nec_exception("ERROR: NP > N IN CONECT()" );
		}
	
		if ( (np == n) && (mp == m) )
			m_ipsym=0;
	} /* if ( ignd != 0) */
	
	if ( n != 0)
	{
		/* Allocate memory to connections */
		icon1.resize((n+m));
		icon2.resize((n+m));
	
		for( i = 0; i < n; i++ )
		{
			iz = i+1;
			xi1= x[i];
			yi1= y[i];
			zi1= z[i];
			xi2= x2[i];
			yi2= y2[i];
			zi2= z2[i];
			slen= sqrt( (xi2- xi1)*(xi2- xi1) + (yi2- yi1) *
				(yi2- yi1) + (zi2- zi1)*(zi2- zi1) ) * SMIN;
		
			/* determine connection data for end 1 of segment. */
			jump = false;
			if ( ignd > 0)
			{
				if ( zi1 <= -slen)
				{
					m_output->nec_printf(
						"\n  GEOMETRY DATA ERROR -- SEGMENT"
						" %d EXTENDS BELOW GROUND", iz );
					stop(-1);
				}
			
				if ( zi1 <= slen)
				{
					icon1[i]= iz;
					z[i]=0.;
					jump = true;	
				} /* if ( zi1 <= slen) */
			} /* if ( ignd > 0) */
		
			if ( ! jump )
			{
				ic= i;
				for( j = 1; j < n; j++)
				{
					ic++;
					if ( ic >= n)
						ic=0;
				
					sep= fabs( xi1- x[ic])+ fabs(yi1- y[ic])+ fabs(zi1- z[ic]);
					if ( sep <= slen)
					{
						icon1[i]= -(ic+1);
						break;
					}
				
					sep= fabs( xi1- x2[ic])+ fabs(yi1- y2[ic])+ fabs(zi1- z2[ic]);
					if ( sep <= slen)
					{
						icon1[i]= (ic+1);
						break;
					}
				
				} /* for( j = 1; j < n; j++) */
			
				if ( ((iz > 0) || (icon1[i] <= PCHCON)) && (sep > slen) )
					icon1[i]=0;
			
			} /* if ( ! jump ) */
		
			/* determine connection data for end 2 of segment. */
			if ( (ignd > 0) || jump )
			{
				if ( zi2 <= -slen)
				{
					m_output->nec_printf(
						"\n  GEOMETRY DATA ERROR -- SEGMENT"
						" %d EXTENDS BELOW GROUND", iz );
					stop(-1);
				}
			
				if ( zi2 <= slen)
				{
					if ( icon1[i] == iz )
					{
						m_output->nec_printf(
							"\n  GEOMETRY DATA ERROR -- SEGMENT"
							" %d LIES IN GROUND PLANE", iz );
						stop(-1);
					}
				
					icon2[i]= iz;
					z2[i]=0.;
					continue;
				} /* if ( zi2 <= slen) */	
			} /* if ( ignd > 0) */
		
			ic= i;
			for( j = 1; j < n; j++ )
			{
				ic++;
				if ( ic >= n)
				ic=0;
			
				sep= fabs(xi2- x[ic])+ fabs(yi2- y[ic])+ fabs(zi2- z[ic]);
				if ( sep <= slen)
				{
					icon2[i]= (ic+1);
					break;
				}
			
				sep= fabs(xi2- x2[ic])+ fabs(yi2- y2[ic])+ fabs(zi2- z2[ic]);
				if ( sep <= slen)
				{
					icon2[i]= -(ic+1);
					break;
				}
			} /* for( j = 1; j < n; j++ ) */
		
			if ( ((iz > 0) || (icon2[i] <= PCHCON)) && (sep > slen) )
			icon2[i]=0;
		
		} /* for( i = 0; i < n; i++ ) */
		
		/* find wire-surface connections for new patches */
		if ( m != 0)
		{
			ix = -1;
			i = 0;
			while( ++i <= m )
			{
				ix++;
				xs= px[ix];
				ys= py[ix];
				zs= pz[ix];
			
				for( iseg = 0; iseg < n; iseg++ )
				{
					xi1= x[iseg];
					yi1= y[iseg];
					zi1= z[iseg];
					xi2= x2[iseg];
					yi2= y2[iseg];
					zi2= z2[iseg];
				
					/* for first end of segment */
					slen=( fabs(xi2- xi1)+ fabs(yi2- yi1)+ fabs(zi2- zi1))* SMIN;
					sep= fabs(xi1- xs)+ fabs(yi1- ys)+ fabs(zi1- zs);
				
					/* connection - divide patch into 4 patches at present array loc. */
					if ( sep <= slen)
					{
						icon1[iseg]=PCHCON+ i;
						ic=0;
						subph( i, ic );
						break;
					}
				
					sep= fabs(xi2- xs)+ fabs(yi2- ys)+ fabs(zi2- zs);
					if ( sep <= slen)
					{
						icon2[iseg]=PCHCON+ i;
						ic=0;
						subph( i, ic );
						break;
					}
				
				} /* for( iseg = 0; iseg < n; iseg++ ) */
			
			} /* while( ++i <= m ) */
		
		} /* if ( m != 0) */
	
	} /* if ( n != 0) */
	
	m_output->nec_printf( "\n\n"
		"     TOTAL SEGMENTS USED: %d   SEGMENTS IN A"
		" SYMMETRIC CELL: %d   SYMMETRY FLAG: %d",
		n, np, m_ipsym );
	
	if ( m > 0)
		m_output->nec_printf(	"\n"
		"       TOTAL PATCHES USED: %d   PATCHES"
		" IN A SYMMETRIC CELL: %d",  m, mp );
	
	
	if (0 == np + mp)
		throw new nec_exception("connect_segments Geometry has zero wires and zero patches.");
	
	iseg=( n+ m)/( np+ mp);
	if ( iseg != 1)
	{
		/*** may be error condition?? ***/
		if ( m_ipsym == 0 )
		{
			nec_error_mode nem(*m_output);
			m_output->endl();
			m_output->line("ERROR: IPSYM=0 IN connect_segments()" );
			stop(-1);
		}
	
		if ( m_ipsym < 0 )
			m_output->nec_printf(
				"\n  STRUCTURE HAS %d FOLD ROTATIONAL SYMMETRY\n", iseg );
		else
		{
			ic= iseg/2;
			if ( iseg == 8)
				ic=3;
			m_output->nec_printf(
				"\n  STRUCTURE HAS %d PLANES OF SYMMETRY\n", ic );
		} /* if ( m_ipsym < 0 ) */
	
	} /* if ( iseg == 1) */
	
	if ( n == 0)
		return;
	
	/* Allocate to connection buffers */
	jco.resize(maxcon);
	
	/* adjust connected seg. ends to exactly coincide.  print junctions */
	/* of 3 or more seg.  also find old seg. connecting to new seg. */
	iseg = 0;
	ipf = false;
	for( j = 0; j < n; j++ )
	{
		jx = j+1;
		iend=-1;
		jend=-1;
		ix= icon1[j];
		ic=1;
		jco[0]= -jx;
		xa= x[j];
		ya= y[j];
		za= z[j];
	
		while( true )
		{
		if ( (ix != 0) && (ix != (j+1)) && (ix <= PCHCON) )
		{
		nsflg=0;
	
		do
		{
			if ( ix == 0 )
			{
				nec_error_mode nem(*m_output);
				m_output->endl();
				m_output->line("\n  CONNECT - SEGMENT CONNECTION ERROR FOR SEGMENT: ");
				m_output->integer(ix);
				stop(-1);
			}
		
			if ( ix < 0 )
				ix= -ix;
			else
				jend= -jend;
		
			jump = false;
		
			if ( ix == jx )
				break;
		
			if ( ix < jx )
			{
				jump = true;
				break;
			}
		
			/* Record max. no. of connections */
			ic++;
			if ( ic >= maxcon )
			{
				maxcon = ic+1;
				jco.resize(maxcon);
			}
			jco[ic-1]= ix* jend;
		
			if ( ix > 0)
				nsflg=1;
		
			ixx = ix-1;
			if ( jend != 1)
			{
				xa= xa+ x[ixx];
				ya= ya+ y[ixx];
				za= za+ z[ixx];
				ix= icon1[ixx];
				continue;
			}
		
			xa= xa+ x2[ixx];
			ya= ya+ y2[ixx];
			za= za+ z2[ixx];
			ix= icon2[ixx];
		} /* do */
		while( ix != 0 );
	
		if ( jump && (iend == 1) )
			break;
		else
			if ( jump )
			{
				iend=1;
				jend=1;
				ix= icon2[j];
				ic=1;
				jco[0]= jx;
				xa= x2[j];
				ya= y2[j];
				za= z2[j];
				continue;
			}
	
		sep= (nec_float)ic;
		xa= xa/ sep;
		ya= ya/ sep;
		za= za/ sep;
	
		for( i = 0; i < ic; i++ )
		{
			ix= jco[i];
			if ( ix <= 0)
			{
				ix=- ix;
				ixx = ix-1;
				x[ixx]= xa;
				y[ixx]= ya;
				z[ixx]= za;
				continue;
			}
		
			ixx = ix-1;
			x2[ixx]= xa;
			y2[ixx]= ya;
			z2[ixx]= za;	
		} /* for( i = 0; i < ic; i++ ) */
	
		if ( ic >= 3)
		{
			if ( ! ipf )
			{
				m_output->nec_printf( "\n\n"
					"    ---------- MULTIPLE WIRE JUNCTIONS ----------\n"
					"    JUNCTION  SEGMENTS (- FOR END 1, + FOR END 2)" );
				ipf = true;
			}
		
			iseg++;
			m_output->nec_printf( "\n   %5d      ", iseg );
		
			for( i = 1; i <= ic; i++ )
			{
				m_output->nec_printf( "%5d", jco[i-1] );
				if ( !(i % 20) )
					m_output->nec_printf( "\n              " );
			}
		} /* if ( ic >= 3) */
	
		} /*if ( (ix != 0) && (ix != j) && (ix <= PCHCON) ) */
	
		if ( iend == 1)
			break;
	
		iend=1;
		jend=1;
		ix= icon2[j];
		ic=1;
		jco[0]= jx;
		xa= x2[j];
		ya= y2[j];
		za= z2[j];
	
		} /* while( true ) */
	
	} /* for( j = 0; j < n; j++ ) */
	
	ax.resize(maxcon);
	bx.resize(maxcon);
	cx.resize(maxcon);
}

/* arc generates segment geometry data for an arc of segment_count segments */
void c_geometry::arc( int tag_id, int segment_count, nec_float rada,
    nec_float ang1, nec_float ang2, nec_float rad )
{
	int istart = n;
	n += segment_count;
	np= n;
	mp= m;
	m_ipsym=0;
	
	if ( segment_count < 1)
		return;
	
	if ( fabs(ang2 - ang1) > 360.0)
	{
		m_output->nec_printf( "\n  ERROR -- ARC ANGLE EXCEEDS 360 DEGREES");
		stop(-1);
	}
	
	/* Reallocate tags buffer */
	segment_tags.resize(n+m);

	/* Reallocate wire buffers */
	x.resize(n);
	y.resize(n);
	z.resize(n);
	x2.resize(n);
	y2.resize(n);
	z2.resize(n);
	bi.resize(n);

	nec_float ang = ang1 * TA;
	nec_float dang = (ang2- ang1)* TA/ segment_count;
	nec_float xs1= rada * cos(ang);
	nec_float zs1= rada * sin(ang);

	for(int i = istart; i < n; i++ )
	{
		ang += dang;
		nec_float xs2 = rada * cos(ang);
		nec_float zs2 = rada * sin(ang);
		x[i]= xs1;
		y[i]=0.;
		z[i]= zs1;
		x2[i]= xs2;
		y2[i]=0.;
		z2[i]= zs2;
		xs1= xs2;
		zs1= zs2;
		bi[i]= rad;
		segment_tags[i]= tag_id;
	} /* for( i = ist; i < n; i++ ) */
}


/*-----------------------------------------------------------------------*/

/* patch generates and modifies patch geometry data */
void c_geometry::patch( int nx, int ny,
    nec_float ax1, nec_float ay1, nec_float az1,
    nec_float ax2, nec_float ay2, nec_float az2,
    nec_float ax3, nec_float ay3, nec_float az3,
    nec_float ax4, nec_float ay4, nec_float az4 )
{
	int mi, ntp, iy, ix;
	nec_float s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., xst=0.;
	nec_float znv, xnv, ynv, xa, xn2, yn2, zn2, salpn, xs, ys, zs, xt, yt, zt;
	
	/* new patches.  for nx=0, ny=1,2,3,4 patch is (respectively) */;
	/* arbitrary, rectagular, triangular, or quadrilateral. */
	/* for nx and ny  > 0 a rectangular surface is produced with */
	/* nx by ny rectangular patches. */
	
	m++;
	mi= m-1;
	
	/* Reallocate patch buffers */
	px.resize(m);
	py.resize(m);
	pz.resize(m);
	t1x.resize(m);
	t1y.resize(m);
	t1z.resize(m);
	t2x.resize(m);
	t2y.resize(m);
	t2z.resize(m);
	pbi.resize(m);
	psalp.resize(m);
	
	if ( nx > 0)
		ntp=2;
	else
		ntp= ny;
	
	if ( ntp <= 1)
	{
		px[mi]= ax1;
		py[mi]= ay1;
		pz[mi]= az1;
		pbi[mi]= az2;
		znv= cos( ax2);
		xnv= znv* cos( ay2);
		ynv= znv* sin( ay2);
		znv= sin( ax2);
		xa= sqrt( xnv* xnv+ ynv* ynv);
	
		if ( xa >= 1.0e-6)
		{
			t1x[mi]=- ynv/ xa;
			t1y[mi]= xnv/ xa;
			t1z[mi]=0.;
		}
		else
		{
			t1x[mi]=1.;
			t1y[mi]=0.;
			t1z[mi]=0.;
		}
	
	} /* if ( ntp <= 1) */
	else
	{
		s1x= ax2- ax1;
		s1y= ay2- ay1;
		s1z= az2- az1;
		s2x= ax3- ax2;
		s2y= ay3- ay2;
		s2z= az3- az2;
	
		if ( nx != 0)
		{
			s1x= s1x/ nx;
			s1y= s1y/ nx;
			s1z= s1z/ nx;
			s2x= s2x/ ny;
			s2y= s2y/ ny;
			s2z= s2z/ ny;
		}
	
		xnv= s1y* s2z- s1z* s2y;
		ynv= s1z* s2x- s1x* s2z;
		znv= s1x* s2y- s1y* s2x;
		xa= sqrt( xnv* xnv+ ynv* ynv+ znv* znv);
		xnv= xnv/ xa;
		ynv= ynv/ xa;
		znv= znv/ xa;
		xst= sqrt( s1x* s1x+ s1y* s1y+ s1z* s1z);
		t1x[mi]= s1x/ xst;
		t1y[mi]= s1y/ xst;
		t1z[mi]= s1z/ xst;
	
		if ( ntp <= 2)
		{
			px[mi]= ax1+.5*( s1x+ s2x);
			py[mi]= ay1+.5*( s1y+ s2y);
			pz[mi]= az1+.5*( s1z+ s2z);
			pbi[mi]= xa;
		}
		else
		{
			if ( ntp != 4)
			{
				px[mi]=( ax1+ ax2+ ax3)/3.;
				py[mi]=( ay1+ ay2+ ay3)/3.;
				pz[mi]=( az1+ az2+ az3)/3.;
				pbi[mi]=.5* xa;
			}
			else
			{
				s1x= ax3- ax1;
				s1y= ay3- ay1;
				s1z= az3- az1;
				s2x= ax4- ax1;
				s2y= ay4- ay1;
				s2z= az4- az1;
				xn2= s1y* s2z- s1z* s2y;
				yn2= s1z* s2x- s1x* s2z;
				zn2= s1x* s2y- s1y* s2x;
				xst= sqrt( xn2* xn2+ yn2* yn2+ zn2* zn2);
				salpn=1./(3.*( xa+ xst));
				px[mi]=( xa*( ax1+ ax2+ ax3)+ xst*( ax1+ ax3+ ax4))* salpn;
				py[mi]=( xa*( ay1+ ay2+ ay3)+ xst*( ay1+ ay3+ ay4))* salpn;
				pz[mi]=( xa*( az1+ az2+ az3)+ xst*( az1+ az3+ az4))* salpn;
				pbi[mi]=.5*( xa+ xst);
				s1x=( xnv* xn2+ ynv* yn2+ znv* zn2)/ xst;
			
				if ( s1x <= 0.9998)
				{
					m_output->nec_printf(
						"\n  ERROR -- CORNERS OF QUADRILATERAL"
						" PATCH DO NOT LIE IN A PLANE" );
					stop(-1);
				}
			} /* if ( ntp != 4) */
		} /* if ( ntp <= 2) */
	
	} /* if ( ntp <= 1) */
	
	t2x[mi]= ynv* t1z[mi]- znv* t1y[mi];
	t2y[mi]= znv* t1x[mi]- xnv* t1z[mi];
	t2z[mi]= xnv* t1y[mi]- ynv* t1x[mi];
	psalp[mi]=1.;
	
	if ( nx != 0)
	{
		m += nx*ny-1;
	
		/* Reallocate patch buffers */
		px.resize(m);
		py.resize(m);
		pz.resize(m);
		t1x.resize(m);
		t1y.resize(m);
		t1z.resize(m);
		t2x.resize(m);
		t2y.resize(m);
		t2z.resize(m);
		pbi.resize(m);
		psalp.resize(m);
	
		xn2= px[mi]- s1x- s2x;
		yn2= py[mi]- s1y- s2y;
		zn2= pz[mi]- s1z- s2z;
		xs= t1x[mi];
		ys= t1y[mi];
		zs= t1z[mi];
		xt= t2x[mi];
		yt= t2y[mi];
		zt= t2z[mi];
	
		for( iy = 0; iy < ny; iy++ )
		{
			xn2 += s2x;
			yn2 += s2y;
			zn2 += s2z;
		
			for( ix = 1; ix <= nx; ix++ )
			{
				xst= (nec_float)ix;
				px[mi]= xn2+ xst* s1x;
				py[mi]= yn2+ xst* s1y;
				pz[mi]= zn2+ xst* s1z;
				pbi[mi]= xa;
				psalp[mi]=1.;
				t1x[mi]= xs;
				t1y[mi]= ys;
				t1z[mi]= zs;
				t2x[mi]= xt;
				t2y[mi]= yt;
				t2z[mi]= zt;
				mi++;
			} /* for( ix = 0; ix < nx; ix++ ) */
		} /* for( iy = 0; iy < ny; iy++ ) */
	
	} /* if ( nx != 0) */
	
	m_ipsym=0;
	np= n;
	mp= m;
}

/*-----------------------------------------------------------------------*/

/*** this function was an 'entry point' (part of) 'patch()' ***/
void c_geometry::subph( int nx, int ny )
{
	int mia, ix, iy, mi;
	nec_float xs, ys, zs, xa, xst, s1x, s1y, s1z, s2x, s2y, s2z, saln, xt, yt;
	
	/* Shift patches to make room for new ones */
	m += 3;
	
	/* Reallocate patch buffers */
	px.resize(m);
	py.resize(m);
	pz.resize(m);
	t1x.resize(m);
	t1y.resize(m);
	t1z.resize(m);
	t2x.resize(m);
	t2y.resize(m);
	t2z.resize(m);
	pbi.resize(m);
	psalp.resize(m);
	
/*
      IF (NY.GT.0) GO TO 10
      IF (NX.EQ.M) GO TO 10
      NXP=NX+1
      IX=LD-M
      DO 9 IY=NXP,M
      IX=IX+1
      NYP=IX-3
      X(NYP)=X(IX)
      Y(NYP)=Y(IX)
      Z(NYP)=Z(IX)
      BI(NYP)=BI(IX)
      SALP(NYP)=SALP(IX)
      T1X(NYP)=T1X(IX)
      T1Y(NYP)=T1Y(IX)
      T1Z(NYP)=T1Z(IX)
      T2X(NYP)=T2X(IX)
      T2Y(NYP)=T2Y(IX)
9     T2Z(NYP)=T2Z(IX)
10 ...

	if (ny <= 0) && (nx != m)
	{
		nxp = nx + 1;
		int ix = ld - m;
		
		for (int iy = nxp; iy < m; iy++)
		{
			ix += 1;
			int nyp =ix-3;
			px(nyp)=px(ix);
			py(nyp)=py(ix);
			pz(nyp)=pz(ix);
			pbi(nyp)=pbi(ix);
			psalp(nyp)=psalp(ix);
			t1x(nyp)=t1x(ix);
			t1y(nyp)=t1y(ix);
			t1z(nyp)=t1z(ix);
			t2x(nyp)=t2x(ix);
			t2y(nyp)=t2y(ix);
			t2z(nyp)=t2z(ix);
		}
	}
*/
	if ( (ny <= 0) && (nx != m) )
	{
//		for( iy = m-1; iy >= nx; iy-- )
		for( iy = m-1; iy > nx; iy-- )
		{
			px[iy]= px[iy-3];
			py[iy]= py[iy-3];
			pz[iy]= pz[iy-3];
			pbi[iy]= pbi[iy-3];
			psalp[iy]= psalp[iy-3];
			t1x[iy]= t1x[iy-3];
			t1y[iy]= t1y[iy-3];
			t1z[iy]= t1z[iy-3];
			t2x[iy]= t2x[iy-3];
			t2y[iy]= t2y[iy-3];
			t2z[iy]= t2z[iy-3];
		}
	} /* if ( (ny <= 0) || (nx != m) ) */
	
	/* divide patch for connection */
	mi= nx-1;
	xs= px[mi];
	ys= py[mi];
	zs= pz[mi];
	xa= pbi[mi]/4.;
	xst= sqrt( xa)/2.;
	s1x= t1x[mi];
	s1y= t1y[mi];
	s1z= t1z[mi];
	s2x= t2x[mi];
	s2y= t2y[mi];
	s2z= t2z[mi];
	saln= psalp[mi];
	xt= xst;
	yt= xst;
	
	if ( ny == 0)
		mia= mi;
	else
	{
		mp++;
		mia= m-1;
	}
	
	for( ix = 1; ix <= 4; ix++ )
	{
		px[mia]= xs+ xt* s1x+ yt* s2x;
		py[mia]= ys+ xt* s1y+ yt* s2y;
		pz[mia]= zs+ xt* s1z+ yt* s2z;
		pbi[mia]= xa;
		t1x[mia]= s1x;
		t1y[mia]= s1y;
		t1z[mia]= s1z;
		t2x[mia]= s2x;
		t2y[mia]= s2y;
		t2z[mia]= s2z;
		psalp[mia]= saln;
	
		if ( ix == 2)
			yt=- yt;
	
		if ( (ix == 1) || (ix == 3) )
			xt=- xt;
	
		mia++;
	}
	
	if ( nx <= mp)
		mp += 3;
	
	if ( ny > 0 )
		pz[mi]=10000.;
}

/*-----------------------------------------------------------------------

	Read Geometry Data from a Card

-------------------------------------------------------------------------*/

void c_geometry::read_geometry_card(FILE* input_fp,  char *gm, int *i1, int *i2, nec_float *x1, nec_float *y1,
    nec_float *z1, nec_float *x2, nec_float *y2, nec_float *z2, nec_float *rad )
{
	char line_buf[134];
	int i, line_idx;
	int n_integer_params = 2, n_float_params = 7;
	int integer_array[2] = { 0, 0 };
	nec_float real_array[7] = { 0., 0., 0., 0., 0., 0., 0. };
	
	/* read a line from input file */
	load_line( line_buf, input_fp );
	
	/* get line length */
	int line_length= strlen( line_buf );
	
	/* abort if card's mnemonic too short or missing */
	if ( line_length < 2 )
	{
		nec_error_mode em(*m_output);
		m_output->endl();
		m_output->line("GEOMETRY DATA CARD ERROR:" );
		m_output->line("CARD'S MNEMONIC CODE TOO SHORT OR MISSING." );
		stop(-1);
	}
	
	/* extract card's mnemonic code */
	strncpy( gm, line_buf, 2 );
	gm[2] = '\0';
	
	/* Exit if "XT" command read (for testing) */
	if ( strcmp( gm, "XT" ) == 0 )
	{
			nec_error_mode em(*m_output);
			m_output->endl();
			m_output->line("Exiting after an \"XT\" command in read_geometry_card()" );
			stop(0);
	}
	
	/* Return if only mnemonic on card */
	if ( line_length == 2 )
	{
		*i1 = *i2 = 0;
		*x1 = *y1 = *z1 = *x2 = *y2 = *z2 = *rad = 0.;
		return;
	}
	
	/* read integers from line */
	line_idx = 1;
	for( i = 0; i < n_integer_params; i++ )
	{
		/* Find first numerical character */
		while( ((line_buf[++line_idx] <  '0')  ||
			(line_buf[  line_idx] >  '9')) &&
			(line_buf[  line_idx] != '+')  &&
			(line_buf[  line_idx] != '-') )
		if ( (line_buf[line_idx] == '\0') )
		{
			*i1= integer_array[0];
			*i2= integer_array[1];
			*x1= real_array[0];
			*y1= real_array[1];
			*z1= real_array[2];
			*x2= real_array[3];
			*y2= real_array[4];
			*z2= real_array[5];
			*rad= real_array[6];
			return;
		}
	
		/* read an integer from line */
		integer_array[i] = atoi( &line_buf[line_idx] );
	
		/* traverse numerical field to next ' ' or ',' or '\0' */
		line_idx--;
		while( 	(line_buf[++line_idx] != ' ') &&
				(line_buf[  line_idx] != ',') &&
				(line_buf[  line_idx] != '\0') )
		{
			/* test for non-numerical characters */
			if ( 	((line_buf[line_idx] <  '0')  ||
					(line_buf[line_idx] >  '9')) &&
					(line_buf[line_idx] != '+')  &&
					(line_buf[line_idx] != '-') )
			{
				m_output->nec_printf(
					"\n  GEOMETRY DATA CARD \"%s\" ERROR:"
					"\n  NON-NUMERICAL CHARACTER '%c' IN INTEGER FIELD AT CHAR. %d\n",
					gm, line_buf[line_idx], (line_idx+1)  );
				stop(-1);
			}
		} /* while( (line_buff[++line_idx] ... */
	
		/* Return on end of line */
		if ( line_buf[line_idx] == '\0' )
		{
			*i1= integer_array[0];
			*i2= integer_array[1];
			*x1= real_array[0];
			*y1= real_array[1];
			*z1= real_array[2];
			*x2= real_array[3];
			*y2= real_array[4];
			*z2= real_array[5];
			*rad= real_array[6];
			return;
		}
	
	} /* for( i = 0; i < n_integer_params; i++ ) */
	
	/* read nec_floats from line */
	for( i = 0; i < n_float_params; i++ )
	{
		/* Find first numerical character */
		while( ((line_buf[++line_idx] <  '0')  ||
			(line_buf[  line_idx] >  '9')) &&
			(line_buf[  line_idx] != '+')  &&
			(line_buf[  line_idx] != '-')  &&
			(line_buf[  line_idx] != '.') )
		if ( (line_buf[line_idx] == '\0') )
		{
			*i1= integer_array[0];
			*i2= integer_array[1];
			*x1= real_array[0];
			*y1= real_array[1];
			*z1= real_array[2];
			*x2= real_array[3];
			*y2= real_array[4];
			*z2= real_array[5];
			*rad= real_array[6];
			return;
		}
	
		/* read a nec_float from line */
		real_array[i] = atof( &line_buf[line_idx] );
	
		/* traverse numerical field to next ' ' or ',' or '\0' */
		line_idx--;
		while(	(line_buf[++line_idx] != ' ') &&
				(line_buf[  line_idx] != ',') &&
				(line_buf[  line_idx] != '\0') )
		{
			/* test for non-numerical characters */
			if ( ((line_buf[line_idx] <  '0')  ||
				(line_buf[line_idx] >  '9')) &&
				(line_buf[line_idx] != '.')  &&
				(line_buf[line_idx] != '+')  &&
				(line_buf[line_idx] != '-')  &&
				(line_buf[line_idx] != 'E')  &&
				(line_buf[line_idx] != 'e') )
			{
				m_output->nec_printf(
					"\n  GEOMETRY DATA CARD \"%s\" ERROR:"
					"\n  NON-NUMERICAL CHARACTER '%c' IN FLOAT FIELD AT CHAR. %d.\n",
					gm, line_buf[line_idx], (line_idx+1) );
				stop(-1);
			}
		} /* while( (line_buff[++line_idx] ... */
	
		/* Return on end of line */
		if ( line_buf[line_idx] == '\0' )
		{
			*i1= integer_array[0];
			*i2= integer_array[1];
			*x1= real_array[0];
			*y1= real_array[1];
			*z1= real_array[2];
			*x2= real_array[3];
			*y2= real_array[4];
			*z2= real_array[5];
			*rad= real_array[6];
			return;
		}
	
	} /* for( i = 0; i < n_float_params; i++ ) */
	
	*i1  = integer_array[0];
	*i2  = integer_array[1];
	*x1  = real_array[0];
	*y1  = real_array[1];
	*z1  = real_array[2];
	*x2  = real_array[3];
	*y2  = real_array[4];
	*z2  = real_array[5];
	*rad = real_array[6];
}



/* compute basis function i */
void c_geometry::tbf( int i, int icap )
{
	int ix, jcox, jcoxx, jend, iend, njun1=0, njun2, jsnop, jsnox;
	nec_float pp, sdh, cdh, sd, omc, aj, pm=0, cd, ap, qp, qm, xxi;
	nec_float d, sig; /*** also global ***/
	
	jsno=0;
	pp=0.;
	ix = i-1;
	jcox= icon1[ix];
	
	if ( jcox > PCHCON)
		jcox= i;
	
	jend=-1;
	iend=-1;
	sig=-1.;
	
	do
	{
		if ( jcox != 0 )
		{
			if ( jcox < 0 )
				jcox=- jcox;
			else
			{
				sig=- sig;
				jend=- jend;
			}
		
			jcoxx = jcox-1;
			jsno++;
			jsnox = jsno-1;
			jco[jsnox]= jcox;
			d= pi()* segment_length[jcoxx];
			sdh= sin( d);
			cdh= cos( d);
			sd=2.* sdh* cdh;
		
			if ( d <= 0.015)
			{
				omc=4.* d* d;
				omc=((1.3888889e-3* omc-4.1666666667e-2)* omc+.5)* omc;
			}
			else
				omc=1.- cdh* cdh+ sdh* sdh;
		
			aj=1./( log(1./( pi()* bi[jcoxx]))-.577215664);
			pp= pp- omc/ sd* aj;
			ax[jsnox]= aj/ sd* sig;
			bx[jsnox]= aj/(2.* cdh);
			cx[jsnox]=- aj/(2.* sdh)* sig;
		
			if ( jcox != i)
			{
				if ( jend == 1)
					jcox= icon2[jcoxx];
				else
					jcox= icon1[jcoxx];
			
				if ( abs(jcox) != i )
				{
					if ( jcox != 0 )
						continue;
					else
					{
						m_output->nec_printf(
							"\n  TBF - SEGMENT CONNECTION ERROR FOR SEGMENT %5d", i );
						stop(-1);
					}
				}
			} /* if ( jcox != i) */
			else
			bx[jsnox] =- bx[jsnox];
		
			if ( iend == 1)
				break;
		} /* if ( jcox != 0 ) */
	
		pm=- pp;
		pp=0.;
		njun1= jsno;
	
		jcox= icon2[ix];
		if ( jcox > PCHCON)
			jcox= i;
	
		jend=1;
		iend=1;
		sig=-1.;
	
	} /* do */
	while( jcox != 0 );
	
	njun2= jsno- njun1;
	jsnop= jsno;
	jco[jsnop]= i;
	d= pi()* segment_length[ix];
	sdh= sin( d);
	cdh= cos( d);
	sd=2.* sdh* cdh;
	cd= cdh* cdh- sdh* sdh;
	
	if ( d <= 0.015)
	{
		omc=4.* d* d;
		omc=((1.3888889e-3* omc-4.1666666667e-2)* omc+.5)* omc;
	}
	else
		omc=1.- cd;
	
	ap=1./( log(1./( pi()* bi[ix]))-.577215664);
	aj= ap;
	
	if ( njun1 == 0)
	{
		if ( njun2 == 0)
		{
			bx[jsnop]=0.;
		
			if ( icap == 0)
				xxi=0.;
			else
			{
				qp= pi()* bi[ix];
				xxi= qp* qp;
				xxi= qp*(1.-.5* xxi)/(1.- xxi);
			}
		
			cx[jsnop]=1./( cdh- xxi* sdh);
			jsno= jsnop+1;
			ax[jsnop]=-1.;
			return;
		} /* if ( njun2 == 0) */
	
		if ( icap == 0)
			xxi=0.;
		else
		{
			qp= pi()* bi[ix];
			xxi= qp* qp;
			xxi= qp*(1.-.5* xxi)/(1.- xxi);
		}
	
		qp=-( omc+ xxi* sd)/( sd*( ap+ xxi* pp)+ cd*( xxi* ap- pp));
		d= cd- xxi* sd;
		bx[jsnop]=( sdh+ ap* qp*( cdh- xxi* sdh))/ d;
		cx[jsnop]=( cdh+ ap* qp*( sdh+ xxi* cdh))/ d;
	
		for( iend = 0; iend < njun2; iend++ )
		{
			ax[iend]=-ax[iend]* qp;
			bx[iend]= bx[iend]* qp;
			cx[iend]=- cx[iend]* qp;
		}
	
		jsno= jsnop+1;
		ax[jsnop]=-1.;
		return;
	} /* if ( njun1 == 0) */
	
	if ( njun2 == 0)
	{
		if ( icap == 0)
			xxi=0.;
		else
		{
			qm= pi()* bi[ix];
			xxi= qm* qm;
			xxi= qm*(1.-.5* xxi)/(1.- xxi);
		}
	
		qm=( omc+ xxi* sd)/( sd*( aj- xxi* pm)+ cd*( pm+ xxi* aj));
		d= cd- xxi* sd;
		bx[jsnop]=( aj* qm*( cdh- xxi* sdh)- sdh)/ d;
		cx[jsnop]=( cdh- aj* qm*( sdh+ xxi* cdh))/ d;
	
		for( iend = 0; iend < njun1; iend++ )
		{
			ax[iend]= ax[iend]* qm;
			bx[iend]= bx[iend]* qm;
			cx[iend]= cx[iend]* qm;
		}
	
		jsno= jsnop+1;
		ax[jsnop]=-1.;
		return;
	
	} /* if ( njun2 == 0) */
	
	qp= sd*( pm* pp+ aj* ap)+ cd*( pm* ap- pp* aj);
	qm=( ap* omc- pp* sd)/ qp;
	qp=-( aj* omc+ pm* sd)/ qp;
	bx[jsnop]=( aj* qm+ ap* qp)* sdh/ sd;
	cx[jsnop]=( aj* qm- ap* qp)* cdh/ sd;
	
	for( iend = 0; iend < njun1; iend++ )
	{
		ax[iend]= ax[iend]* qm;
		bx[iend]= bx[iend]* qm;
		cx[iend]= cx[iend]* qm;
	}
	
	jend= njun1;
	for( iend = jend; iend < jsno; iend++ )
	{
		ax[iend]=- ax[iend]* qp;
		bx[iend]= bx[iend]* qp;
		cx[iend]=- cx[iend]* qp;
	}
	
	jsno= jsnop+1;
	ax[jsnop]=-1.;
}

/*-----------------------------------------------------------------------*/

/*-----------------------------------------------------------------------*/

/* compute the components of all basis functions on segment j */
void c_geometry::trio( int j )
{
	int jcox, jcoxx, jsnox, jx, jend=0, iend=0;
	
	jsno=0;
	jx = j-1;
	jcox= icon1[jx];
	jcoxx = jcox-1;
	
	if ( jcox <= PCHCON)
	{
		jend=-1;
		iend=-1;
	}
	
	if ( (jcox == 0) || (jcox > PCHCON) )
	{
		jcox= icon2[jx];
		jcoxx = jcox-1;
	
		if ( jcox <= PCHCON)
		{
			jend=1;
			iend=1;
		}
	
		if ( jcox == 0 || (jcox > PCHCON) )
		{
			jsnox = jsno;
			jsno++;
		
			/* Allocate to connections buffers */
			if ( jsno >= maxcon )
			{
				maxcon = jsno +1;
				jco.resize(maxcon);
				ax.resize(maxcon);
				bx.resize(maxcon);
				cx.resize(maxcon);
			}
		
			sbf( j, j, &ax[jsnox], &bx[jsnox], &cx[jsnox]);
			jco[jsnox]= j;
			return;
		}
	} /* if ( (jcox == 0) || (jcox > PCHCON) ) */
	
	do
	{
		if ( jcox < 0 )
			jcox=- jcox;
		else
			jend=- jend;
		jcoxx = jcox-1;
	
		if ( jcox != j)
		{
			jsnox = jsno;
			jsno++;
		
			/* Allocate to connections buffers */
			if ( jsno >= maxcon )
			{
				maxcon = jsno +1;
				jco.resize(maxcon );
				ax.resize(maxcon);
				bx.resize(maxcon);
				cx.resize(maxcon);
			}
		
			sbf( jcox, j, &ax[jsnox], &bx[jsnox], &cx[jsnox]);
			jco[jsnox]= jcox;
		
			if ( jend != 1)
				jcox= icon1[jcoxx];
			else
				jcox= icon2[jcoxx];
		
			if ( jcox == 0 )
			{
				m_output->nec_printf(
					"\n  TRIO - SEGMENT CONNENTION ERROR FOR SEGMENT %5d", j );
				stop(-1);
			}
			else
				continue;
		} /* if ( jcox != j) */
	
		if ( iend == 1)
			break;
	
		jcox= icon2[jx];
	
		if ( jcox > PCHCON)
			break;
	
		jend=1;
		iend=1;
	} /* do */
	while( jcox != 0 );
	
	jsnox = jsno;
	jsno++;
	
	/* Allocate to connections buffers */
	if ( jsno >= maxcon )
	{
		maxcon = jsno +1;
		jco.resize(maxcon );
		ax.resize(maxcon);
		bx.resize(maxcon);
		cx.resize(maxcon);
	}
	
	sbf( j, j, &ax[jsnox], &bx[jsnox], &cx[jsnox]);
	jco[jsnox]= j;
}



/* compute component of basis function i on segment is. */
void c_geometry::sbf( int i, int is, nec_float *aa, nec_float *bb, nec_float *cc )
{
	int ix, jsno, june, jcox, jcoxx, jend, iend, njun1=0, njun2;
	nec_float d, sig, pp, sdh, cdh, sd, omc, aj, pm=0, cd, ap, qp, qm, xxi;
	
	*aa=0.;
	*bb=0.;
	*cc=0.;
	june=0;
	jsno=0;
	pp=0.;
	ix=i-1;
	
	jcox= icon1[ix];
	if ( jcox > PCHCON)
		jcox= i;
	jcoxx = jcox-1;
	
	jend=-1;
	iend=-1;
	sig=-1.;
	
	do
	{
		if ( jcox != 0 )
		{
			if ( jcox < 0 )
				jcox=- jcox;
			else
			{
				sig=- sig;
				jend=- jend;
			}
		
			jcoxx = jcox-1;
			jsno++;
			d= pi()* segment_length[jcoxx];
			sdh= sin( d);
			cdh= cos( d);
			sd=2.* sdh* cdh;
		
			if ( d <= 0.015)
			{
				omc=4.* d* d;
				omc=((1.3888889e-3* omc -4.1666666667e-2)* omc +.5)* omc;
			}
			else
				omc=1.- cdh* cdh+ sdh* sdh;
		
			aj=1./( log(1./( pi()* bi[jcoxx]))-.577215664);
			pp -= omc/ sd* aj;
		
			if ( jcox == is)
			{
				*aa= aj/ sd* sig;
				*bb= aj/(2.* cdh);
				*cc=- aj/(2.* sdh)* sig;
				june= iend;
			}
		
			if ( jcox != i )
			{
				if ( jend != 1)
					jcox= icon1[jcoxx];
				else
					jcox= icon2[jcoxx];
			
				if ( abs(jcox) != i )
				{
					if ( jcox == 0 )
					{
						m_output->nec_printf(
						"\n  SBF - SEGMENT CONNECTION ERROR FOR SEGMENT %d", i);
						stop(-1);
					}
					else
						continue;
				}	
			} /* if ( jcox != i ) */
			else
			if ( jcox == is)
				*bb=- *bb;
		
			if ( iend == 1)
				break;
		} /* if ( jcox != 0 ) */
	
		pm=- pp;
		pp=0.;
		njun1= jsno;
	
		jcox= icon2[ix];
		if ( jcox > PCHCON)
			jcox= i;
	
		jend=1;
		iend=1;
		sig=-1.;
	
	} /* do */
	while( jcox != 0 );
	
	njun2= jsno- njun1;
	d= pi()* segment_length[ix];
	sdh= sin( d);
	cdh= cos( d);
	sd=2.* sdh* cdh;
	cd= cdh* cdh- sdh* sdh;
	
	if ( d <= 0.015)
	{
		omc=4.* d* d;
		omc=((1.3888889e-3* omc -4.1666666667e-2)* omc +.5)* omc;
	}
	else
		omc=1.- cd;
	
	ap=1./( log(1./( pi()* bi[ix])) -.577215664);
	aj= ap;
	
	if ( njun1 == 0)
	{
		if ( njun2 == 0)
		{
			*aa =-1.;
			qp= pi()* bi[ix];
			xxi= qp* qp;
			xxi= qp*(1.-.5* xxi)/(1.- xxi);
			*cc=1./( cdh- xxi* sdh);
				return;
		}
	
		qp= pi()* bi[ix];
		xxi= qp* qp;
		xxi= qp*(1.-.5* xxi)/(1.- xxi);
		qp=-( omc+ xxi* sd)/( sd*( ap+ xxi* pp)+ cd*( xxi* ap- pp));
	
		if ( june == 1)
		{
			*aa=- *aa* qp;
			*bb=  *bb* qp;
			*cc=- *cc* qp;
			if ( i != is)
				return;
		}
	
		*aa -= 1.;
		d = cd - xxi * sd;
		*bb += (sdh + ap * qp * (cdh - xxi * sdh)) / d;
		*cc += (cdh + ap * qp * (sdh + xxi * cdh)) / d;
		return;
	
	} /* if ( njun1 == 0) */
	
	if ( njun2 == 0)
	{
		qm= pi()* bi[ix];
		xxi= qm* qm;
		xxi= qm*(1.-.5* xxi)/(1.- xxi);
		qm=( omc+ xxi* sd)/( sd*( aj- xxi* pm)+ cd*( pm+ xxi* aj));
	
		if ( june == -1)
		{
			*aa= *aa* qm;
			*bb= *bb* qm;
			*cc= *cc* qm;
			if ( i != is)
				return;
		}
	
		*aa -= 1.;
		d= cd- xxi* sd;
		*bb += ( aj* qm*( cdh- xxi* sdh)- sdh)/ d;
		*cc += ( cdh- aj* qm*( sdh+ xxi* cdh))/ d;
		return;	
	} /* if ( njun2 == 0) */
	
	qp= sd*( pm* pp+ aj* ap)+ cd*( pm* ap- pp* aj);
	qm=( ap* omc- pp* sd)/ qp;
	qp=-( aj* omc+ pm* sd)/ qp;
	
	if ( june != 0 )
	{
		if ( june < 0 )
		{
			*aa= *aa* qm;
			*bb= *bb* qm;
			*cc= *cc* qm;
		}
		else
		{
			*aa=- *aa* qp;
			*bb= *bb* qp;
			*cc=- *cc* qp;
		}
	
		if ( i != is)
			return;
	} /* if ( june != 0 ) */
	
	*aa -= 1.;
	*bb += ( aj* qm+ ap* qp)* sdh/ sd;
	*cc += ( aj* qm- ap* qp)* cdh/ sd;
}


/*
	get_current_coefficients computes coefficients of the:
		constant 	[air, aii]
		sine		[bir, bii]
		cosine		[cir, cii]
	terms in the current interpolation functions for the current vector curx.
 */
void c_geometry::get_current_coefficients(nec_float wavelength, complex_array& curx,
			real_array& air, real_array& aii,
			real_array& bir, real_array& bii,
			real_array& cir, real_array& cii,
			complex_array& vqds, int nqds,
			int_array& iqds)
{
	static nec_complex s_CCJ(0.0,-0.01666666667);
	
	nec_float ar, ai, sh;
	nec_complex cs1, cs2;
	
	if ( n != 0)
	{
		for(int i = 0; i < n; i++ )
		{
			air[i]=0.;
			aii[i]=0.;
			bir[i]=0.;
			bii[i]=0.;
			cir[i]=0.;
			cii[i]=0.;
		}
	
		for (int i = 0; i < n; i++ )
		{
			ar= real( curx[i]);
			ai= imag( curx[i]);
			tbf( i+1, 1 );
		
			for (int jx = 0; jx < jsno; jx++ )
			{
				int j = jco[jx]-1;
				air[j] += ax[jx]* ar;
				aii[j] += ax[jx]* ai;
				bir[j] += bx[jx]* ar;
				bii[j] += bx[jx]* ai;
				cir[j] += cx[jx]* ar;
				cii[j] += cx[jx]* ai;
			}
		} /* for( i = 0; i < n; i++ ) */
	
		for (int is = 0; is < nqds; is++ )
		{
			int i= iqds[is]-1;
			int jx= icon1[i];
			icon1[i]=0;
			tbf(i+1,0);
			icon1[i]= jx;
			sh = segment_length[i]*.5;
			
			
			nec_complex curd = s_CCJ * 
				vqds[is]/( (log(2.* sh/ bi[i])-1.)*
				(bx[jsno-1]* cos(two_pi() * sh)+ cx[jsno-1]* sin(two_pi() * sh))* wavelength );
			ar = real( curd);
			ai = imag( curd);
		
			for ( jx = 0; jx < jsno; jx++ )
			{
				int j = jco[jx]-1;
				air[j] += ax[jx]* ar;
				aii[j] += ax[jx]* ai;
				bir[j] += bx[jx]* ar;
				bii[j] += bx[jx]* ai;
				cir[j] += cx[jx]* ar;
				cii[j] += cx[jx]* ai;
			}
		
		} /* for( is = 0; is < nqds; is++ ) */
	
		for (int i = 0; i < n; i++ )
			curx[i]= nec_complex( air[i]+cir[i], aii[i]+cii[i] );
	
	} /* if ( n != 0) */
	
	if ( m == 0)
		return;
	
	/* convert surface currents from */
	/* t1,t2 components to x,y,z components */
	int jco1 = n_plus_2m;
	int jco2 = jco1 + m;
	
	for (int i = 1; i <= m; i++ )
	{
		jco1 -= 2;
		jco2 -= 3;
		cs1= curx[jco1];
		cs2= curx[jco1+1];
		curx[jco2]  = cs1* t1x[m-i]+ cs2* t2x[m-i];
		curx[jco2+1]= cs1* t1y[m-i]+ cs2* t2y[m-i];
		curx[jco2+2]= cs1* t1z[m-i]+ cs2* t2z[m-i];
	}
}


void c_geometry::frequency_scale(nec_float freq_mhz, real_array& xtemp, real_array& ytemp,
	real_array& ztemp, real_array& sitemp, real_array& bitemp)
{
	nec_float fr = freq_mhz/ CVEL;

	if ( n != 0)
	{
		for (int i = 0; i < n; i++ )
		{
			x[i]= xtemp[i]* fr;
			y[i]= ytemp[i]* fr;
			z[i]= ztemp[i]* fr;
			segment_length[i]= sitemp[i]* fr;
			bi[i]= bitemp[i]* fr;
		}
	}

	if ( m != 0)
	{
		nec_float fr2 = fr* fr;
		for (int i = 0; i < m; i++ )
		{
			int j = i+n;
			px[i]= xtemp[j]* fr;
			py[i]= ytemp[j]* fr;
			pz[i]= ztemp[j]* fr;
			pbi[i]= bitemp[j]* fr2;
		}
	}
}

void c_geometry::fill_temp_geom(int *ifrtmw, int *ifrtmp, real_array& xtemp, 
	real_array& ytemp, real_array& ztemp, real_array& sitemp, 
	real_array& bitemp)
{
	if ( (n != 0) && (*ifrtmw != 1) )
	{
		*ifrtmw=1;
		for (int i = 0; i < n; i++ )
		{
			xtemp[i]= x[i];
			ytemp[i]= y[i];
			ztemp[i]= z[i];
			sitemp[i]= segment_length[i];
			bitemp[i]= bi[i];
		}
	}

	if ( (m != 0) && (*ifrtmp != 1) )
	{
		*ifrtmp=1;
		for (int i = 0; i < m; i++ )
		{
			int j = i+n;
			xtemp[j]= px[i];
			ytemp[j]= py[i];
			ztemp[j]= pz[i];
			bitemp[j]= pbi[i];
		}
	}
}

