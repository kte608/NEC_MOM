/* Translation to C++ by Tim Molteno

	Based on the C port by N. Kyriazis
 	Including pieces from additional work by
		Jeroen Vreeken <pe1rxq@amsat.org>
	
	Fixed a few bugs in the process.
	
	Using the std vector library in preparation
	for moving to ATLAS for the matrix and vector
	operations.

		Debian Build Instructions
		
	apt-get install atlas3-base atlas3-headers atlas3-base-dev
	apt-get install refblas3-dev lapack3-dev lapack3-doc
	
	For more information on using LAPACK for doing 
	efficient computation, see
		http://seehuhn.de/comp/linear.html
*/

/* Original disclaimer that came with the FORTRAN code */

/******* Translated to the C language by N. Kyriazis  20 Aug 2003 *******/
/*									*/
/* Program NEC(input,tape5=input,output,tape11,tape12,tape13,tape14,	*/
/* tape15,tape16,tape20,tape21)						*/
/*									*/
/* Numerical Electromagnetics Code (NEC2)  developed at Lawrence	*/
/* Livermore lab., Livermore, CA.  (contact G. Burke at 415-422-8414	*/
/* for problems with the NEC code. For problems with the vax implem- 	*/
/* entation, contact J. Breakall at 415-422-8196 or E. Domning at 415 	*/
/* 422-5936) 								*/
/* file created 4/11/80. 						*/
/*									*/
/*                ***********Notice********** 				*/
/* This computer code material was prepared as an account of work 	*/
/* sponsored by the United States government.  Neither the United 	*/
/* States nor the United States Department Of Energy, nor any of 	*/
/* their employees, nor any of their contractors, subcontractors, 	*/
/* or their employees, makes any warranty, express or implied, or	*/
/* assumes any legal liability or responsibility for the accuracy, 	*/
/* completeness or usefulness of any information, apparatus, product 	*/
/* or process disclosed, or represents that its use would not infringe 	*/
/* privately-owned rights. 						*/
/*									*/
/************************************************************************/

#include "getopt.h"

#include "nec2cpp.h"

#include <signal.h>
#include <vector>
#include <string>

using namespace std;

#include "nec_context.h"

#ifndef __MINGW32__
/* Signal handler */
static void sig_handler( int signal );
#endif


/*-------------------------------------------------------------------*/

int nec_main( int argc, char **argv, nec_output_file& s_output );

/*	New main() function

	This places an exception handler around the old main loop to 
	allow errors to be nicely caught!
*/
int main( int argc, char **argv )
{
	nec_output_file s_output;
	
	try
	{
		nec_main(argc, argv, s_output);
	}
	catch (const char* message)
	{
		nec_error_mode nem(s_output);
		s_output.line("NEC++ Runtime Error: ");
		s_output.line(message);
		exit(1);
	}
#ifdef NEC_ERROR_CHECK
	catch (nec_exception* nex)
	{
		nec_error_mode nem(s_output);
		s_output.line("NEC++ Runtime Error: ");
		s_output.line(nex->get_message().c_str());
		exit(1);
	}
#endif
}

#include "c_geometry.h"
int benchmark();

int benchmark()
{
/*
CM Center feed Horizontal Half-Wave Dipole over excellent ground
CE
GW 1,299,-139.,0, 6.,+139.,0, 6., .001,
GE 0,
GN 1,
FR 0,0,0,0, 0.54,
EX 0, 1,150,0,1., 0.,
RP 1, 1, 1,0000, 1.5, 0., 0., 0., 1000.,
EN



*/
	cout << "The nec2++ benchmark." << endl << endl;
	
	nec_float start_timer, stop_timer;
	
	secnds( &start_timer );
	for (int i=0; i<20; i++)
	{
		{
			/*
			CMEXAMPLE 2. CENTER FED LINEAR ANTENNA. 
			CM           CURRENT SLOPE DISCONTINUITY SOURCE. 
			CM           1. THIN PERFECTLY CONDUCTING WIRE 
			CE           2. THIN ALUMINUM WIRE 
			GW 0 8 0. 0. -.25 0. 0. .25 .00001 
			GE 
			FR 0 3 0 0 200. 50. 
			EX 5 0 5 1 1. 0. 50. 
			XQ 
			LD 5 0 0 0 3.720E+07 
			FR 0 1 0 0 300. 
			EX 5 0 5 0 1. 
			GN 1
			XQ 
			EN
			*/

			nec_context nec;
			nec.set_results_stdout(false);
			nec.set_gain_only(true);
			nec.initialize();
			
			c_geometry& geo = nec.get_geometry();
			geo.wire(0,8, 0.0, 0.0, -0.25,
					0.0, 0.0, 0.25,
					0.00001,
					1.0, 1.0);
					
			nec.geometry_complete(0,0);
			
			nec.fr_card(0, 3, 200.0, 50.0);
			nec.ex_card(5, 0, 5, 1, 1.0, 0.0, 50.0, 0.0, 0.0, 0.0);
			nec.xq_card(0);
			nec.ld_card(5, 0, 0,0, 3.72e7, 0.0, 0.0);
			nec.fr_card(0, 1, 300.0, 0.0);
			nec.ex_card(5, 0, 5, 0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
			nec.gn_card(1, 0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0);
			
			// nec.xq_card(0);
	/*	This line kills the version on Wine! Very strange! I think that the gn card is bad
			fixme:msvcrt:MSVCRT_signal (11 (nil)):stub
			err:seh:EXC_DefaultHandling Unhandled exception code c0000005 flags 0 addr 0x404e2380
			Wine failed with return code 1
	*/
		}
		{
			/*
				An example showing how to evaluate for maximum gain.
				
			CM GA - NEC FILE
			CE go blue ! 
			GW 0 36 0 0 0 -0.042 0.008 0.017 0.001
			GW 0 21 -0.042 0.008 0.017 -0.048 0.021 -0.005 0.001
			GW 0 70 -0.048 0.021 -0.005 0.039 0.032 -0.017 0.001
			GW 0 70 -0.048 0.021 -0.005 0.035 0.043 0.014 0.001
			GW 0 50 -0.042 0.008 0.017 0.017 -0.015 0.014 0.001
			GW 0 66 0.017 -0.015 0.014 -0.027 0.04 -0.031 0.001
			GW 0 85 -0.027 0.04 -0.031 0.046 -0.01 0.028 0.001
			GW 0 47 0.046 -0.01 0.028 -0.013 -0.005 0.031 0.001
			GW 0 70 0.017 -0.015 0.014 -0.048 -0.038 -0.04 0.001
			GW 0 77 -0.048 -0.038 -0.04 0.049 -0.045 -0.04 0.001
			GE 0
			GN -1
			LD 5 0 0 0 3.720E+07 
			FR 0 1 0 0 2400
			PT -1
			EX 1 1 1 0 0 0 0 0 0 0 0
			RP 0 1 1 0500 90 90 0 0
			EN
			
			*/
			nec_context nec;
			nec.set_results_stdout(false);
			nec.set_gain_only(false);
			nec.initialize();
			
			c_geometry& geo = nec.get_geometry();
			geo.wire(0, 36, 0, 0, 0, -0.042, 0.008, 0.017, 0.001, 1.0, 1.0);
			geo.wire(0, 21, -0.042, 0.008, 0.017, -0.048, 0.021, -0.005, 0.001, 1.0, 1.0);
			geo.wire(0, 70, -0.048, 0.021, -0.005, 0.039, 0.032, -0.017, 0.001, 1.0, 1.0);
			geo.wire(0, 70, -0.048, 0.021, -0.005, 0.035, 0.043, 0.014, 0.001, 1.0, 1.0);
			geo.wire(0, 50, -0.042, 0.008, 0.017, 0.017, -0.015, 0.014, 0.001, 1.0, 1.0);
			geo.wire(0, 66, 0.017, -0.015, 0.014, -0.027, 0.04, -0.031, 0.001, 1.0, 1.0);
			geo.wire(0, 85, -0.027, 0.04, -0.031, 0.046, -0.01, 0.028, 0.001, 1.0, 1.0);
			geo.wire(0, 47, 0.046, -0.01, 0.028, -0.013, -0.005, 0.031, 0.001, 1.0, 1.0);
			geo.wire(0, 70, 0.017, -0.015, 0.014, -0.048, -0.038, -0.04, 0.001, 1.0, 1.0);
			geo.wire(0, 77, -0.048, -0.038, -0.04, 0.049, -0.045, -0.04, 0.001, 1.0, 1.0);
			nec.geometry_complete(0,0);
			
			nec.gn_card(-1,0,0.0, 0.0, 0.0,0.0, 0.0, 0.0);
			nec.ld_card(5,0,0,0,3.72e7,0.0,0.0);
			nec.pt_card(-1, 0, 0, 0);
			nec.ex_card(1, 1, 1, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
			nec.fr_card(0, 2, 2400.0, 100.0);
			nec.rp_card(0, 1, 1, 500, 90.0, 90.0, 0.0, 0.0, 0.0, 0.0);
			double g = nec.get_maximum_gain();
			cout << "Gain @ 2.400 GHz is " << g << " dB" << endl;
		}
	}
	/* time the process */
	secnds( &stop_timer );
	stop_timer -= start_timer;
	cout << endl << endl;
	nec_float sec = stop_timer / 1000.0;
	
	cout << "Benchmark completed in " << sec << " seconds." << endl;
	cout << "Your computer's score is: " << 70.0 / sec << " NEC's" << endl;
	return 0;
}

/*-------------------------------------------------------------------*/


int nec_main( int argc, char **argv, nec_output_file& s_output )
{
	nec_output_flags s_output_flags;
	FILE *input_fp=NULL;
	FILE *output_fp=NULL;
	
	string input_filename, output_filename;
	
	char ain[3], line_buf[81];
	
	/* input card mnemonic list */
	/* "XT" stands for "exit", added for testing */
	#define CMD_NUM  20
	char *atst[CMD_NUM] =
	{
		"FR", "LD", "GN", "EX", "NT", "TL", \
		"XQ", "GD", "RP", "NX", "PT", "KH", \
		"NE", "NH", "PQ", "EK", "CP", "PL", \
		"EN", "WG"
	};
	
	// int ifrtmw, ifrtmp,  ix11;
	int mpcnt;
	int itmp3, itmp2, itmp4;
	
	int ain_num;    /* ain mnemonic as a number */
	
	nec_float tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
	nec_float ex_timer;
	
	/* getopt() variables */
	extern char *optarg;
	int option;
	
#ifndef __MINGW32__
	/*** signal handler related code ***/
	/* new and old actions for sigaction() */
	struct sigaction sa_new, sa_old;
	
	
	/* initialize new actions */
	sa_new.sa_handler = sig_handler;
	sigemptyset( &sa_new.sa_mask );
	sa_new.sa_flags = 0;
	
	/* register function to handle signals */
	sigaction( SIGINT,  &sa_new, &sa_old );
	sigaction( SIGSEGV, &sa_new, 0 );
	sigaction( SIGFPE,  &sa_new, 0 );
	sigaction( SIGTERM, &sa_new, 0 );
	sigaction( SIGABRT, &sa_new, 0 );
#endif
	
	/*** command line arguments handler ***/
	if ( argc == 1 )
	{
		usage();
		exit(-1);
	}
	
	bool results_to_stdout = false;
	
	/* process command line options */
	while( (option = getopt(argc, argv, "i:o:hvsgb") ) != -1 )
	{
		switch( option )
		{
		case 'i' : /* specify input file name */
			if ( strlen(optarg) > 75 )
			abort_on_error(-1);
			input_filename = optarg;
			break;
		
		case 'o' : /* specify output file name */
			if ( strlen(optarg) > 75 )
			abort_on_error(-2);
			output_filename = optarg;
			break;
		
		case 'g': /* return only the maximum gain to stdout */
			s_output_flags.set_gain_only(true);
			break;

		case 's': /* return output to stdout */
			results_to_stdout = true;
			break;
		
		case 'h' : /* print usage and exit */
			usage();
			exit(0);
		
		case 'v' : /* print nec2++ version */
			puts( "nec2++ " nec_version );
			exit(0);

		case 'b' : /* Run benchmark */
			benchmark();
			exit(0);
		
		default: /* print usage and exit */
			usage();
			exit(-1);
		
		} /* end of switch( option ) */	
	} /* while( (option = getopt(argc, argv, "i:o:hv") ) != -1 ) */

	/*** open input file ***/
	if ( (input_fp = fopen(input_filename.c_str(), "r")) == NULL )
	{
		string mesg = "nec2++: "  + input_filename;
		
		perror( mesg.c_str() );
		exit(-1);
	}

	/* make an output file name if not */
	/* specified by user on invocation */
	if ( output_filename  == "" )
	{
		/* strip the input file name extension if there is one */
		output_filename = input_filename.substr(0, input_filename.find(".",0)) + ".out";
	}
	
	/* open output file */
	if ( (output_fp = fopen(output_filename.c_str(), "w")) == NULL )
	{
		string mesg = "nec2++: "  + output_filename;
		
		perror( mesg.c_str() );
		exit(-1);
	}
	
	s_output.set_file(output_fp);
    
	secnds( &ex_timer );

	// allocate a new nec_context;
	nec_context s_context;
	s_context.set_output(s_output, s_output_flags);
	s_context.set_results_stdout(results_to_stdout);
	s_context.initialize();
    
  /* main execution loop, exits at various points */
  /* depending on error conditions or end of jobs */
	while( true )
	{
	
		s_output.end_section();
		s_output.set_indent(31);
		s_output.line(" __________________________________________");
		s_output.line("|                                          |");
		s_output.line("| NUMERICAL ELECTROMAGNETICS CODE (nec2++) |");
		s_output.line("| Translated to 'C++' in Double Precision  |");
		s_output.line("|              Version " nec_version "              |");
		s_output.line("|__________________________________________|");
	
		/* read a line from input file */
		if ( load_line(line_buf, input_fp) == EOF )
			abort_on_error(-3);
		
		/* separate card's id mnemonic */
		strncpy( ain, line_buf, 2 );
		ain[2] = '\0';
	
		/* If its an "XT" card, exit (used for debugging) */
		if ( strcmp(ain, "XT") == 0 )
		{
			nec_error_mode em(s_output);
			s_output.end_section();
			s_output.line("nec2++: Exiting after an \"XT\" command in main()");
			stop(0);
		}
	
		/* if its a "cm" or "ce" card start reading comments */
		if ( 	(strcmp(ain, "CM") == 0) ||
			(strcmp(ain, "CE") == 0) )
		{
			s_output.end_section();
			s_output.set_indent(31);
			s_output.line("---------------- COMMENTS ----------------");
			s_output.line(&line_buf[2]);
			while( strcmp(ain, "CM") == 0 )
			{
				/* read a line from input file */
				if ( load_line(line_buf, input_fp) == EOF )
					abort_on_error(-3);
			
				/* separate card's id mnemonic */
				strncpy( ain, line_buf, 2 );
				ain[2] = '\0';
			
				/* write comment to output file */
				s_output.line(&line_buf[2]);
			} /* while( strcmp(ain, "CM") == 0 ) */
		
			/* no "ce" card at end of comments */
			if ( strcmp(ain, "CE") != 0 )
			{
				s_output.end_section();
				s_output.line("  ERROR: INCORRECT LABEL FOR A COMMENT CARD");
				abort_on_error(-4);
			}
		} /* if ( strcmp(ain, "CM") == 0 ... */
		else
		{
			rewind( input_fp );
		}
	
		/* initializations etc from original fortran code */
		mpcnt=0;
		
		/* set up geometry data in subroutine parse_geometry */
		c_geometry& geo = s_context.get_geometry();
		geo.parse_geometry(&s_context, input_fp);
		
		s_context.calc_prepare();
		s_output.end_section();
		
/*		xtemp.resize(n_plus_m);
		ytemp.resize(n_plus_m);
		ztemp.resize(n_plus_m);
		sitemp.resize(n_plus_m);
		bitemp.resize(n_plus_m);
*/		
		
		/*
			Main input section, exits at various points
			depending on error conditions or end of job.
			This is called the card input loop.
		*/
				
		bool next_job = false; /* start next job (next structure) flag */
		while( ! next_job )
		{
			int itmp1;
			
			//jmp_iloop = false;
			/* main input section - standard read statement - jumps */
			/* to appropriate section for specific parameter set up */
			readmn(input_fp, output_fp, ain, &itmp1, &itmp2, &itmp3, &itmp4,
				&tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6 );
			
			/* If its an "XT" card, exit */
			if ( strcmp(ain, "XT" ) == 0 )
			{
				nec_error_mode em(s_output);
				s_output.endl();
				s_output.line("nec2++: Exiting after an \"XT\" command in main()" );
				stop(0);
			}
			
			mpcnt++;
			fprintf( output_fp,
				"\n*****  DATA CARD N0. %3d "
				"%s %3d %5d %5d %5d %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E",
				mpcnt, ain, itmp1, itmp2, itmp3, itmp4,
				tmp1, tmp2, tmp3, tmp4, tmp5, tmp6 );
			
			/* identify card id mnemonic (except "ce" and "cm") */
			for( ain_num = 0; ain_num < CMD_NUM; ain_num++ )
				if ( strncmp( ain, atst[ain_num], 2) == 0 )
					break;
		
			/* take action according to card id mnemonic */
			switch( ain_num )
			{
			case 0: /* "fr" card, frequency parameters */
				s_context.fr_card(itmp1, itmp2, tmp1, tmp2);
				continue;
		
			case 1: /* "ld" card, loading parameters */
				s_context.ld_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3);
				continue;
		
			case 2: /* "gn" card, ground parameters under the antenna */
				s_context.gn_card(itmp1, itmp2, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				continue; 
		
			case 3: /* "ex" card, excitation parameters */
				s_context.ex_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				continue; /* continue card input loop */
		
			case 4: /* "nt" card, network parameters */
				s_context.nt_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				continue; /* continue card input loop */
		
			case 5: /* "tl" card, network parameters */
				s_context.tl_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				continue; /* continue card input loop */
		
			case 6: /* "xq" execute card - calc. including radiated fields */
				s_context.xq_card(itmp1);
				continue;
		
			case 7: /* "gd" card, ground representation */
				s_context.gd_card(tmp1, tmp2, tmp3, tmp4);
				continue; /* continue card input loop */
		
			case 8: /* "rp" card, standard observation angle parameters */
				s_context.rp_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				continue; /* was break; followed by special code */
		
			case 9: /* "nx" card, do next job */
				next_job = true;
				continue; /* continue card input loop */
		
			case 10: /* "pt" card, print control for current */
				s_context.pt_card(itmp1, itmp2, itmp3, itmp4);
				continue; /* continue card input loop */
		
			case 11: /* "kh" card, matrix integration limit */
				s_context.kh_card(tmp1);
				continue; /* continue card input loop */
		
			case 12:  /* "ne" card, near field calculation parameters */
				s_context.ne_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				continue;
			
			case 13:  /* "nh" card, near field calculation parameters */
				s_context.nh_card(itmp1, itmp2, itmp3, itmp4, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);
				continue;
		
			case 14: /* "pq" card, print control for charge */
				s_context.pq_card(itmp1, itmp2, itmp3, itmp4);
				continue; /* continue card input loop */
		
			case 15: /* "ek" card,  extended thin wire kernel option */
				s_context.ek_card(itmp1);
				continue; /* continue card input loop */
		
			case 16: /* "cp" card, maximum coupling between antennas */
				s_context.cp_card(itmp1, itmp2, itmp3, itmp4);
				continue; /* continue card input loop */
		
			case 17: /* "pl" card, plot flags */
				{
					std::string ploutput_filename(input_filename);
					ploutput_filename += ".plt";
					try
					{
						s_context.pl_card(ploutput_filename.c_str(), itmp1, itmp2, itmp3, itmp4);
					}
					catch(int x)
					{
						char mesg[88] = "nec2c: ";
					
						strcat( mesg, ploutput_filename.c_str() );
						perror( mesg );
						exit(-1);
					}
				}
				continue; /* continue card input loop */
		
			case 19: /* "wg" card, not supported */
				abort_on_error(-5);
		
			default:
				if ( ain_num != 18 ) // EN card
				{
					fprintf( output_fp,
						"\n\n  FAULTY DATA CARD LABEL AFTER GEOMETRY SECTION" );
					stop(-1);
				}
			
				/******************************************************
				*** normal exit of nec2++ when all jobs complete ok ***
				******************************************************/
				s_context.all_jobs_completed();
				
				/* time the process */
				secnds( &tmp1 );
				tmp1 -= ex_timer;
				fprintf( output_fp, "\n\n  TOTAL RUN TIME: %d msec", (int)tmp1 );
				stop(0);
			} /* switch( ain_num ) */
		
			/*
				End of the main input section. 
			
				far_field_flag is true if last card was XQ or RP
				
				This is no longer used, but I am leaving it in here while I iron this out
				properly. simulate() should be called by the xq card and the rp card.
			*/
			bool far_field_flag = (ain_num == 6) || (ain_num == 8);
			ASSERT(false == far_field_flag);
			s_context.simulate(false);
		} /* while( ! next_job ) */
	
	} /* while(true)  */

	return(0);
} /* end of nec_main() */


/*-----------------------------------------------------------------------*/

void readmn(FILE* input_fp, FILE* output_fp, char *gm, int *i1, int *i2, int *i3, int *i4,
    nec_float *f1, nec_float *f2, nec_float *f3,
    nec_float *f4, nec_float *f5, nec_float *f6 )
{
	char line_buf[134];
	int line_idx;
	int n_integer_params = 4, n_float_params = 6;
	int integer_array[4] = { 0, 0, 0, 0 };
	nec_float real_array[6] = { 0., 0., 0., 0., 0., 0. };
	
	/* read a line from input file */
	int eof = load_line( line_buf, input_fp );
	
	/* get line length */
	int line_length = strlen( line_buf );
	
	/* abort if card's mnemonic too short or missing */
	if ( line_length < 2 )
	{
		if (EOF == eof)
		{
			// insert and EN card if we get to an end of file
			strncpy( gm, "EN", 2 );
			return;
		}
		else
		{
			fprintf( output_fp,
				"\n  COMMAND DATA CARD ERROR:"
				"\n  CARD'S MNEMONIC CODE TOO SHORT OR MISSING." );
			stop(-1);
		}
	}
	
	/* extract card's mnemonic code */
	strncpy( gm, line_buf, 2 );
	gm[2] = '\0';
	
	/* Exit if "XT" command read (for testing) */
	if ( strcmp( gm, "XT" ) == 0 )
	{
		fprintf( stderr,
			"\nnec2++: Exiting after an \"XT\" command in read_geometry_card()\n" );
		fprintf( output_fp,
			"\n\n  nec2++: Exiting after an \"XT\" command in read_geometry_card()" );
		stop(0);
	}
	
	/* Return if only mnemonic on card */
	if ( line_length == 2 )
	{
		*i1 = *i2 = *i3 = *i4 = 0;
		*f1 = *f2 = *f3 = *f4 = *f5 = *f6 = 0.0;
		return;
	}
	
	/* read integers from line */
	line_idx = 1;
	for (int i = 0; i < n_integer_params; i++ )
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
			*i3= integer_array[2];
			*i4= integer_array[3];
			*f1= real_array[0];
			*f2= real_array[1];
			*f3= real_array[2];
			*f4= real_array[3];
			*f5= real_array[4];
			*f6= real_array[5];
			return;
		}
		
		/* read an integer from line */
		integer_array[i] = atoi( &line_buf[line_idx] );
		
		/* traverse numerical field to next ' ' or ',' or '\0' */
		line_idx--;
		while( (line_buf[++line_idx] != ' ') &&
			(line_buf[  line_idx] != ',') &&
			(line_buf[  line_idx] != '\0') )
		{
			/* test for non-numerical characters */
			if ( ((line_buf[line_idx] <  '0')  ||
				(line_buf[line_idx] >  '9')) &&
				(line_buf[line_idx] != '+')  &&
				(line_buf[line_idx] != '-') )
			{
				fprintf( output_fp,
				"\n  COMMAND DATA CARD \"%s\" ERROR:"
				"\n  NON-NUMERICAL CHARACTER '%c' IN INTEGER FIELD AT CHAR. %d\n",
				gm, line_buf[line_idx], (line_idx+1) );
				stop(-1);
			}		
		} /* while( (line_buff[++line_idx] ... */
		
		/* Return on end of line */
		if ( line_buf[line_idx] == '\0' )
		{
			*i1= integer_array[0];
			*i2= integer_array[1];
			*i3= integer_array[2];
			*i4= integer_array[3];
			*f1= real_array[0];
			*f2= real_array[1];
			*f3= real_array[2];
			*f4= real_array[3];
			*f5= real_array[4];
			*f6= real_array[5];
			return;
		}		
	} /* for( i = 0; i < n_integer_params; i++ ) */
	
	/* read nec_floats from line */
	for (int i = 0; i < n_float_params; i++ )
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
			*i3= integer_array[2];
			*i4= integer_array[3];
			*f1= real_array[0];
			*f2= real_array[1];
			*f3= real_array[2];
			*f4= real_array[3];
			*f5= real_array[4];
			*f6= real_array[5];
			return;
		}
		
		/* read a nec_float from line */
		real_array[i] = atof( &line_buf[line_idx] );
		
		/* traverse numerical field to next ' ' or ',' */
		line_idx--;
		while( (line_buf[++line_idx] != ' ') &&
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
				fprintf( output_fp,
				"\n  COMMAND DATA CARD \"%s\" ERROR:"
				"\n  NON-NUMERICAL CHARACTER '%c' IN FLOAT FIELD AT CHAR. %d\n",
				gm, line_buf[line_idx], (line_idx+1) );
				stop(-1);
			}		
		} /* while( (line_buff[++line_idx] ... */
		
		/* Return on end of line */
		if ( line_buf[line_idx] == '\0' )
		{
			*i1= integer_array[0];
			*i2= integer_array[1];
			*i3= integer_array[2];
			*i4= integer_array[3];
			*f1= real_array[0];
			*f2= real_array[1];
			*f3= real_array[2];
			*f4= real_array[3];
			*f5= real_array[4];
			*f6= real_array[5];
			return;
		}		
	} /* for( i = 0; i < n_float_params; i++ ) */
	
	*i1= integer_array[0];
	*i2= integer_array[1];
	*i3= integer_array[2];
	*i4= integer_array[3];
	*f1= real_array[0];
	*f2= real_array[1];
	*f3= real_array[2];
	*f4= real_array[3];
	*f5= real_array[4];
	*f6= real_array[5];
}


/*-----------------------------------------------------------------------*/


#ifndef __MINGW32__
static void sig_handler(int signal )
{
	switch( signal )
	{
		case SIGINT :
			fprintf(stderr, "nec2++: exiting via user interrupt" );
			exit( signal );
	
		case SIGSEGV :
			fprintf(stderr, "nec2++: segmentation fault" );
			exit( signal );
	
		case SIGFPE :
			fprintf(stderr, "nec2++: floating point exception" );
			exit( signal );
	
		case SIGABRT :
			fprintf(stderr, "nec2++: abort signal received" );
			exit( signal );
	
		case SIGTERM :
			fprintf(stderr, "nec2++: termination request received" );
			stop( signal );
	}
}
#endif

