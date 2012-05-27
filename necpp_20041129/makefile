###########################################################################
#	
#	MAKEFILE for the C++ port of NEC
#	
#	Author: Tim Molteno.
#
#
#
#	Copyright (C) 2004  Timothy C.A. Molteno
#	
#	This program is free software; you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation; either version 2 of the License, or
#	(at your option) any later version.
#	
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	
#	You should have received a copy of the GNU General Public License
#	along with this program; if not, write to the Free Software
#	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#

###################################################################################
#
#	ATLAS Instructions
#		
#	NEC2++ depends on ATLAS. It is possible to build it without
#	atlas at the moment -- but this will disappear.
#
#	apt-get install atlas3-base atlas3-headers atlas3-base-dev
#	apt-get install refblas3-dev lapack3-dev lapack3-doc
#
#
#	UNIX install instructions
#
#	make clean nec2++
#	sudo make install
#

HDRS 	= src/nec2cpp.h src/misc.h src/c_ggrid.h src/c_evlcom.h \
	src/c_plot_card.h src/math_util.h src/common.h src/matrix_algebra.h \
	src/nec_output.h src/c_geometry.h src/nec_results.h

SRC 	= src/nec_context.cpp src/nec2cpp.cpp \
	src/misc.cpp src/c_ggrid.cpp src/c_evlcom.cpp src/c_plot_card.cpp \
	src/matrix_algebra.cpp src/nec_output.cpp src/c_geometry.cpp \
	src/nec_ground.cpp src/nec_radiation_pattern.cpp 
# net_solve.cpp

OBJS 	= $(SRC:.cpp=.o)

all:	nec2cpp

###############################################################################
#
#	The following is the output of the 'make depends' target.
#
#
nec2cpp.o: src/nec2cpp.cpp src/nec2cpp.h src/math_util.h src/safe_array.h \
	src/misc.h src/nec_context.h src/common.h src/c_ggrid.h src/c_evlcom.h \
	src/c_geometry.h src/nec_output.h src/matrix_algebra.h src/electromag.h \
	src/nec_radiation_pattern.h src/nec_results.h src/nec_ground.h \
	src/c_plot_card.h

TARGET	= nec2cpp
NEC2dx	= nec2

LDATLAS =-L /usr/lib/atlas -lblas  -latlas  -llapack_atlas

#	Uncomment Following to use Gcc 3.4
#CC=gcc-3.4
#CXX=g++-3.4

CXXFLAGS= -Wall -g $(CATLAS) -DNEC_ERROR_CHECK

CATLAS =
LDFLAGS=
#CATLAS =-DATLAS
#LDFLAGS=$(LDATLAS)

nec2cpp:	$(OBJS) $(HDRS)
	$(CC) -o $(TARGET) $(OBJS) -lstdc++ $(LDFLAGS) 

necDiff:	src/necDiff.cpp src/AntennaInput.h src/RadiationInput.h src/PowerBudget.h src/BaseInput.h
	$(CC) $(CXXFLAGS)  -o necDiff src/necDiff.cpp -lstdc++

nec2:
	g77 -o $(NEC2dx) nec2dx.f

.PHONY : clean
clean:
	rm -f $(OBJS) nec2_prof nec2++ nec2++.exe $(TARGET) src/*.o

depends:
	$(CC) -MM src/nec2cpp.cpp


############################################################################
#
#	Target for builing with the release version
#
# 	NOTE: Optimization flags will produce incorrect code on -O3.
#			Optimization -02 now appears to work fine.
#
#	To use ALTAS, uncomment the second line

nec2++: $(SRC) $(HDRS)
	make nec2cpp 'TARGET=nec2++' 'CXXFLAGS= -Wall -O2 -march=athlon-xp'
#	make nec2cpp 'TARGET=nec2++' 'CXXFLAGS= -DATLAS -Wall -O2 -march=athlon-xp' 'LDFLAGS=$(LDATLAS)'

test_fast:	
	make test 'TARGET=nec2++'
	
install:
	cp ./nec2++ /usr/local/bin

	
###########################################################################################################
#
#	This builds a source code distribution (and a Windows executable) and packages it into
#	an archive. Also uploads it to the server
#

BUILDDIR=necpp_`date +%Y%m%d`
DEST=elec@www.physics.otago.ac.nz:electronics/nec
ARCHIVE=${BUILDDIR}.tar.bz2

dist:
	rm -rf ${BUILDDIR}
	mkdir -p ${BUILDDIR}/src
	cp src/*.cpp ${BUILDDIR}/src
	cp src/*.h ${BUILDDIR}/src
	cp makefile ${BUILDDIR}
	cp README ${BUILDDIR}
	mkdir -p ${BUILDDIR}/test_data
	cp test_data/*.nec ${BUILDDIR}/test_data
	tar -cf ${BUILDDIR}.tar ${BUILDDIR}
	rm -f ${ARCHIVE}
	bzip2 -9 ${BUILDDIR}.tar
#
	rm -rf www
	mkdir -p www
	cp ${ARCHIVE} www/necpp_latest.tar.bz2
	mv ${ARCHIVE} www
#	
	mv makefile makefile.old
	make -f Makefile.cvs
	./configure
	make -f Makefile dist
	mv makefile.old makefile
	mv necpp-1.0.4.tar.gz www
#
	make -f makefile clean
	make -f makefile nec2++.exe
	zip nec2++.zip nec2++.exe
	mv nec2++.zip www
	cp docs/*.html www
#	cp docs/users_guide.pdf www
	rsync --verbose --rsh=ssh  --cvs-exclude --recursive www/* ${DEST}


	
############################################################################
#
#	Target for builing a Win32 version with MinGW32
#
MINGW_FLAGS = -D__MINGW32__ -D__WIN32__=1  -mwindows

nec2++.exe:  $(SRC) $(HDRS)
	make nec2cpp 'TARGET=nec2++.exe' 'CC=i586-mingw32msvc-gcc' 'CXX=i586-mingw32msvc-c++' 'CXXFLAGS= -Wall -O2 $(MINGW_FLAGS)'
	strip --strip-all nec2++.exe

############################################################################
#
#	Target for builing profiling
#
#
# NOTE: Optimization flags will produce incorrect code on  O3.
#	Optimization -02 now appears to work fine. -DATLAS  removed
nec2_prof:
	echo "Making Profiling Target"
	make nec2cpp 'TARGET=nec2_prof' 'CXXFLAGS= -pg -Wall -O2 -mcpu=athlon-xp' 'LDFLAGS= -pg $(LDATLAS)'

test_prof:	
#	make test 'TARGET=nec2_prof'
#	./nec2_prof -i test_data/buoy.nec -o test_data/buoy.outcpp
	./nec2_prof -i test_data/test299.nec -o test_data/test299.outcpp
	gprof nec2_prof gmon.out > gprof.txt
	kwrite gprof.txt	

############################################################################
#	
#	TESTHARNESS
#
#	This tests the output of the C++ port against the nec2c, C 
#	port as well as the original FORTRAN
#

NEC_TESTS = test_data/buoy.nec \
	test_data/GA487.nec \
	test_data/ga_pjw_0.nec \
	test_data/ga_pjw_1.nec \
	test_data/Yg_4x6_b.nec \
	test_data/36dip.nec \
	test_data/K5332187.nec \
	test_data/Gs_8d_bb.nec \
	test_data/example2.nec \
	test_data/example3.nec \
	test_data/example4.nec \
	test_data/example5.nec \
	test_data/example6.nec \
	test_data/sommerfield1.nec \
	test_data/yagi16.nec \
	test_data/yagi.nec \
	test_data/test299.nec

SHORT_TESTS = test_data/example3.nec

############################################################################
#	
#	Rewrite Rules are used to automate the generation
#	and comparison of output files.

#	Generate an output file from the FORTRAN nec2 (not working)
%.nec %.out:
	./nec2 -i $*.nec -o $(*F)
	
# 	Generate an output file from nec2c
%.nec %.outc:
	- ./nec2c -i $*.nec -o temp.$(*F)
	mv temp.$(*F).out $@
	
# 	Generate an output file from nec2cpp
%.nec %.outcpp:
	- ./$(TARGET) -i $*.nec -o $@

#	Compare two outputs (c and cpp)
%.nec %.diff: %.outc %.outcpp
	./necDiff $*.outc $*.outcpp > $@

#	Compare two outputs (FORTRAN and cpp)
%.nec %.dif: %.out %.outcpp
	./necDiff $*.out $*.outcpp > $@

test:
	rm -f *.diff
	rm -f $(NEC_TESTS:.nec=.outcpp)
	rm -f $(NEC_TESTS:.nec=.outc)
	make test_aux

test_short:
	./$(TARGET) -j 1 -i test_data/example3.nec -o test_data/example3.outcpp


# 	Change the following line to
#	test_aux:	necDiff nec2cpp $(NEC_TESTS:.nec=.outcpp) $(NEC_TESTS:.nec=.outc) $(NEC_TESTS:.nec=.diff) $(NEC_TESTS:.nec=.dif)
#	in order to include the FORTRAN code in the test suite.

test_aux:	necDiff $(TARGET) $(NEC_TESTS:.nec=.outcpp) $(NEC_TESTS:.nec=.outc) $(NEC_TESTS:.nec=.diff)
	echo "Test Complete"
	

# A target for testing nec2cpp against the FORTRAN original.
	
test_f:	necDiff nec2cpp nec2
	./$(TARGET) -i yagi.nec -o yagi++
	./$(NEC2dx) -i yagi.nec -o yagi
	./$(TARGET) -i yagi16.nec -o yagi16++
	./$(NEC2dx) -i yagi16.nec -o yagi16
	./$(TARGET) -i TEST299.NEC -o TEST299++
	./$(NEC2dx) -i TEST299.NEC -o TEST299
	./$(TARGET) -i example2.nec -o example2++
	./$(NEC2dx) -i example2.nec -o example2
	./$(TARGET) -i example3.nec -o example3++
	./$(NEC2dx) -i example3.nec -o example3
	./necDiff example2++.out example2.out
	./necDiff example3++.out example3.out
	./necDiff yagi++.out yagi.out
	./necDiff yagi16++.out yagi16.out
	./necDiff TEST299++.out TEST299.out
# DO NOT DELETE
