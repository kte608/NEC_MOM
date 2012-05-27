#!/usr/bin/env python
from KStandard import attributesFromDict
from numpy import *
from necpp import *


def Example1():
    # Page 83 of the manual
    ls=LineSeg(Vector((0,0,-0.25)),Vector((0,0,0.25)))
    #Visualization([ls])
    mod=NECModel(description="EXAMPLE 1. CENTER FED LINEARA ANTENNA \nPAGE 83 OF THE MANUAL",
                 geometry={"wiresegs":[WireSeg(ls.startvect.xyz(),ls.endvect.xyz(),
                                               radius=0.001,
                                               numSplit=7,
                                               load=SeriesRLCLoad(R=10,L=3e-9,C=5.3e-11,LDTAGF=4))]},
                 excitation=VoltageSource(0,segsplit=4,ImpedanceNormalization=None),
                 freq=299.8e6)
    mod.runSim(printCardDeck=True,printRawResult=True)

def Example2(plot="a"):
    from KPlot import KBodePlot
    # Page 87 of the manual
    ls=LineSeg(Vector((0,0,-0.25)),Vector((0,0,0.25)))
    #Visualization([ls])
    descrip=        "EXAMPLE 2. CENTER FED LINEAR ANTENNA.\n"
    descrip=descrip+"CURRENT CLOPE DISCONTINUITY SOURCE\n"
    descrip=descrip+"THIN PERFECTLY CONDUCTING WIRE\nPAGE 87"
    mod=NECModel(description=descrip,
                 geometry={"wiresegs":[WireSeg(ls.startvect.xyz(),ls.endvect.xyz(),
                                               radius=0.00001,
                                               numSplit=8,
                                               load=None)]},
                 excitation=VoltageSource(0,segsplit=5,
                                          calcMaxAdmittanceMatrixAssym=1,
                                          extype="CURRENTSLOPE"),
                 freq=LinearFrequencyRange(freq=200e6,dfreq=10e6,num=20))
    mod.runSim()
    #print mod.RESULT.freqResults[0][1]["ANTENNA INPUT PARAMETERS"]['TAG{0}:SEG{5}']
    tmp=[[v[0],v[1]["ANTENNA INPUT PARAMETERS"]['TAG{0}:SEG{5}']["IMPEDANCE"]] for v in mod.RESULT.freqResults]
    freqs,impedances=transpose(tmp)
    if plot=="a":
        KBodePlot([(freqs,impedances,"input impedance")],scaletype="value_linear")
    #print mod.ProduceCardDeck()
    descrip=        "EXAMPLE 2. CENTER FED LINEAR ANTENNA.\n"
    descrip=descrip+"CURRENT CLOPE DISCONTINUITY SOURCE\n"
    descrip=descrip+"THIN ALUMINIUM WIRE\nPAGE 87"
    mod.set_description(descrip)
    mod.geometry["wiresegs"][0].set_load(WireConductivityLoad(3.720e7))
    mod.set_freq(LinearFrequencyRange(freq=300e6,dfreq=50e6,num=1))
    mod.runSim()
    currents=[v[2] for v in  mod.RESULT.freqResults[0][1]["CURRENTS AND LOCATION"]]
    positions=[i for i in range(len(currents))]
    #print currents
    if plot=="b":
        KBodePlot([(positions,currents,"currents")],linearfreqs=True,scaletype="value_linear")
    #print mod.ProduceCardDeck()

def Example3():
    # page 92
    pass

def Example4():
    # page 98
    pass

def Example5():
    # page 102
    descrip=        "EXAMPLE 5: 12 ELEMENT LOG PERIODIC ANTENNA IN FREE SPACE.\n"
    descrip=descrip+"78 SEGMENTS. SIGMA=O/L   RECEIVING AND TRANS. PATTERNS.\n"
    descrip=descrip+"DIPOLE LENGTH TO DIAMETER RATIO=150.\n"
    descrip=descrip+"TAU=0.93.  SIGMA=0.70.  BOOM IMPEDANCE=50. OHMS. \nPAGE 102"

    linesegs=[Cylinder(Vector((0      ,-1     ,0))  ,Vector((0     ,1     ,0))    ,0.00667),
              Cylinder(Vector((-.7527 ,-1.0753,0))  ,Vector((-.7527,1.0753,0))    ,0.00717),
              Cylinder(Vector((-1.562 ,-1.1562,0))  ,Vector((-1.562,1.1562,0))    ,0.00771),
              Cylinder(Vector((-2.4323,-1.2432,0))  ,Vector((-2.4323,1.2432,0))   ,0.00829),
              Cylinder(Vector((-3.368,-1.3368,0))   ,Vector((-3.368,1.3368,0))    ,0.00891),
              Cylinder(Vector((-4.3742,-1.4374,0))  ,Vector((-4.3742,1.4374,0))   ,0.00958),
              Cylinder(Vector((-5.4562,-1.5456,0))  ,Vector((-5.4562,1.5456,0))   ,0.0103),
              Cylinder(Vector((-6.6195,-1.6619,0))  ,Vector((-6.6195,1.6619,0))   ,0.01108),
              Cylinder(Vector((-7.8705,-1.787,0))   ,Vector((-7.8705,1.787,0))    ,0.01191),
              Cylinder(Vector((-9.2156,-1.9215,0))  ,Vector((-9.2156,1.9215,0))   ,0.01281),
              Cylinder(Vector((-10.6619,-2.0662,0)) ,Vector((-10.6619,2.0662,0))  ,0.01377),
              Cylinder(Vector((-12.2171,-2.2217,0)) ,Vector((-12.2171,2.2217,0))  ,0.01481)
              ]
    Zline=-50
    tlines=[TransmissionLineSeg(0,1,Port1SEG=3,Port2SEG=3,Z0=Zline),
            TransmissionLineSeg(1,2,Port1SEG=3,Port2SEG=3,Z0=Zline),
            TransmissionLineSeg(2,3,Port1SEG=3,Port2SEG=3,Z0=Zline),
            TransmissionLineSeg(3,4,Port1SEG=3,Port2SEG=3,Z0=Zline),
            TransmissionLineSeg(4,5,Port1SEG=3,Port2SEG=4,Z0=Zline),
            TransmissionLineSeg(5,6,Port1SEG=4,Port2SEG=4,Z0=Zline),
            TransmissionLineSeg(6,7,Port1SEG=4,Port2SEG=4,Z0=Zline),
            TransmissionLineSeg(7,8,Port1SEG=4,Port2SEG=4,Z0=Zline),
            TransmissionLineSeg(8,9,Port1SEG=4,Port2SEG=4,Z0=Zline),
            TransmissionLineSeg(9,10,Port1SEG=4,Port2SEG=5,Z0=Zline),
            TransmissionLineSeg(10,11,Port1SEG=5,Port2SEG=5,Z0=Zline,Port2ShuntAdmittance=complex(0,0.02))
            ]
            
    #Visualization(linesegs)
    splits=[5,5,5,5,5,7,7,7,7,7,9,9]
    wiresegs=[WireSeg(l.startvect.xyz(),l.endvect.xyz(),radius=l.radius,numSplit=s) for l,s in zip(linesegs,splits)]
    
    mod=NECModel(description=descrip,
                 geometry={"wiresegs":wiresegs,"transmissionlines":tlines},
                 excitation=VoltageSource(0,segsplit=3,
                                          calcMaxAdmittanceMatrixAssym=1,
                                          calcImpedanceTable=0),
                 freq=LinearFrequencyRange(freq=46.29e6,dfreq=0,num=0),
                 radiationPattern=RadiationPatternSpherical(numTheta=37,numPhi=1,thetaS=pi/2,
                                                            thetaE=-pi/2,phiS=0,phiE=0,
                                                            gainNormalization="MAJOR",
                                                            gainType="DIRECTIVE",
                                                            gainAveraging=True))
    mod.runSim(printCardDeck=True,printRawResult=True)
def Example6():
    # page 109
    pass

def Example7and8():
    # page 117
    pass



def _ParallelRLC(plotvis=False,frequency=298e6,verbose=True):
    R=50    # Ohms
    L=9e-9  # Henries
    C=2e-12 # Farads
    wireradius=0.00002
    seglength=0.0002
    # The voltage source (Seg 0)
    segs=[Cylinder(Vector((0.0       ,-seglength*0.5,0)),
                   Vector((0.0       , seglength*0.5,0)),
                   wireradius,color=(1,0,1))]
    # connecting wires
    segs.append(Cylinder(Vector((0.0       ,seglength*0.5,0)),
                         Vector((seglength ,seglength*0.5,0)),
                         wireradius,color=(1,1,1)))
    segs.append(Cylinder(Vector((seglength ,-seglength*0.5,0)),
                         Vector((0.0       ,-seglength*0.5,0)),
                         wireradius,color=(1,1,1)))
    
    # the resistor (Seg 3)
    segs.append(Cylinder(Vector((seglength ,seglength*0.5,0)),
                         Vector((seglength ,-seglength*0.5,0)),
                         wireradius,color=(1,1,1)))
    
    # connecting wires
    segs.append(Cylinder(Vector((seglength ,seglength*0.5,0)),
                         Vector((2*seglength ,seglength*0.5,0)),
                         wireradius,color=(1,1,1)))
    segs.append(Cylinder(Vector((2*seglength ,-seglength*0.5,0)),
                         Vector((seglength       ,-seglength*0.5,0)),
                         wireradius,color=(1,1,1)))
    

    # the Inductor (Seg 6)
    segs.append(Cylinder(Vector((2*seglength ,seglength*0.5,0)),
                         Vector((2*seglength ,-seglength*0.5,0)),
                         wireradius,color=(1,1,1)))
    
    
    # connecting wires
    segs.append(Cylinder(Vector((2*seglength ,seglength*0.5,0)),
                         Vector((3*seglength ,seglength*0.5,0)),
                         wireradius,color=(1,1,1)))
    segs.append(Cylinder(Vector((3*seglength ,-seglength*0.5,0)),
                         Vector((2*seglength       ,-seglength*0.5,0)),
                         wireradius,color=(1,1,1)))
    
    # the Capacitor (Seg 9)
    segs.append(Cylinder(Vector((3*seglength ,seglength*0.5,0)),
                         Vector((3*seglength ,-seglength*0.5,0)),
                         wireradius,color=(1,1,1)))
    
    wiresegs=[]
    for i,l in enumerate(segs):
        if i==3: # The resistor
            load=SeriesRLCLoad(R=50)#SimpleImpedanceLoad(complex(50.0,0.0))#WireConductivityLoad(2)#
            l=Cylinder(l.startvect,l.endvect,l.radius,color=(1,0,0))
            segs[i]=l
        elif i==6: # The inductor
            load=SeriesRLCLoad(L=26.7e-9)#SimpleImpedanceLoad(complex(0.0,50.0))#SeriesRLCLoad(R=50)
            l=Cylinder(l.startvect,l.endvect,l.radius,color=(0,1,0))
            segs[i]=l
        elif i==9: # The capacitor
            load=SeriesRLCLoad(C=10.68e-12)#SimpleImpedanceLoad(complex(0.0,-50.0))#SeriesRLCLoad(R=50)
            l=Cylinder(l.startvect,l.endvect,l.radius,color=(0,0,1))
            segs[i]=l
        else:
            load=None
        print i,load
        wiresegs.append(WireSeg(l.startvect.xyz(),l.endvect.xyz(),radius=l.radius,
                                load=load))
        
    
    #mod.geometry["wiresegs"][0].set_load(WireConductivityLoad(3.720e7))
    if plotvis: Visualization(segs)

        
    mod=NECModel(geometry={"wiresegs":wiresegs},
                 excitation=VoltageSource(0),
                 freq=LinearFrequencyRange(freq=frequency))
                 
    mod.runSim(printCardDeck=verbose,printRawResult=verbose)
    #print mod.RESULT.freqResults
    currents=[v[2] for v in  mod.RESULT.freqResults[0][1]["CURRENTS AND LOCATION"]]
    print "\n################"
    print currents
    positions=[i for i in range(len(currents))]
    
    
    from KPlot import KBodePlot
    KBodePlot([(positions,currents,"currents")],linearfreqs=True,scaletype="value_linear")
    #return currents[:loop2NumStart],currents[loop2NumStart+1:]


if __name__ == '__main__': 
    Example1()
    _ParallelRLC()
    #Example2(plot="b")
    #Example5()

