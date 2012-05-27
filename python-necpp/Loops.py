#!/usr/bin/env python
from KStandard import attributesFromDict
from numpy import *
from necpp import *

def BuildArc(phiStart,phiEnd,num,loopradius=0.05,wireradius=0.001):
    result=[]
    dphi=(phiEnd-phiStart)/float(num)
    for i in range(num):
        result.append(Cylinder(Vector((loopradius,phiStart+dphi*i,0),vtype="CYLINDRICAL"),
                               Vector((loopradius,phiStart+dphi*(i+1),0),vtype="CYLINDRICAL"),
                               wireradius))
    return result

def KarlLoop():
    # I want to make a simple loop analysis 
    # with the capacitors and whatnot.
    wireradius=0.001
    loopRadius=0.05
    
    def AdmittanceSeriesZ(Z):
        """
        Despite how it looks:
        Ymn does not always equal Zmn
        """
        result=[[1/Z,-1/Z],
                [-1/Z,1/Z]]
        return result

    
    capSize=0.005
    capArc=atan(capSize/float(loopRadius))
    lineSegs=[Cylinder(Vector((loopRadius        ,2*pi-capArc*0.5,0),vtype="CYLINDRICAL"),
                       Vector((loopRadius,0,0),vtype="CYLINDRICAL"),
                       wireradius),
              Cylinder(Vector((loopRadius        ,0,0),vtype="CYLINDRICAL"),
                       Vector((loopRadius,capArc*0.5,0),vtype="CYLINDRICAL"),
                       wireradius),
              
              Cylinder(Vector((loopRadius        ,pi*0.5-capArc*0.5,0),vtype="CYLINDRICAL"),
                       Vector((loopRadius-capSize,pi*0.5-capArc*0.5,0),vtype="CYLINDRICAL"),
                       wireradius),
              Cylinder(Vector((loopRadius        ,pi*0.5+capArc*0.5,0),vtype="CYLINDRICAL"),
                       Vector((loopRadius-capSize,pi*0.5+capArc*0.5,0),vtype="CYLINDRICAL"),
                       wireradius),
              
              Cylinder(Vector((loopRadius        ,pi-capArc*0.5,0),vtype="CYLINDRICAL"),
                       Vector((loopRadius-capSize,pi-capArc*0.5,0),vtype="CYLINDRICAL"),
                       wireradius),
              Cylinder(Vector((loopRadius        ,pi+capArc*0.5,0),vtype="CYLINDRICAL"),
                       Vector((loopRadius-capSize,pi+capArc*0.5,0),vtype="CYLINDRICAL"),
                       wireradius),

              Cylinder(Vector((loopRadius        ,1.5*pi-capArc*0.5,0),vtype="CYLINDRICAL"),
                       Vector((loopRadius-capSize,1.5*pi-capArc*0.5,0),vtype="CYLINDRICAL"),
                       wireradius),
              Cylinder(Vector((loopRadius        ,1.5*pi+capArc*0.5,0),vtype="CYLINDRICAL"),
                       Vector((loopRadius-capSize,1.5*pi+capArc*0.5,0),vtype="CYLINDRICAL"),
                       wireradius)
              ]
    lineSegs.extend(BuildArc(capArc*0.5,pi*0.5-capArc*0.5,10,wireradius=wireradius,loopradius=loopRadius))
    lineSegs.extend(BuildArc(pi*0.5+capArc*0.5,pi-capArc*0.5,10,wireradius=wireradius,loopradius=loopRadius))
    lineSegs.extend(BuildArc(pi+capArc*0.5,pi*1.5-capArc*0.5,10,wireradius=wireradius,loopradius=loopRadius))
    lineSegs.extend(BuildArc(pi*1.5+capArc*0.5,2*pi-capArc*0.5,10,wireradius=wireradius,loopradius=loopRadius))
    frequency=298.06e6  # Hz
    Capacitance=3.3e-12 # Farads
    Y=AdmittanceSeriesZ(complex(0,-2*pi*frequency*Capacitance))
    


    
    
    Visualization(lineSegs)                             
    wiresegs=[WireSeg(l.startvect.xyz(),l.endvect.xyz(),radius=l.radius) for l in lineSegs]
    
    capacitors=[#NetworkSeg(0,1,Y),
                NetworkSeg(2,3,Y),
                NetworkSeg(4,5,Y),
                NetworkSeg(6,7,Y)
                ]

    mod=NECModel(geometry={"wiresegs":wiresegs,
                           "networksegs":capacitors},
                 excitation=VoltageSource(0),
                 freq=LinearFrequencyRange(freq=frequency),
                 radiationPattern=RadiationPatternSpherical(numTheta=37,numPhi=1,thetaS=pi/2,
                                                            thetaE=-pi/2,phiS=0,phiE=0,
                                                            gainNormalization="MAJOR",
                                                            gainType="DIRECTIVE",
                                                            gainAveraging=True))
    mod.runSim(printCardDeck=True,printRawResult=True)
    #print mod.RESULT.freqResults
    currents=[v[2] for v in  mod.RESULT.freqResults[0][1]["CURRENTS AND LOCATION"]]
    positions=[i for i in range(len(currents))]
    #print currents
    from KPlot import KBodePlot
    KBodePlot([(positions[8:],currents[8:],"currents")],linearfreqs=True,scaletype="value_linear")

def DualLoops(offsetFactor=3.0,wireradius=0.0005,
              loopradius=0.05,frequency=298.06e6,verticalOffsetFactor=4.0,verbose=False):
    from copy import deepcopy
    
    nip=pi/30.0
    numelements=16

    offset=offsetFactor*loopradius
    
    loop1=[Cylinder(Vector((loopradius,2*pi-nip,0),vtype="CYLINDRICAL"),
                    Vector((loopradius,     nip,0),vtype="CYLINDRICAL"),wireradius,color=(0,1,0))]
    loop1.extend(BuildArc(nip,2*pi-nip,numelements,loopradius=loopradius,wireradius=wireradius))
    loop2NumStart=len(loop1)
    loop2=[Cylinder(Vector((loopradius,2*pi-nip,0),vtype="CYLINDRICAL"),
                    Vector((loopradius,     nip,0),vtype="CYLINDRICAL"),wireradius,color=(0,1,0))]
    loop2.extend(BuildArc(nip,2*pi-nip,numelements,loopradius=loopradius,wireradius=wireradius))
    
    # Now let's move the loops
    loop1=[ls.translate(Vector((0,-offset*0.5,-wireradius*verticalOffsetFactor))) for ls in loop1]
    loop2=[ls.translate(Vector((0, offset*0.5, wireradius*verticalOffsetFactor))) for ls in loop2]

    objs=loop1
    objs.extend(loop2)
    

    wiresegs=[]
    for i,l in enumerate(objs):
        if i==0 or i==loop2NumStart:
            load=SimpleImpedanceLoad(complex(50.0,0.0))#SeriesRLCLoad(R=50)
            l=Cylinder(l.startvect,l.endvect,l.radius,color=(0,0,1))
            objs[i]=l
        else:
            load=None
        wiresegs.append(WireSeg(l.startvect.xyz(),l.endvect.xyz(),radius=l.radius,
                                load=load))
        
    if verbose: Visualization(objs)
        
    mod=NECModel(geometry={"wiresegs":wiresegs},
                 excitation=VoltageSource(0),
                 freq=LinearFrequencyRange(freq=frequency))
                 
    mod.runSim(printCardDeck=verbose,printRawResult=verbose)
    #print mod.RESULT.freqResults
    currents=[v[2] for v in  mod.RESULT.freqResults[0][1]["CURRENTS AND LOCATION"]]
    positions=[i for i in range(len(currents))]
    #print currents
    if verbose:
        from KPlot import KBodePlot
        KBodePlot([(positions,currents,"currents")],linearfreqs=True,scaletype="value_linear")
    return currents[:loop2NumStart],currents[loop2NumStart+1:]


def ComputeMultipleDualLoops():
    #DualLoops(verbose=True)
    #raise Exception
    from MRMaterials import Nuclei
    B0=1.0 # teslas
    loopRadius=0.05   # 0.01 m radius at 1.0 T gives a nice mutual inductance dip.
    wireradius=0.001
    verticaloffsetFactor=10.0


    Proton=Nuclei["1_H"]
    frequency=Proton.frequencyHz(B0) # B0 in teslas
    print "frequency:",frequency*1e-6,"MHz"
    
    allCurrentsloop1=[]
    allCurrentsloop2=[]
    plots=[]
    dOffset=0.1
    offsetFactors=arange(0.0,3+dOffset,dOffset)
    #print offsetFactors

    #def _worker(args):
    #    oF=args[0]
    #    currentsloop1,currentsloop2=DualLoops(offsetFactor=oF,loopradius=loopRadius,wireradius=wireradius,
    #                                          frequency=frequency,verbose=False,verticalOffsetFactor=verticaloffsetFactor)
    #    return currentsloop1,currentsloop2

    #from KFarm import DoPool
    #stuff=DoPool(_worker,[[o] for o in offsetFactors],[],resultlen=2)

    for oF in offsetFactors:
        currentsloop1,currentsloop2=DualLoops(offsetFactor=oF,loopradius=loopRadius,wireradius=wireradius,
                                             frequency=frequency,verbose=False,verticalOffsetFactor=verticaloffsetFactor)
        allCurrentsloop1.append(currentsloop1)
        allCurrentsloop2.append(currentsloop2)
        positions=[i for i in range(len(currentsloop1))]
        #plots.append((positions,currentsloop1,"oF:"+repr(oF)))
        #print currents
        #print positions
        print oF
    
    result=[B0,frequency,loopRadius,wireradius,dOffset,verticaloffsetFactor,allCurrentsloop1,allCurrentsloop2]
    import cPickle
    fname="out.cPickle"
    fout=open(fname,"w+")
    cPickle.dump(result,fout)
    fout.close()

    PlotDualLoopResult(fname=fname)
    #from KPlot import KBodePlot
    #KBodePlot(allCurrents,linearfreqs=True,scaletype="value_linear")

    

def PlotDualLoopResult(fname=None):
    if fname==None:
        import sys
        fname=sys.argv[1]
    import cPickle
    fin=open(fname,"r")
    result=cPickle.load(fin)
    fin.close()
    B0,frequency,loopRadius,wireradius,dOffset,verticaloffsetFactor,allCurrentsloop1,allCurrentsloop2=result
    print frequency

    def myfunc1(x,y):
        return abs(allCurrentsloop1[int(x)][int(y)])

    def myfunc2(x,y):
        return abs(allCurrentsloop2[int(x)][int(y)])

    from KPlot import kplot
    print len(allCurrentsloop1),len(allCurrentsloop1[0])
    #kplot(myfunc1,          0,len(allCurrentsloop1)-1,1,0,len(allCurrentsloop1[0])-1,1,type="3D")
    print len(allCurrentsloop2),len(allCurrentsloop2[0])
    kplot(myfunc2,          0,len(allCurrentsloop2)-1,1,0,len(allCurrentsloop2[0])-1,1,type="3D")


if __name__ == '__main__': 
    import sys
    if len(sys.argv)!=1:
        PlotDualLoopResult()
    else:
        ComputeMultipleDualLoops()

