#!/usr/bin/env python
from KStandard import attributesFromDict
from KGeomVisualization import *
from numpy import *
import commands
from math import *
import fileaccess
#
# A Python interface to the NECcpp code.
#
#

class Excitation():
    """
    Page 48 of the manual.
    """
    def __init__(self):
        attributesFromDict(locals())

    def Generate(self):
        """
        returns the EX line
        """
        return "EX"

class VoltageSource(Excitation):
    """
    page 48
    """
    def __init__(self,
                 segnum,     # I2
                 segsplit=1, # I3
                 extype="APPLIED_EFIELD",
                 calcMaxAdmittanceMatrixAssym=0, # I4 first param
                 calcImpedanceTable=1,           # I4 second param
                 VReal=1.0,
                 VImag=0.0,
                 ImpedanceNormalization=None
                 ):
        """
        APPLIED_EFIELD or CURRENTSLOPE
        """
        attributesFromDict(locals())
    def Generate(self):
        result=Excitation.Generate(self)
        if self.extype=="APPLIED_EFIELD":
            result=result+" 0"
        elif self.extype=="CURRENTSLOPE":
            result=result+" 5"
        else:
            raise Exception
        result=result+" "+repr(self.segnum)
        result=result+" "+repr(self.segsplit)
        result=result+" "+repr(self.calcMaxAdmittanceMatrixAssym)
        result=result+repr(self.calcImpedanceTable)
        result=result+" "+repr(self.VReal)
        result=result+" "+repr(self.VImag)
        if self.calcImpedanceTable==1 and self.ImpedanceNormalization!=None:
            result=result+" "+repr(self.ImpedanceNormalization)
        result=result+"\n"
        return result

class IncidentWaveSource(Excitation):
    """
    page 48
    """
    def __init__(self,extype="LINEAR"):
        """
        LINEAR,RIGHTELLIPTIC,LEFTELLIPTIC
        """
        attributesFromDict(locals())

class CurrentSource(Excitation):
    """
    page 48
    """
    def __init__(self):
        attributesFromDict(locals())
    
class WireConductivityLoad():
    """
    page 58 of the manual
    """
    def __init__(self,
                 conductivity, # mhos/meter
                 LDTAGF=0,LDTAGT=None
                 ):
        attributesFromDict(locals())
    def Generate(self,num):
        result="LD"+" 5" #Five for wire conductivity
        result=result+" "+repr(num)+" "+repr(self.LDTAGF)+" "
        if self.LDTAGT!=None:
            result=result+repr(self.LDTAGT)
        else:
            result=result+repr(self.LDTAGF) # Again
        result=result+" "+repr(self.conductivity)
        result=result+"\n"
        return result


class SimpleImpedanceLoad():
    """
    page 58 of the manual
    """
    def __init__(self,
                 impedance, # complex number (ohms)
                 LDTAGF=0,LDTAGT=None
                 ):
        attributesFromDict(locals())
    def Generate(self,num):
        result="LD"+" 4" #Impedance
        result=result+" "+repr(num)+" "+repr(self.LDTAGF)+" "
        if self.LDTAGT!=None:
            result=result+repr(self.LDTAGT)
        else:
            result=result+repr(self.LDTAGF) # Again
        result=result+" "+repr(self.impedance.real)+" "+repr(self.impedance.imag)
        result=result+"\n"
        return result


class SeriesRLCLoad(WireConductivityLoad):
    """
    page 58 of the manual
    """
    def __init__(self,
                 R=0, # Resistance Ohms
                 L=0, # Inductance Henries
                 C=0, # Capacitance Farads
                 LDTAGF=0,LDTAGT=None):
        attributesFromDict(locals())
    def Generate(self,num):
        result="LD"+" 0"#Zero for SeriesRLC
        result=result+" "+repr(num)+" "+repr(self.LDTAGF)+" "
        if self.LDTAGT!=None:
            result=result+repr(self.LDTAGT)
        else:
            result=result+repr(self.LDTAGF) # Again
        result=result+" "+repr(self.R)+" "+repr(self.L)+" "+repr(self.C)
        result=result+"\n"
        return result

class LinearFrequencyRange():
    """
    Page 52

    We probably don't want to use this. Just put in numerical frequencies and
    run the loop in python...
    """
    def __init__(self,freq=299.8e6,dfreq=1e6,num=1):
        attributesFromDict(locals())
    def Generate(self):
        result="FR 0 "+repr(self.num)+" 0 0"
        result=result+" "+repr(self.freq*1e-6)+" "+repr(self.dfreq*1e-6)+"\n"
        return result

class MultiplicativeFrequencyRange(LinearFrequencyRange):
    """
    Page 52

    We probably don't want to use this. Just put in numerical frequencies and
    run the loop in python...
    """
    def __init__(self,freq=299.8e6,dfreq=1.1,num=1):
        attributesFromDict(locals())
    def Generate(self):
        result="FR 1 "+repr(self.num)+" 0 0"
        result=result+" "+repr(self.freq*1e-6)+" "+repr(self.dfreq)+"\n"
        return result

class WireSeg():
    """
    page 26
    """
    def __init__(self,xyzStart,xyzEnd,
                 radius=0.001,          
                 numSplit=1,            # I2
                 load=None              # A load as on page 58 of the manual.
                 ):
        """
        The tag number is calculated automatically in NECModel
        """
        attributesFromDict(locals())
    
    def set_load(self,load):
        self.load=load

    def Generate(self,num):
        result="GW "+repr(num)+","+repr(self.numSplit)
        result=result+","+repr(self.xyzStart[0])+","+repr(self.xyzStart[1])+","+repr(self.xyzStart[2])
        result=result+","+repr(self.xyzEnd[0])+","+repr(self.xyzEnd[1])+","+repr(self.xyzEnd[2])
        result=result+","+repr(self.radius)+"\n"
                     
        return result

class TaperedWireSeg(WireSeg):
    """
    page 26. Tapered wire not yet implemented.
    """
    pass

class NetworkSeg():
    """
    page 62
    """
    def __init__(self,Port1TAG,Port2TAG,YMatrix,Port1SEG=1,Port2SEG=1):
        attributesFromDict(locals())
    def Generate(self):
        if self.YMatrix[0][1] != self.YMatrix[1][0]: raise Exception
        result="NT"
        result=result+" "+repr(self.Port1TAG)+" "+repr(self.Port1SEG)
        result=result+" "+repr(self.Port2TAG)+" "+repr(self.Port2SEG)
        result=result+" "+repr(self.YMatrix[0][0].real)+" "+repr(self.YMatrix[0][0].imag)
        result=result+" "+repr(self.YMatrix[0][1].real)+" "+repr(self.YMatrix[0][1].imag)
        result=result+" "+repr(self.YMatrix[1][1].real)+" "+repr(self.YMatrix[1][1].imag)
        result=result+"\n"
        return result

class TransmissionLineSeg(NetworkSeg):
    """
    page 73
    """
    def __init__(self,Port1TAG,Port2TAG,Port1SEG=1,Port2SEG=1,Z0=50.0,length=None,
                 Port1ShuntAdmittance=complex(0,0),Port2ShuntAdmittance=complex(0,0)):
        attributesFromDict(locals())
    def Generate(self):
        result="TL"
        result=result+" "+repr(self.Port1TAG)+" "+repr(self.Port1SEG)
        result=result+" "+repr(self.Port2TAG)+" "+repr(self.Port2SEG)
        result=result+" "+repr(self.Z0)+" "
        if self.length!=None:
            result=result+repr(self.length)
        result=result+","+repr(self.Port1ShuntAdmittance.real)+","+repr(self.Port1ShuntAdmittance.imag)
        result=result+","+repr(self.Port2ShuntAdmittance.real)+","+repr(self.Port2ShuntAdmittance.imag)
        result=result+"\n"
        return result


class RadiationPatternSpherical():
    """
    page 69
    """
    def __init__(self,
                 thetaS=0,thetaE=pi,numTheta=8,
                 phiS=0,phiE=2*pi,numPhi=8,
                 r=None,
                 outputFormat="VERTHORIZ", #VERTHORIZ or MAJORMINOR                   | X in XNDA
                 gainNormalization="NONE", #NONE, MAJOR, MINOR, VERT, HORIZ, or TOTAL | N in XNDA
                 gainType="POWER",         #POWER or DIRECTIVE                        | D in XNDA
                 gainAveraging=False,      #True or False. We forget A=2.             | A in XNDA 
                 ):
        dTheta=(thetaE-thetaS)/float(numTheta-1)
        dPhi=(phiE-phiS)/float(numPhi)   # No minus one because of the wrap around.
        attributesFromDict(locals())
    def Generate(self):
        result="RP"
        result=result+" "+"0"
        result=result+" "+repr(self.numTheta)+" "+repr(self.numPhi)+" "
        if self.outputFormat=="VERTHORIZ":
            result=result+"1"
        elif self.outputFormat=="MAJORMINOR":
            result=result+"0"
        else:
            raise Exception
        if self.gainNormalization=="NONE":
            result=result+"0"
        elif self.gainNormalization=="MAJOR":
            result=result+"1"
        elif self.gainNormalization=="MINOR":
            result=result+"2"
        elif self.gainNormalization=="VERT":
            result=result+"3"
        elif self.gainNormalization=="HORIZ":
            result=result+"4"
        elif self.gainNormalization=="TOTAL":
            result=result+"5"
        else:
            raise Exception
        if self.gainType=="POWER":
            result=result+"0"
        elif self.gainType=="DIRECTIVE":
            result=result+"1"
        else:
            raise Exception
        if self.gainAveraging:
            result=result+"0"
        else:
            result=result+"1"
        result=result+" "+repr(self.thetaS*180/pi)+" "+repr(self.phiS*180/pi)
        result=result+" "+repr(self.dTheta*180/pi)+" "+repr(self.dPhi*180/pi)
        if self.r!=None:
            result=result+" "+repr(self.r)
        result=result+"\n"
        return result

class NECModel():
    """
    Some geometry to be modelled.
    """
    def __init__(self,description="My Model\n Is here.",
                 geometry={"wiresegs":[]   # A straight list of the wiresegments
                           },
                 excitation=None,
                 freq=LinearFrequencyRange(),
                 radiationPattern=None,
                 computeCharges=True):
        """
        The geometry is a list of KGeomVisualization Objects.
        These objects can be viewed and can be simulated.
        """
        attributesFromDict(locals())
        
    def set_freq(self,freq):
        self.freq=freq
    def set_description(self,descrip):
        self.description=descrip
    def set_excitation(self,exc):
        self.excitation=exc

    def _genTMPFILE(self):
        return "./out"

    def runSim(self,printCardDeck=False,printRawResult=False):
        CardDeck=self.ProduceCardDeck()
        if printCardDeck:
            print CardDeck
        deckF=self._genTMPFILE()+".txt"
        resultF=self._genTMPFILE()+".out.txt"
        fout=open(deckF,"w+")
        fout.write(CardDeck)
        fout.close()
        command="/home/edlerk/Projects/ElectronicsUtilities/RF/NEC_MOM/python-necpp/nec2cpp -i "+deckF+" -o "+resultF

        commands.getoutput(command)
        fin=open(resultF,"r")
        self.resultLines=fileaccess.getlineslist(resultF)
        if printRawResult:
            for l in self.resultLines: print l,
        self.RESULT=NECResult(self.resultLines)
        #print NR
    def ProduceCardDeck(self):
        """
        Returns a string representing a necpp card deck.
        """
        result=""
        # First we put in the description as comments 
        # (page 15 of the manual)
        descriplines=self.description.split('\n')
        for dl in descriplines[:-1]:
            result=result+"CM "+dl+"\n"
        result=result+"CE "+descriplines[-1]+"\n"
        
        # Read numerical green's function file "GF"
        # (page 19 and 78)

        # Wire line segments
        # (page 26 of the manual)
        # Note that there is a radius tapering option that 
        # is not yet implemented here.
        for i,ws in enumerate(self.geometry["wiresegs"]):
            result=result+ws.Generate(i)

        # Wire Arcs (page 16 of the manual)

        # Wire helix/spiral (page 20 of the manual)

        # Coordinate transformation (page 21 of the manual)
        # I am not sure how useful this will be for my needs.

        # Generate Cylindrical Structure 
        # (page 23 of the manual)

        # Reflection in Coordinate Planes "GX" (page 27)

        # Surface Patch
        # (page 30 of the manual)

        # Multiple Patch Surface
        # (page 34 of the manual)

        # Scale Geometry
        # (page 25 of the manual)
        #result=result+"GS\t"+repr(self.scaleFactor)+"\n"

        # End Geometry (page 17 of the manual)
        result=result+"GE\n"

        # Now for the analysis stuff.
        # (Program control starting on page 43 of the manual)

        # Maximum coupling "CP" (page 45)

        # Extended Thin Wire Kernel "EK" (page 46)

        
        # Loading "LD" (page 58)
        for i,ws in enumerate(self.geometry["wiresegs"]):
            if ws.load!=None:
                result=result+ws.load.Generate(i)

        # Networks "NT" (page 62)
        if self.geometry.has_key("networksegs"):
            for nt in self.geometry["networksegs"]:
                result=result+nt.Generate()

        # Transmission Line "TL" (page 73)
        if self.geometry.has_key("transmissionlines"):
            for tl in self.geometry["transmissionlines"]:
                result=result+tl.Generate()

        # Excitation "EX" (page 48)
        result=result+self.excitation.Generate()
        
        
        # Frequency "FR" (page 52)
        if isinstance(self.freq,LinearFrequencyRange):
            result=result+self.freq.Generate()
        else:
            f=LinearFrequencyRange(freq=self.freq)
            result=result+f.Generate()

        # Next Structure "NX" (page 65)
        # I don't think we care about this.

        # Interaction Approximation Range "KH" (page 57)

        # Ground Parameteres "GN" (page 55)

        # Additional Ground Parameters "GD" (page 53)

        # Print control for charge on wires "PQ" (page 66)
        # Print Charge densities on everything
        if self.computeCharges: result=result+"PQ 0\n"

        # Near Fields "NE" "NH" (page 60)

        # Radiation Pattern "RP" (page 69)
        if self.radiationPattern!=None:
            result=result+self.radiationPattern.Generate()

        # Write NGF File "WG" (page 75)

        # Execute "XQ" (page 76)
        # I don't quite understand what this does yet
        result=result+"XQ\n"


        # End of run "EN" (page 46)
        result=result+"EN\n"
        return result

class NECResult():
    """
    The result of running a simulation. 
    There may be several sorts of results burried
    in this data structure.
    """
    def _buildFrequencyResult(self,i,resultLinesList,verbose=False):
        result={}
        numLines=len(resultLinesList)
        l=resultLinesList[i]
        freq=float( l.strip().split()[1] )*1e6 #Hz
        if verbose: print "FREQ:",freq
        tmp=""
        while(i<numLines and tmp!="--------- FREQUENCY --------"):
            while(tmp[0:4]!="----" and i<numLines):
                l=resultLinesList[i]
                tmp=l.strip()
                i=i+1
            # We found something
            
            if tmp=="--------- FREQUENCY --------": break
            if tmp=="--------- ANTENNA INPUT PARAMETERS ---------":
                if verbose: print "Found:",tmp
                # We found the antenna input parameters
                i=i+2
                AIPLines={}
                l=resultLinesList[i]
                tmp=tmp=l.strip()
                while(tmp[0:4]!="----"):
                    line=tmp.split()
                    TAG=line[0]
                    SEG=line[1]
                    #print line
                    VOLTAGE=complex(float(line[2]),float(line[3]))
                    CURRENT=complex(float(line[4]),float(line[5]))
                    IMPEDANCE=complex(float(line[6]),float(line[7]))
                    ADMITTANCE=complex(float(line[8]),float(line[9]))
                    POWER=float(line[10])
                    AIPLines["TAG{"+TAG+"}:SEG{"+SEG+"}"]={"VOLTAGE":VOLTAGE,"CURRENT":CURRENT,"IMPEDANCE":IMPEDANCE,"ADMITTANCE":ADMITTANCE,"POWER":POWER}
                    i=i+1
                    l=resultLinesList[i]
                    tmp=l.strip()
                    
                result["ANTENNA INPUT PARAMETERS"]=AIPLines
            if tmp=="-------- CURRENTS AND LOCATION --------":
                if verbose: print "Found:",tmp
                i=i+4
                l=resultLinesList[i]
                tmp=l.strip()
                CURRENTRESULTS=[]
                while(tmp[0:4]!="----"):
                    line=tmp.split()
                    #print line
                    TAG=line[1]
                    SEG=line[0]
                    XC=float(line[2])
                    YC=float(line[3])
                    ZC=float(line[4])
                    LL=float(line[5])
                    CUR=complex(float(line[6]),float(line[7]))
                    CURRENTRESULTS.append([TAG,SEG,CUR])
                    i=i+1
                    l=resultLinesList[i]
                    tmp=l.strip()
                    
                result["CURRENTS AND LOCATION"]=CURRENTRESULTS
            if tmp=="------ CHARGE DENSITIES ------":
                if verbose: print "Found:",tmp
                i=i+4
                l=resultLinesList[i]
                tmp=l.strip()
                CHARGERESULTS=[]
                while(tmp[0:4]!="----"):
                    line=tmp.split()
                    #print line
                    TAG=line[1]
                    SEG=line[0]
                    XC=float(line[2])
                    YC=float(line[3])
                    ZC=float(line[4])
                    LL=float(line[5])
                    CHARGE=complex(float(line[6]),float(line[7]))
                    CHARGERESULTS.append([TAG,SEG,CHARGE])
                    i=i+1
                    l=resultLinesList[i]
                    tmp=l.strip()
                    
                result["CURRENTS AND LOCATION"]=CURRENTRESULTS
            if tmp=="---------- RADIATION PATTERNS -----------":
                if verbose: print "Found:",tmp
                i=i+4
                l=resultLinesList[i]
                tmp=l.strip()
                RADPATTERN=[]
                while(tmp[0:4]!="----"):
                    line=tmp.split()
                    #print line
                    THETA=float(line[0])*pi/180
                    PHI=float(line[1])*pi/180
                    XTHETAM=float(line[8])
                    XTHETAP=float(line[9])*pi/180
                    XPHIM=float(line[10])
                    XPHIP=float(line[11])*pi/180
                    
                    RADPATTERN.append([THETA,PHI,
                                       complex(XTHETAM*cos(XTHETAP),XTHETAM*sin(XTHETAP)),
                                       complex(XPHIM*cos(XPHIP),XPHIM*sin(XPHIP)),
                                       ])
                    i=i+1
                    l=resultLinesList[i]
                    tmp=l.strip()
                    
                result["RADIATION PATTERN"]=RADPATTERN
            if tmp=="---------- POWER BUDGET ---------":
                PB={}
                i=i+1
                if verbose:
                    print "Found:",tmp
                l=resultLinesList[i]
                tmp=l.strip()
                #print tmp.split()
                PB["INPUT POWER"]=float(tmp.split()[-2])
                i=i+1
                l=resultLinesList[i]
                tmp=l.strip()
                #print tmp.split()
                PB["RADIATED POWER"]=float(tmp.split()[-2])
                i=i+1
                l=resultLinesList[i]
                tmp=l.strip()
                #print tmp.split()
                PB["STRUCTURE LOSS"]=float(tmp.split()[-2])
                i=i+1
                l=resultLinesList[i]
                tmp=l.strip()
                #print tmp.split()
                PB["NETWORK LOSS"]=float(tmp.split()[-2])
                i=i+1
                l=resultLinesList[i]
                tmp=l.strip()
                #print tmp.split()
                PB["EFFICIENCY"]=float(tmp.split()[-2])*0.01
                result["POWER BUDGET"]=PB
                #print PB
            if i>=numLines: break
            i=i+1
            l=resultLinesList[i]
            tmp=l.strip()
            
        self.freqResults.append([freq,result])
        return i



    def __init__(self,resultLinesList):
        #attributesFromDict(locals())
        # First we look for the comments
        numLines=len(resultLinesList)
        i=0
        tmp=""
        while(tmp!="---------------- COMMENTS ----------------" and i<numLines):
            l=resultLinesList[i]
            tmp=l.strip()
            i=i+1
        # We must have found the comments
        #print "Found comments",i
        comments=""
        tmp=""
        while(tmp!="-------- STRUCTURE SPECIFICATION --------" and i<numLines):
            comments=comments+tmp+"\n"
            l=resultLinesList[i]
            tmp=l.strip()
            i=i+1
        #print "######################"
        #print comments
        self.comments=comments
        # We have now found the structure specification
        structureSpec=""
        tmp=""
        l=""
        while(tmp!="--------- FREQUENCY --------" and i<numLines):
            structureSpec=structureSpec+l
            l=resultLinesList[i]
            tmp=l.strip()
            i=i+1
        #print "######################"
        #print structureSpec
        self.structureSpec=structureSpec
        
        # Now we should enter a loop looking at the different frequencies.
        self.freqResults=[]
        while(i<numLines):
            i=self._buildFrequencyResult(i,resultLinesList)




if __name__ == '__main__':
    pass
    
    #KarlLoop()
    #DualLoops()
    #PlotDualLoopResult()
