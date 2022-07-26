from PathFinder.output import CloneFrequencyAnalizer
from PathFinder.output.CloneFrequencyAnalizer import CloneFrequencyAnalizer
#from output.CloneFrequencyAnalizer import CloneFrequencyAnalizer
#from output.CloneFrequencyAnalizer import remove_minor_clone
from PathFinder.alignments.MegaAlignment import MegaAlignment
import sys
import os

Align=MegaAlignment()
MAO='infer_ML_nucleotide.mao'
CF=CloneFrequencyAnalizer()
ItMeg=sys.argv[1]   ###CloneFinder output meg file
TuAnno=sys.argv[2]  ## tumor annotation table, e.g., m8Mseed76CFin_TuAnno.txt (input to the pipeline)/ "P" is primary, or adjust the PathFinder command line. 
CloFre=ItMeg[:-4]+'.txt' ## CloneFinder output txt file (clone frequency)
Cut=0.05  ## minimum inferred clone frequency to be present. Can you allow users to change the cutoff? default is 0.02
cwd=os.getcwd()


CF.remove_minor_clone(CloFre,Cut,ItMeg) 
In=open(ItMeg,'r').readlines()
CloLs,Clo2Seq=Align.name2seq(In)
Samp2Clo,HitCloLs,T2C2F=CF.GetCloHitForTu( ItMeg[:-4]+'_'+str(Cut)+'.txt',0)
Samp2Tu={}
NewT2C2F={}
TuAnno=open(TuAnno,'r').readlines()[1:]
for i in TuAnno:
    i=i.strip().split('\t')
    Samp=i[0].replace('T-', '')
    Tu=i[1].replace('T-', '')
    Samp2Tu[Samp]=Tu
    NewT2C2F[Tu]={}	
	

print (Samp2Tu)	
print (Samp2Clo)
TotCloLs=[]
for Samp in Samp2Clo:
    Samp = Samp
    Tu=Samp2Tu[Samp]
    CloLsS=Samp2Clo[Samp]
    for Clo in CloLsS:
       NewT2C2F[Tu][Clo.replace('#','')]=1
    TotCloLs+=CloLsS         
print (NewT2C2F)	   
CF.UpCloFreqTa1(NewT2C2F,HitCloLs,ItMeg[:-4]+'_'+str(Cut)+'.txt')	 

Align.remove_identical_clones(Clo2Seq,NewT2C2F,ItMeg[:-4]+'_PathFinder_'+str(Cut))
CloLs,Clo2Seq=Align.name2seq(open(ItMeg[:-4]+'_PathFinder_'+str(Cut)+'.meg','r').readlines())
out=''
for Clo in Clo2Seq:
    if TotCloLs.count(Clo.replace('#',''))!=0:
       out+='>'+Clo.replace('#','')+'\n'+Clo2Seq[Clo]+'\n'
out+='>Normal\n'+('A'*len(Clo2Seq[Clo]))+'\n'	
OutF=open(ItMeg[:-4]+'_PathFinder_'+str(Cut)+'.fas','w')
OutF.write(out)
OutF.close()	
