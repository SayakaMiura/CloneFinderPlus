from resample.generate_replicate import generate_replicate
from alignments.MegaAlignment import MegaAlignment
import pandas as pd
import sys
import os
import shutil
import glob
import numpy as np
import pandas as pd
from scipy import stats
import Functions3

CFin=sys.argv[1]
TuAnno=sys.argv[2]
BooC=int(sys.argv[3])
 


def Clean(Tar):
    Ls=glob.glob(Tar)
    for i in Ls: 
       os.remove(i)
 
Clean('FullData*')
if os.path.exists(TuAnno)!=True: 
   print (TuAnno+' not found')
elif os.path.exists(CFin)!=True: 
   print (CFin+' not found')  
DriMap='n'
if DriMap=='y':
    ID=CFin.split('\\')[-1][:-7]
    DriMet2File={}
    DriMet2File['CHASMpri']='C:\\Users\\kumarlab\\Desktop\\CloneFinderPlus\\Empirical\\CRAVATout\\Townsend\\AnyCancer\\'+ID+'_CRAVATout.txt'
    DriMet2File['PriOncodriverMUT']='C:\\Users\\kumarlab\\Desktop\\CloneFinderPlus\\Empirical\\CGIout\\AnyCancer\\'+ID+'_CGIout.txt'	
    PreType='AnyCancer'#'SpecificCancer'
    DriverLabDic={'CHASMpri':'6,0.05,Cut','PriOncodriverMUT':'4,driver:passenger:non-protein_affecting'}
Wdir=os.getcwd()
CFPpath=Wdir+os.sep+'CloneFinderPlus'+os.sep
CFPcommand='RunCloneFinder.py'
PPmigProCut=0
python='python'
Rscript='Rscript'
CloneMinDiff=1
CloneBooCut=0.1
PPoutdir=CFPpath+os.sep+'PathFinderRes'
Summary=''
Align=MegaAlignment()
shutil.copy2(CFin, 'FullData.txt')
print ('Generate clone prediction')
Summary+="ClonePrediction\tyes\n"

if os.path.exists(CFPpath+CFPcommand)!=True:
        Summary+=CFPcommand+" should be at "+CFPpath+'\n'
else:
        os.chdir(CFPpath)
        FileLs=glob.glob('FastCloneRes'+os.sep+'*')	
        for i in FileLs:
            InLs=glob.glob(i+os.sep+'*')
            for ii in InLs: os.remove(ii)	
        FileLs=glob.glob('Data'+os.sep+'FastCloneData'+os.sep+'*.tsv')		
        for i in FileLs: os.remove(i)		
        FileLs=glob.glob('PathFinderRes'+os.sep+'*.txt')+glob.glob('PathFinderRes'+os.sep+'*.png')+glob.glob('PathFinderRes'+os.sep+'*.dot')
        for i in FileLs: os.remove(i)
        FileLs=glob.glob('PathFinderRes'+os.sep+'*')	
        for i in FileLs:
            InLs=glob.glob(i+os.sep+'*')
            for ii in InLs: os.remove(ii)	         
        os.system(python+" "+CFPcommand+" "+Wdir+'\\FullData.txt '+TuAnno)	
        print (python+" "+CFPcommand+" "+Wdir+'\\FullData.txt '+TuAnno)
        MigInf=glob.glob('PathFinderRes'+os.sep+'*_Mig.txt')	
        MigInf.append('NA')
        print("******************" + MigInf[0]) 
        print(os.getcwd())
        PFout=MigInf[0]	
        AncSeqTaLs=glob.glob('PathFinderRes\\ancestral_inference_logging\\nuc_anc_seqs_*_eps.csv')		
        if os.path.exists(MigInf[0])==True:
            if DriMap=='y':		 
                for DriPreMet in DriMet2File:
                    DriAnno=DriMet2File[DriPreMet]
                    DriverLab=DriverLabDic[DriPreMet]	
                    MapDrCom=	python+' MapDriverCount.py '+DriverLab+' '+DriAnno+' '+PFout+' '+CFin+' '+AncSeqTaLs[0]+ ' '+DriPreMet+'_'+PreType		
                    print (MapDrCom)					
                    os.system(MapDrCom)
                    OutF=PFout[:-4]+'_'+DriPreMet+'_'+PreType+'.txt'
                    print (OutF)									
                    shutil.copy2(OutF,os.path.join(os.path.dirname(os.getcwd()),CFin[:-4]+'_'+DriPreMet+'_'+PreType+'_Mig.txt')) 
            else: shutil.copy2(MigInf[0],CFin[:-4]+'_Mig.txt')          		
        os.chdir(Wdir)				
		

print (Summary)
print ('generate bootstrap replicates')
print (BooC)
BooG=generate_replicate()
print ('generate inference on each bootstrap replicates')	
BooTreeFileLs=[]
FailLs=[]
outCloFre=['Rep\tTumor\tRef clone\trep clone frequency']
MigBooOut=[]
if os.path.exists(CFPpath+CFPcommand)!=True:
        Summary+=CFPcommand+" should be at "+CFPpath+'\n'
else:		
        os.chdir(CFPpath)	
        Boo=1
        PFsuccessC=0
        Mig2C={}		
        while Boo<=BooC:
           BooG.sample_read_1rep(Wdir+'\\FullData.txt',Boo)			
           FileLs=glob.glob('FastCloneRes'+os.sep+'*')	
           for i in FileLs:
               InLs=glob.glob(i+os.sep+'*')
               for ii in InLs: os.remove(ii)
           FileLs=glob.glob('Data'+os.sep+'FastCloneData'+os.sep+'*.tsv')		
           for i in FileLs: os.remove(i)	
           FileLs=glob.glob('PathFinderRes'+os.sep+'*.txt')+glob.glob('PathFinderRes'+os.sep+'*.png')+glob.glob('PathFinderRes'+os.sep+'*.dot')
           for i in FileLs: os.remove(i)
           FileLs=glob.glob('PathFinderRes'+os.sep+'*')	
           for i in FileLs:
               InLs=glob.glob(i+os.sep+'*')
               for ii in InLs: os.remove(ii)		
                           
           os.system(python+" "+CFPcommand+" "+Wdir+'\\FullData_Rep'+str(Boo)+'.txt '+TuAnno)
           if os.path.exists(Wdir+'\\FullData_Rep'+str(Boo)+'snv_CloneFinderPlus.meg')!=True:
               FailLs.append(Boo)	
           MigInf=glob.glob('PathFinderRes'+os.sep+'*_Mig.txt')	
           MigInf.append('NA')	
           if os.path.exists(MigInf[0])==True:
               BooMig=open(MigInf[0],'r').readlines()[1:]
               for Path in BooMig:
                    Path=Path.split('\t')[0]	
                    Mig2C[Path]=Mig2C.get(Path,0)+1					
               PFsuccessC+=1	             			   
               PFout=MigInf[0]	
               AncSeqTaLs=glob.glob('PathFinderRes\\ancestral_inference_logging\\nuc_anc_seqs_*_eps.csv')		
               if DriMap=='y':	
			   
                for DriPreMet in DriMet2File:
                    DriAnno=DriMet2File[DriPreMet]
                    DriverLab=DriverLabDic[DriPreMet]	
                    MapDrCom=	python+' MapDriverCount.py '+DriverLab+' '+DriAnno+' '+PFout+' '+CFin+' '+AncSeqTaLs[0]+ ' '+DriPreMet+'_'+PreType		
                    print (MapDrCom)					
                    os.system(MapDrCom)
                    OutF=PFout[:-4]+'_'+DriPreMet+'_'+PreType+'.txt'
                    print (OutF)					
                    MigOut=open(OutF,'r').readlines()[1:]
                    for i in MigOut:
                       i=i.split('\t')
                       Ed=i[0].split('[')[0]
                       if Ed[0]=='P': Type='P->M'
                       elif Ed.split('->')[1][0]=='P': Type='M->P'
                       else: Type='M->M'					   
                       MigBooOut+=[str(Boo)+'\t'+DriPreMet+'\t'+Ed+'\t'+Type+'\t'+i[4]+'\t'+i[6]+'\n']								
               Boo+=1
           if os.path.exists(TuAnno)!=True:                 
               if os.path.exists(Wdir+'\\FullData_Rep'+str(Boo)+'snv_CloneFinderPlus.meg')==True: 
                    print (Wdir+'\\FullData_Rep'+str(Boo)+'snv_CloneFinderPlus.meg')               
                    Boo+=1

        os.chdir(Wdir)
if DriMap=='y':
    out='RepID\tMethod\tPath\tPathType\tDriverCount\tTotalMutCount\n'+''.join(MigBooOut)
    Functions3.GetOut(CFin[:-4]+'_AllMig.txt',out)
	
if Mig2C!={}:
    print ('score migration paths')			
    out='Edge\tBootstrap Confidence\n'
    for Mig in Mig2C:
       Bc=Mig2C[Mig]
       BcPro=1.0*Bc/PFsuccessC	   
       out+=Mig+'\t'+str(BcPro)+'\n'
    Align.GetOut(CFin[:-4]+'_Mig_consensus.txt',out)
    os.system(python+' MakeDot.py '+CFin[:-6]+'_Mig_consensus.txt')	
print (PFsuccessC)
print ('score and aggregate for clone')
if BooC>0:
    os.system(python+' CloneGenoPreAbsBooSum1plus.py '+str(BooC)) 
if os.path.exists('FullDatasnv_CloneFinderPlus.meg')==True:
    shutil.copy2('FullDatasnv_CloneFinderPlus.meg',CFin[:-4]+'_CloneFinderPlus.meg')
    shutil.copy2('FullDatasnv_CloneFinderPlus.txt',CFin[:-4]+'_CloneFinderPlus.txt')
if os.path.exists('FullDatasnv_CloneFinderPlus_All.fasta')==True:
    shutil.copy2('FullDatasnv_CloneFinderPlus_All.fasta',CFin[:-4]+'_All.fasta')
    shutil.copy2('FullDatasnv_CloneFinderPlus_All.txt',CFin[:-4]+'_All.txt')

    Functions3.MakeConsenGeno(CFin[:-4]+'_All.fasta',CloneMinDiff,CloneBooCut,BooC) #1 (diff),0.5,100 (boo count),1 (diff),0,100 (boo count) 
if DriMap=='y':
   ConGen=CFin[:-4]+'_All_consen.txt'
   MetOr=list(DriMet2File.keys())  
   Pos2Pre={}
   Pos2Gene={}   
   for Met in MetOr:
      Anno=CFin[:-4]+'_'+Met+'_AnyCancer.txt'
      Anno=open(Anno,'r').readlines()
      Last=Anno[0].split('\t')[-1].strip()
      if Last=='Gene': Col=-2	
      else: Col=-1	  
      Anno=Anno[1:]
      Pos=1	  
      for i in Anno:
         i=i.split('\t')
         Mut=i[Col].strip()
         Pos2Pre[Pos]=Pos2Pre.get(Pos,[])+[Mut]		 
         if Last=='Gene': Pos2Gene[Pos]=i[-1].strip()
         Pos+=1	
   ConGen=open(ConGen,'r').readlines()
   out=ConGen[0].strip()+ '\t'+'\t'.join(MetOr)+'\tConsensus\tGene\n'
   ConGen=ConGen[1:]
   for i in ConGen:
      if i[0]=='>': out+=i[1:].strip()+'\t'
      else: out+=i.strip()+'\t'	  
      Pos=int(i.split('\t')[1])
      PreLs=Pos2Pre[Pos]
      if PreLs==['driver','driver']: Consen='driver'
      else: Consen='Nondriver'
      out+='\t'.join(PreLs)+'\t'+Consen+'\t'+Pos2Gene[Pos]+'\n'
   Functions3.GetOut(CFin[:-4]+'_All_consenDri.txt',out)	  

Clean('FullData*.meg')	
Clean('FullData*.txt')
Clean('FullData*.nwk') 
Clean('FullData*.fas')
Clean('FullData*.fasta') 
#Clean('FullData*summary.txt') 
#Clean('FullData_Rep*.txt')
#Clean('FullData*_0.02.fas')
#Clean('FullData*_0.02.meg')
#Clean('FullData*_0.02.txt')
#os.remove('Tree_bootstrap.nwk')
#os.remove('FullData.txt')
#Clean('FullDatasnv_CloneFinder_All.*')
Clean(CFin[:-4]+'_All_consen*_partitions.txt')
os.remove(CFin[:-4]+'_All.fasta')
os.remove(CFin[:-4]+'_All.txt')
os.remove(CFin[:-4]+'_All_consen_CloneAnno.txt')
if os.path.exists(TuAnno)!=True: 
   print ('PathFinder was not performed because '+TuAnno+' not found')