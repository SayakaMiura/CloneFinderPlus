import sys
import Functions3
import os

#C:\Users\kumarlab\Desktop\PathFinder_2.1\pathfinder-dev\pathfinder-dev>python38 pathfinder.py C:\Users\kumarlab\Desktop\CloneFinderPlus\Empirical\Data\Townsend\ATP402Papersnv_CloneFinder.fas C:\Users\kumarlab\Desktop\CloneFinderPlus\Empirical\Data\Townsend\ATP402Papersnv_CloneFinder.txt --primary PrimT  -o C:\Users\kumarlab\Desktop\CloneFinderPlus\Empirical\Data\Townsend --log_ancestral_inference
DriverLabIn=sys.argv[1]#'4,Yes' #4,0.05,Cut #Col,Val
DriverLabLs=DriverLabIn.split(',')
DriCol=int(DriverLabLs[0])
OutExt=sys.argv[6]
Class='n'
#print (DriverLabLs)
#open('a','r').readlines()
if len(DriverLabLs)==3:
    DriverLab=float(DriverLabLs[1])
    Cut='y'	
elif len(DriverLabLs)==2:	
    DriverLab=DriverLabLs[1]
    Cut='n'
    if DriverLab.find(':')!=-1:
        Class='y'
        CatLs=DriverLab.replace('_',' ').split(':')		
        Dri=[CatLs[0]]
        Pass=CatLs[1:]
        print (Dri,Pass)
      #  open('a','r').readlines()		
else: open('a','r').readlines()	
#ID='ATP400Paper'
DriAnno=sys.argv[2]#'C:\\Users\\kumarlab\\Desktop\\CloneFinderPlus\\Empirical\\IntogenTownsend\\IntoGen\\'+ID+'_IntoGen.txt'#Chr Pos Wild Mut Driver prediction
PFout=sys.argv[3]#'C:\\Users\\kumarlab\\Desktop\\CloneFinderPlus\\Empirical\\Data\\Townsend\\'+ID+'snv_CloneFinder_Mig.txt'
Out=PFout[:-4]+'_'+OutExt+'.txt'
print (Out)
#open('a','r').readlines()
CFin=sys.argv[4]#'C:\\Users\\kumarlab\\Desktop\\CloneFinderPlus\\Empirical\\Data\\Townsend\\'+ID+'.txt'#Chr Pos Wild Mut
AncSeqTa=sys.argv[5]#'C:\\Users\\kumarlab\\Desktop\\CloneFinderPlus\\Empirical\\Data\\Townsend\\ancestral_inference_logging\\nuc_anc_seqs_564551_eps.csv'
#NodeMap=AncSeqTa[:-8]+'_nodeMap.txt' Need this file

NodeID=AncSeqTa[:-4]+'_nodeMap_nodeNames.txt'
python='python'
MEGnode2PFnode={}
PFnode2MEGnode={}
NodeID=open(NodeID,'r').readlines()
for i in NodeID:
    i=i.split('\t')
    MEGnode2PFnode[i[1].strip()]=i[0].strip()
    PFnode2MEGnode[i[0].strip()]=i[1].strip()	
print ('get Driver')
DriAnno=open(DriAnno,'r').readlines()
Head=DriAnno[0].split('\t')
GeneInfo='n'
if Head[5]=='AA' and Head[6]=='Gene': GeneInfo='y'
DriAnno=DriAnno[1:]
Pos2Dri={}
Pos2Gene={}
print (DriverLab)
for i in DriAnno:
    i=i.split('\t')
    Pos='\t'.join(i[:4])
    D=i[DriCol].strip()	
    Pos2Dri[Pos]=D

    if GeneInfo=='y': Pos2Gene[Pos]=i[6]+':'+i[5]
        	
print (len(Pos2Dri))	
print ('get mutation order')
OutCF=CFin[:-4]+'_'+OutExt+'.txt'
CFin=open(CFin,'r').readlines()
PosOrder=[]
out=CFin[0].strip()+'\tMutType'
if GeneInfo=='y': out+='\tGene'
out+='\n'  

CFin=CFin[1:]
for i in CFin:
    out+=i.strip()+'\t'
    i=i.split('\t')
    Pos='\t'.join(i[:2])+'\t'+i[5]+'\t'+i[6]	
    PosOrder.append(Pos)
    if Cut=='n':out+=Pos2Dri[Pos].split(' ')[0]
    elif Cut=='y':
        D=Pos2Dri[Pos]	
        if float(D)<DriverLab: out+='driver'
        else: out+='passenger'		
    else: open('a','r').readlines()		
    	
    if GeneInfo=='y': out+='\t'+Pos2Gene[Pos]
    out+='\n'
Functions3.GetOut(OutCF,out)	
print (len(PosOrder))	
	
print ('make ancestral sequences')
os.system(python+' MakeAncSeq.py '+AncSeqTa)#WithAnc.fasta

print ('map driver count')
NodeLs,Node2Seq=Functions3.ReadFasSeq(AncSeqTa[:-4]+'WithAnc.fasta')
PFout=open(PFout,'r').readlines()
out=PFout[0].strip()+'\tDriver count\tPassenger count\tTotalForward\tDriverlist\n'
PFout=PFout[1:]
for i in PFout:
    out+=i.strip()+'\t'
    Path=i.split('\t')[3].split('->')
    AncSeq=Node2Seq['>'+PFnode2MEGnode[Path[0].replace('(','')]]
    DecSeq=Node2Seq['>'+PFnode2MEGnode[Path[1].split(')')[0]]]
    
    MutPos=Functions3.ListMutPos(AncSeq,DecSeq)
    DifC=Functions3.CountDifNum(AncSeq,DecSeq)	
    DriverLs=[]
    PassLs=[]	
    for M in MutPos:
        Pos=PosOrder[M]
        D=Pos2Dri[Pos]
        if Class=='y':
            D=D.split(' (')[0]
            if Dri.count(D)!=0: DriverLs.append(Pos.replace('\t','_'))
            elif Pass.count(D)!=0: PassLs.append(Pos.replace('\t','_'))
            elif D.strip()!='': 
               print (D)
               open('a','r').readlines()	
           # else:	
            #   print (D)
             #  open('a','r').readlines()				
        elif Cut=='n':		
            if D==DriverLab: DriverLs.append(Pos.replace('\t','_'))
        elif Cut=='y':
            D=D.strip()
            if D!='NA' and D!='' and D!='9999':
               D=float(D)
               if D>1: pass			   
               elif D<DriverLab: 	DriverLs.append(Pos.replace('\t','_'))
               else: PassLs.append(Pos.replace('\t','_'))			   
        else: open('a','r').readlines()		
    NonDC=len(MutPos)-len(DriverLs)
    out+=str(len(DriverLs))+'\t'+str(len(PassLs))+'\t'+str(len(MutPos))+'\t'+';'.join(DriverLs)+'\n'
Functions3.GetOut(Out,out)	
	