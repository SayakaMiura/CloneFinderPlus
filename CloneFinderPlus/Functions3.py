import os
import shutil
import glob
from shutil import copy
import numpy
#from ete2 import Tree
from Bio import Phylo
import pandas as pd
import numpy as np
from itertools import combinations

from Bio.Phylo.Consensus import _BitString

def Clean(Target):
	filelist = glob.glob(Target)
	for f in filelist:
		os.remove(f)
		
def ListMutPos(ASeq,DSeq):
            Len=len(ASeq)		
            Dif=[]
            c=0
            while c<Len:
                if ASeq[c]=='A' and DSeq[c]=='T': Dif.append(c)
                c+=1
            return Dif
def ReadFasSeq(Fas):
 ID2Seq={}
 Fas=open(Fas,'r').readlines()
 IDls=[]
 for i in Fas:
  # print (i)
   if i[0]=='>':
       #if Seq!='start': ID2Seq[ID]=Seq
       ID=i.strip()
       IDls.append(ID)	   
       ID2Seq[ID]=''	   
      # Seq=''

   else: 
     Seq=i.strip()
     ID2Seq[ID]+=Seq
 return IDls,ID2Seq	 
	

def ReadMegSeq(Meg): 
  Meg=open(Meg,'r').readlines()
  Read='s'
  out2=''
  NameOrder=[]
  Name2Seq={}
  for i in Meg:
    if i[0]=='#' and i.strip()!='#MEGA' and i.strip()!='#mega' :
        Read='n'
        Name=i.strip()
        NameOrder.append(Name)
        Name2Seq[Name]=''
    elif Read=='n': Name2Seq[Name]+=i.strip()
    elif Read=='s': out2+=i
  return NameOrder, Name2Seq, out2
  


def UpMeg(Name2Seq,NameLs,AA,Out):
    out='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n\n'
    for Name in NameLs:
        if Name[0]!='#': Name='#'+Name
        out+=Name+'\n'+Name2Seq[Name]+'\n'
    GetOut(Out,out)	
def all_parents(tree):
    parents = {}
    for clade in tree.find_clades(order='level'):
        for child in clade:
            parents[child] = clade
    return parents	

		
def Meg2MP(Meg,RootTaxa,Rpath):
    MegID=Meg[:-4]
    CloLs,Clo2Seq,out=ReadMegSeq(Meg)
    if CloLs.count('#hg19')==0: Clo2Seq['#hg19']='A'*(len(Clo2Seq[CloLs[0]]))
    for Clo in Clo2Seq:
        out+=Clo+'\n'+Clo2Seq[Clo]+'\n'
    GetOut(MegID+'_WithOut.meg',out)		
    os.system('megacc -a infer_MP_nucleotide.mao -d '+MegID+'_WithOut.meg'+' -o '+MegID+'.nwk') 
  
    if os.path.exists(MegID+'.nwk')==True:      
        os.system('python MakeTopologeTree.py ' +MegID+'.nwk '+Rpath)
    if os.path.exists(MegID+'multi.nwk')==True:  			 
       Nwkstr=open(MegID+'multi.nwk','r').readlines()[0]
              
       Tree=RootTree(Nwkstr, RootTaxa[1:])
       os.remove(MegID+'.nwk')	 
       os.remove(MegID+'multi.nwk')	
       os.remove(MegID+'_summary.txt')	
       AncFile=glob.glob(MegID+'_ancestral_states_*.txt')
       for File in AncFile:
           os.remove(File) 		   
    else: Tree='NA\n'	   
    os.remove(MegID+'_WithOut.meg')	

	   
    return Tree 		
def do_MP(InMeg,Rpath):
    Tree=''
    if os.path.exists(InMeg[:-4]+'.nwk')!=True:
      CloLs, Clo2Seq, AA=ReadMegSeq(InMeg)
      Clo2Seq['#Normal']='A'*len(Clo2Seq[CloLs[0]]) 
      UpMeg(Clo2Seq,Clo2Seq.keys(),'AA',InMeg[:-4]+'_withNarmal.meg')
      os.system('megacc -a infer_MP_nucleotide.mao -d '+InMeg[:-4]+'_withNarmal.meg'+' -o '+InMeg[:-4]+'.nwk')

    if os.path.exists(InMeg[:-4]+'.nwk')==True and os.path.exists(InMeg[:-4]+'multi.nwk')!=True:# or InfRenameMeg.find('CloneFinder')!=-1:		
             os.system('python MakeTopologeTree.py ' +InMeg[:-4]+'.nwk'+' '+Rpath)
    if os.path.exists(InMeg[:-4]+'multi.nwk')==True:			 
             Tree_str=open(InMeg[:-4]+'multi.nwk','r').readlines()[0]
             Tree=RootTree(Tree_str, 'Normal')
    if os.path.exists(InMeg[:-4]+'.nwk')==True: os.remove(InMeg[:-4]+'.nwk')
    if os.path.exists(InMeg[:-4]+'_withNarmal.meg')==True:  os.remove(InMeg[:-4]+'_withNarmal.meg')	
    if os.path.exists(InMeg[:-4]+'multi.nwk')==True:  os.remove(InMeg[:-4]+'multi.nwk')	
    return Tree
	
def RootTree(OriNwk, RootTaxa):	
        OutF=open('test.nwk','w')
        OutF.write(OriNwk)
        OutF.close()		
        trees = list(Phylo.parse('test.nwk', 'newick'))
        for tree in trees:
           tree = tree.root_with_outgroup({'name': RootTaxa})
        Phylo.write(trees, 'test1.nwk', "newick")	

        return open('test1.nwk','r').readlines()[0]	



def GetSharePosi1(Ori2Seq,ShareNuc):
  for i in Ori2Seq:
    Name=i
  SharePosi=[]
  Len=len(Ori2Seq[Name])
  c=0
  while c<Len:
    AllMut='y'
    for Ori in Ori2Seq:
      if Ori!='#hg19' and Ori!='#Normal':
       Nuc=Ori2Seq[Ori][c]
       if Nuc!=ShareNuc: AllMut='n'
    if AllMut=='y': 
        SharePosi.append(c)
    c+=1
  return SharePosi

def CountDifNum(Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif


	
def GetHead(Head):
    Head=Head.strip().split('\t')
    Len=len(Head)
    c=0
    Name2Col={}
    NameOrder=[]	
    while c<Len:
        Name2Col[Head[c]]=c
        NameOrder.append(Head[c])		
        c+=1
    return NameOrder,Name2Col 	
def GetHead_csv(Head):
    Head=Head.strip().split(',')
    Len=len(Head)
    c=0
    Name2Col={}
    NameOrder=[]	
    while c<Len:
        Name2Col[Head[c]]=c
        NameOrder.append(Head[c])		
        c+=1
    return NameOrder,Name2Col 	

def configure(input_csv):
  df = pd.read_csv(input_csv, delim_whitespace=True)
  if 'mutation_id' not in df.columns:
    df.insert(0, "mutation_id", np.nan)
    for x in range(len(df)):
      df.loc[x,'mutation_id'] = "Mut_{}".format(str(x))
    df.to_csv(input_csv, sep='\t', index = False)
    
def configureAnno(annoFile):
  df2 = pd.read_csv(annoFile, delim_whitespace=True)
  df2.columns = ['Sample', "Site"]
  for x in range(len(df2)):
    tu1 = df2.loc[x,'Sample']
    if(tu1[:2] != "T-"):
      df2.loc[x,'Sample'] = "T-{}".format(tu1)
  df2.to_csv(annoFile, sep='\t', index = False)    
 
def ListColStr(File):
  File=open(File,'r').readlines()
  NameOrder,Name2Col=GetHead(File[0])
  File=File[1:]
  Tu2Freq={}
  for Tu in NameOrder:
    Tu2Freq[Tu]=[]
  for i in File:
    i=i.strip().split('\t')
    for Tu in Name2Col:
        Tu2Freq[Tu].append(i[Name2Col[Tu]])
  return Tu2Freq	  
  
def ListColStr_csv(File):
  File=open(File,'r').readlines()
  NameOrder,Name2Col=GetHead_csv(File[0])
  File=File[1:]
  Tu2Freq={}
  for Tu in NameOrder:
    Tu2Freq[Tu]=[]
  for i in File:
    i=i.strip().split(',')
    for Tu in Name2Col:
        Tu2Freq[Tu].append(i[Name2Col[Tu]])
  return Tu2Freq	 					  
			
def PruneTree(KeepCloSet, Nwk_str):
    from ete2 import Tree
    t=Tree(Nwk_str)
    t.prune(KeepCloSet)
    t1=t.write(format=1)	

    return t1
	  

def GetOut(File,In):
    OutF=open(File,'w')
    OutF.write(In)
    OutF.close()	
	 