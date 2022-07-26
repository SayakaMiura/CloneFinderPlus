#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 19:02:05 2022

@author: jaredhuzar
"""

import sys
import Functions3
import numpy as np
import pandas as pd
import csv
import os
from os.path import exists
from FastClone_GuanLab.fastclone import __main__
#import subprocess
from alignments.MegaAlignment import MegaAlignment
import shutil
from PathFinder.output import CloneFrequencyAnalizer
from PathFinder.output.CloneFrequencyAnalizer import CloneFrequencyAnalizer
import os

Cut = 0.05

Dir = os.getcwd()
Dir5 = os.path.join(Dir,'Data')
CFin = sys.argv[1]
tuAnnoPath = sys.argv[2]

os.rename(os.path.join(Dir,CFin), os.path.join(Dir5,CFin))
if os.path.exists(tuAnnoPath)==True:
     os.rename(os.path.join(Dir,tuAnnoPath), os.path.join(Dir5,tuAnnoPath))
     Functions3.configureAnno(os.path.join(Dir5,tuAnnoPath))
    	  
Functions3.configure(os.path.join(Dir5,CFin))

Col2In=Functions3.ListColStr(os.path.join(Dir5,CFin))
MutIDor=Col2In['mutation_id']
TuLs=[]
TCut=0.99
ACut=0.01#1-TCut
ColLs=list(Col2In.keys())
for i in ColLs:
    if i.find(':ref')!=-1: TuLs.append(i.split(':ref')[0])
print (TuLs)

def CountClo(CloLs,ExLs):
    C=[]
    for i in CloLs:
        if ExLs.count(i)==0: C.append(i)	
    return C
def GetBase(Val,ACut,TCut):			  		  
              if Val<ACut: 
                   Seq='A'
              elif Val>=TCut: Seq='T'
              else: Seq='?'
              return Seq			  
print ('make FastClone sequence')



Dir2= os.path.join(Dir,'Data', 'FastCloneData', '')

ResDir = os.path.join(Dir, 'FastCloneRes', '')



def cloneFinder_to_fastClone(x,ref,alt):
    data = pd.read_table(str(x))
    newData = pd.DataFrame(np.nan, index=range(0,len(data.index)), columns=['mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'minor_cn', 'major_cn', 'variant_gene', 'variant_freq', 'variant_coding', 'variant_effect'])
    newData.iloc[:,0] = data.iloc[:,0]
    newData.iloc[:,1] = data.iloc[:,ref]
    newData.iloc[:,2] = data.iloc[:,alt]
    newData.iloc[:,3] = 2
    newData.iloc[:,4] = 0
    newData.iloc[:,5] = 2
    newData.iloc[:,6] = newData.iloc[:,0]
    freq = newData.iloc[:,2]/(newData.iloc[:,1]+newData.iloc[:,2])
    newData.iloc[:,7] = freq
    newData.iloc[:,8] = newData.iloc[:,0]
    newData.iloc[:,9] = "missense"
    return newData

def seperateSamples(z):
    data = pd.read_table(z)
    SampLs=[]
    for i in data.columns:
        if i.find(':')!=-1: SampLs.append(i.split(':')[0])
    SampLs=list(set(SampLs))

    for sample in SampLs:	

        newData = cloneFinder_to_fastClone(z,list(data.columns).index(sample+':ref'),list(data.columns).index(sample+':alt'))

        newData.to_csv(os.path.join(Dir2, sample + ".tsv"), sep = "\t", index = False)
  		



seperateSamples(os.path.join(Dir5,CFin))
filesList = os.listdir(Dir2)
for x in range(0,len(filesList)):
    instance1 = __main__.Entrypoint()
    instance1.load_pyclone(filesList[x].split(".tsv")[0], Dir2+ filesList[x], None)
    instance1.solve(ResDir + filesList[x].split(".tsv")[0])



out=''
Allout='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n\n'
CloDic={}	
CloAnno='CloneID\tFastCloneID\n'
CloneID=1
DoneDic={}
for Tu in TuLs:
    print ('read tree')

    FCdir=os.path.join(ResDir,Tu)

    Tree=os.path.join(FCdir,'phylogeny_proportion.csv')
    if os.path.exists(Tree)==True:	
     Tree=open(Tree,'r').readlines()[1:]
     Dec2Anc={}
     AncLs=[]
     for i in Tree:
       i=i.split(',')
       Anc=i[1]
       Dec=i[2]
       Dec2Anc[Dec]=Anc	
       AncLs.append(Anc)
     AncLs=list(set(AncLs))	   
     TipLs=[]
     for Dec in Dec2Anc:
        if AncLs.count(Dec)==0: TipLs.append(Dec)

     print ('read SNV assignment')
    else:
     Dec2Anc={}
     AncLs=[]	
     TipLs=['0']
    Geno=os.path.join(FCdir,'scores.csv')
    Col2Geno=Functions3.ListColStr_csv(Geno)

    PosLs=Col2Geno['']
    CloLs=AncLs+TipLs	

    ExcLs=[] 
    
    for Clo in CloLs:
        Clo=Clo.strip()	
	
        ValLs=Col2Geno[Clo]	
        if Clo in Dec2Anc: AncValLs=Col2Geno[Dec2Anc[Clo]]
        else: AncValLs=[]		
        Seq=''	
        AmbP=[]
        DetP=[]		
        for Mut in MutIDor:
           if PosLs.count(Mut)==0:
               	Seq+='A'
           else:
              P=PosLs.index(Mut)
		  
              Val=float(ValLs[P])
			  
              B=GetBase(Val,ACut,TCut)
              if B=='?': AmbP.append(Mut)
              else: DetP.append(Mut)			  
			  
              if Clo in Dec2Anc: AB=GetBase(float(AncValLs[P]),ACut,TCut)	
              else: AB=''
              if B=='T': Seq+=B
              else:
                if Clo in Dec2Anc: Seq+=AB
                else: Seq+=B				


        if len(AmbP)>len(DetP):
           print (Tu,'mostly ambiguous, so use tumor seq')
           Seq=''
           AllMut=AmbP+DetP		   
           for Mut in MutIDor:
                if AllMut.count(Mut)!=0: Seq+='T'
                else: Seq+='A'				
        out+='>'+Tu+'_'+Clo+'\n'+Seq+'\n'
        CloDic[Tu+'_'+Clo]=Seq		
        if Seq.find('T')==-1: ExcLs.append(Clo)
    print (Tu,'exc',ExcLs)

    Tc=CountClo(TipLs,ExcLs)
    Ac=CountClo(AncLs,ExcLs)	
    print ('count tip',Tc,Ac)
    if len(Tc)==1 and len(Ac)==0:
        Seq=CloDic[Tu+'_'+Tc[0]]	
        DoneDic[Seq]=DoneDic.get(Seq,0)+1
        SeqC=DoneDic[Seq]	
        if SeqC==1:		
           Allout+= '#Clone'+str(CloneID)+'\n'+Seq+'\n'
           CloAnno+='Clone'+str(CloneID)+'\t'+Tu+'_'+Tc[0]+'\n'
           CloneID+=1		
    elif len(Tc)==1 and len(Ac)>=1:
        Anc=Dec2Anc[Tc[0]]
        Seq=CloDic[Tu+'_'+Anc]	
        DoneDic[Seq]=DoneDic.get(Seq,0)+1
        SeqC=DoneDic[Seq]	
        if SeqC==1:		
            Allout+= '#Clone'+str(CloneID)+'\n'+Seq+'\n'
            CloAnno+='Clone'+str(CloneID)+'\t'+Tu+'_'+Anc+'\n'
            CloneID+=1
    elif len(Tc)>1:
        for T in Tc:
          Seq=CloDic[Tu+'_'+T]	
          DoneDic[Seq]=DoneDic.get(Seq,0)+1
          SeqC=DoneDic[Seq]	
          if SeqC==1:			
            Allout+= '#Clone'+str(CloneID)+'\n'+Seq+'\n'
            CloAnno+='Clone'+str(CloneID)+'\t'+Tu+'_'+T+'\n'
            CloneID+=1
    else:
        print ('all clones are removed, so tumor seq was added',Ac,Tc)
        Seq=''		
        for Mut in MutIDor:
           if PosLs.count(Mut)==0:
               	Seq+='A'
           else: Seq+='T'
        DoneDic[Seq]=DoneDic.get(Seq,0)+1
        SeqC=DoneDic[Seq]		
        if SeqC==1:					   
            Allout+= '#Clone'+str(CloneID)+'\n'+Seq+'\n'
            CloAnno+='Clone'+str(CloneID)+'\t'+Tu+'_tumorSeq\n'
            CloneID+=1		
			

Functions3.GetOut(os.path.join(Dir5,CFin[:-4])+'_FastClone.fasta',out) 		

Allout+='#hg19\n'+('A'*len(MutIDor))+'\n'

Functions3.GetOut(os.path.join(Dir5,CFin[:-4])+'_IniClone.meg',Allout)	

Functions3.GetOut(os.path.join(Dir5,CFin[:-4])+'_IniCloneFastCloneAnno.txt',CloAnno)

FCLs,FCSeq=Functions3.ReadFasSeq(os.path.join(Dir5,CFin[:-4])+'_FastClone.fasta') 


os.system('python clonefinderPlus1.py snv {} {}'.format(os.path.join(Dir5,CFin), os.path.join(Dir5,CFin[:len(CFin) - 4] + "_IniClone.meg")))

makePathFinder = CFin[:len(CFin) - 4] + "snv_CloneFinderPlus" + ".meg"

def editPathFinderIn(fasFile,txtFile):
    file1 = open(fasFile, 'r')
    Lines = file1.readlines()
    file2 = open(txtFile, 'r')
    Lines2 = file2.readlines()
    tumorsIn=Lines2[0].replace('\n', '').split('\t')[1::]
    tumorsIn.append('Normal')
    cloneTot = []
    for line in Lines:
        if 'Clone' in line:
            cloneTot.append(line.replace('>', '').replace('\n', ''))
    if(len(list(set(cloneTot).difference(set(tumorsIn)))) > 0):
        for x in range(len(Lines)):
            if list(set(cloneTot).difference(set(tumorsIn)))[0] in Lines[x]:
                Lines.remove(Lines[x])
                Lines.remove(Lines[x])
                break
        file1.close()
        os.remove(fasFile)
        file3 = open(fasFile, 'w')
        for i in range(len(Lines)):
            file3.write(Lines[i])
        file3.close()

if os.path.exists(os.path.join(os.path.dirname(os.getcwd()),tuAnnoPath))==True:

    os.system('python make_PathFinder_input2.py {} {}'.format(os.path.join(os.path.dirname(os.getcwd()),makePathFinder),os.path.join(os.path.dirname(os.getcwd()),tuAnnoPath)))
    
    pathFinderIn = makePathFinder[:len(makePathFinder) - 4] + "_PathFinder_" + str(Cut) + ".fas"
    
    editPathFinderIn(os.path.join(Dir5,pathFinderIn),os.path.join(Dir5,pathFinderIn[:len(pathFinderIn) - 4]+".txt"))

    PFcom='python pathfinder.py '+os.path.join(Dir5,pathFinderIn)+' '+ os.path.join(Dir5,pathFinderIn[:len(pathFinderIn) - 4])+".txt"+' --primary P -o '+os.getcwd()+'\\PathFinderRes --relax_threshold  --max_graphs_per_tree 100  --log_ancestral_inference'	
    print (PFcom)
    os.chdir('PathFinder')	
    os.system(PFcom)
    os.chdir(os.getcwd())    	
 





