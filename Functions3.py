import os
import shutil
import glob
from shutil import copy
import numpy
from Bio import Phylo


def Clean(Target):
	filelist = glob.glob(Target)
	for f in filelist:
		os.remove(f)
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

def GV2MutOr(Clones,Len): #'_map0.gv',Len: MutNum
  #  CellOrderF=Clones[:-13]+'_CellOrder.txt'
   # CellOrder=open(CellOrderF,'r').readlines()
    #Num=len(CellOrder)
    #print (Num)  
 #   OutSeq=Clones[:-7]+'.meg'#CellSeqFol+ID+'SCITE.meg'   
    Clones=open(Clones,'r').readlines()
 #   Cell2Clo={}
 #   Clo2Cell={}
 #   CloLs=[]
 #   Dec2Anc={}
   # out='Cell,Clone\n'
   # outSeq='#mega\n!Title Cell;\n!Format DataType=DNA indel=-;\n'
  #  MutPosi=[]
    MutOrLs=[]
    for i in Clones:
      if i.find('-> ssnv')!=-1:

       pass
      elif i.find('-> ')!=-1:
       MutOrLs.append(i.strip().replace(';','').replace(' ',''))
    #   i=i.strip().split('->')
    #   Dec2Anc['C'+i[1].strip().replace(';','')]='C'+i[0].strip()
     #  MutPosi+=[int(i[0].strip())-1,int(i[1].replace(';','').strip())-1]   
       
  # MutPosi=list(set(MutPosi))
 # GetOut(OutSeq,outSeq) 
    return MutOrLs 	 
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

  return NameOrder, Name2Seq

def Meg2PreAbs(Meg): #OriMeg[:-4]+'_PreAbs.txt' from 1
  CloLs,Clo2Seq=ReadMegSeq(Meg)
  out='Clone\tPos\tAssign\n'
  Clo2MutPosOri={}
  SNVc=len(Clo2Seq[CloLs[0]])
  for Clo in CloLs:
     c=0
     MutLs=[]	 
     while c<SNVc:
         Ass=Clo2Seq[Clo][c]	 
         out+=Clo.replace('#','')+'\t'+str(c+1)+'\t'+Ass+'\n'
         if Ass=='T': MutLs.append(c+1)
         c+=1
     Clo2MutPosOri[Clo]=MutLs	
  GetOut(Meg[:-4]+'_PreAbs.txt',out)
  return  Clo2MutPosOri 
def Invert(Dic):
    IDic={}
    for i in Dic:
      IDic[Dic[i]]=IDic.get(Dic[i],[])+[i]	
    return IDic	  
def GetHeadCloneFreq(Head):
        Head=Head.strip().split('\t')
        Len=len(Head)
        c=1
        Name2Col={}
        while c<Len:
            Name2Col[Head[c]]=c
            c+=1
        return Name2Col
def GetCloHitForTu(File,Cut):
        File=open(File,'r').readlines()
        Name2Col=GetHeadCloneFreq(File[0])
        File=File[1:]
        Tu2Clo={}
        T2C2F={}
        HitCloLs=[]		
        for i in File:
            i=i.strip().split('\t')	
            Tu=i[0]
            Tu2Clo[Tu]=[]
            T2C2F[Tu]={}		
            for Name in Name2Col:
                Fre0=i[Name2Col[Name]]		
                Fre=float(Fre0)
                if Fre0.find('-')==-1 and Fre>Cut:
                     T2C2F[Tu][Name]=Fre
                     Tu2Clo[Tu].append(Name)
                     HitCloLs.append(Name)					 
                else: T2C2F[Tu][Name]=0			
        HitCloLs=list(set(HitCloLs))       
        return Tu2Clo,HitCloLs,T2C2F
def ListColStr(File):
  File=open(File,'r').readlines()
  NameOrder,Name2Col=GetHead(File[0])
  File=File[1:]
  Tu2Freq={}
  for Tu in NameOrder:
    Tu2Freq[Tu]=[]
  for i in File:
    if i.find('\t')!=-1: i=i.strip().split('\t')
    else: i=i.strip().split(' ')	
    for Tu in Name2Col:
        Tu2Freq[Tu].append(i[Name2Col[Tu]])
  return Tu2Freq	
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
def GetHeadObsFreqTaHead(Head): 
  Head=Head.strip().split('\t')
  SampOrder=[]
  Name2Col={}
  c=0
  Len=len(Head)
  while c<Len:
    i=Head[c]
    if i.find(':')!=-1:
        Samp=i.split(':')[0]
        Code=Samp in SampOrder
        if Code!=True:
             SampOrder.append(Samp)
        Name2Col[i]=c
    c+=1
  return SampOrder,Name2Col		
def ObsFreq2SampSeq(Freq): 
  SNVOut=Freq.replace('.txt','_tumor.fas')
  Freq=open(Freq,'r').readlines() 
  SampOrder,Name2Col=GetHeadObsFreqTaHead(Freq[0])
  Samp2FreqIn={}
  Samp2TotRead={}
  for Samp in SampOrder:
      Samp2FreqIn[Samp]=[]
      Samp2TotRead[Samp]=[]	  
  Freq=Freq[1:]
  outS=''
  SNVnum=0
  SeqOut={}
  for i in Freq:
    i=i.strip().split('\t')
    TMP={}
    for Samp in SampOrder:
        
        if SNVnum==0: outS+=Samp+'\t'
        MutC=int(i[Name2Col[Samp+':alt']])
        RefC=int(i[Name2Col[Samp+':ref']])
        if (MutC+RefC)==0: MutFreq='NA'		
        else: MutFreq=1.0*MutC/(MutC+RefC)
        		
        Tot=MutC+RefC
        Samp2FreqIn[Samp].append(MutFreq)
        Samp2TotRead[Samp].append(Tot)
        if MutC==0: SeqOut[Samp]=SeqOut.get(Samp,'')+'A'
        else: SeqOut[Samp]=SeqOut.get(Samp,'')+'T'
     		
    SNVnum+=1
  out='>Normal\n'+('A'*SNVnum)+'\n'
  for S in SeqOut:
      out+='>'+S+'\n'+SeqOut[S]+'\n'
  GetOut(SNVOut,out)  
  outS=outS[:-1]+'\n'
  c=0
  while c<SNVnum:
    for Samp in SampOrder:
        outS+=str(Samp2FreqIn[Samp][c])+'\t'
    outS=outS[:-1]+'\n'
    c+=1
 # GetOut(SNVOut.split('\\')[-1],outS)
  return SampOrder, Samp2FreqIn	,SNVnum,Samp2TotRead  
  
def CountDifNum(Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif
def CountDifNum_excMis(Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c] and Seq0[c]!='?' and Seq1[c]!='?': Dif+=1
                c+=1
            return Dif			
def Pair_clone(MegT,MegI):
      OriLs,Ori2Seq=ReadMegSeq(MegT)
      BooLs,Boo2Seq=ReadMegSeq(MegI)
      Old2NewID={}
      for Ori in OriLs:
         if Ori!='#Normal' and Ori!='#hg19':
             OSeq=Ori2Seq[Ori]
             BestC=9999
             BestSeq=''			 
             for Boo in BooLs:
                 if Boo!='#Normal' and Boo!='#hg19':
                     BSeq=Boo2Seq[Boo]				 
                     DifC=CountDifNum(OSeq,BSeq)
                     if DifC<BestC: 
                        BestC=DifC
                        BestSeq=Boo	
             #out+=Ori+'\n'+Boo2Seq[BestSeq]+'\n'
             Old2NewID[BestSeq.replace('#','')]=Old2NewID.get(BestSeq.replace('#',''),[])+[Ori.replace('#','')]
      print ('clone annotation',Old2NewID)
      return Old2NewID #inf2true	  
def ReadBooTree(Tree):
    Tree=open(Tree,'r').readlines()[0]
    BooLen=[]
    Boo=''
    Blen=''	
    c=0
    Len=len(Tree)
    Read='s'	
    while c<Len:
     #   print (c,Read,Boo,Blen)	
	   
        if Read=='b':
           if Tree[c]==':':
              Read='l'	
              Blen=''
           else:
              Boo+=Tree[c]
        elif Read=='l':
           if Tree[c]==',' or Tree[c]==')':
             Read='s'
             BooLen.append([float(Boo),float(Blen)])
             Boo=''
             Blen=''
             if Tree[c]==')': 
                Read='b'			 
           else:
             Blen+=Tree[c]
        elif Tree[c]==')':
           Read='b'
           Boo=''				 
        c+=1			 
    if Blen!='': BooLen.append([float(Boo),float(Blen)])
    return BooLen	
def GetSimClo(Seq0,SeqDic,Exc,Max):
   IdenLs=[]
   for ID in SeqDic:
       if Exc.count(ID)==0:     
          Seq1=SeqDic[ID]
          DifC=CountDifNum_excMis(Seq0,Seq1)
          if DifC<=Max:
             IdenLs.append(ID)
   return IdenLs
def CountBaseAssign(Seq,Dic,C):
    Len=len(Seq)
    c=0	
    while c<Len:
       if c not in Dic: Dic[c]={}
       Dic[c][Seq[c]]=Dic[c].get(Seq[c],0)+C	
       c+=1
    return Dic	   
def MakeConsenGeno(Fas,MaxDiff,MinPro,TotBoo):	
    ConAssinMin=0.9
    BooClo2BooLs={}
    File=Fas[:-6]+'.txt'
    CloFreFile='FullDatasnv_CloneFinder_AllCloFre.txt'#Fas[:-6]+'CloFre.txt'	
    File=open(File,'r').readlines()[1:]
    for i in File:
       i=i.split('\t')
       BooClo=i[0]	   
       BooClo2BooLs[BooClo]=BooClo2BooLs.get(BooClo,[])+[i[1].strip()]	   
    BooCloLs0,BooClo2Seq0=ReadFasSeq(Fas)
    Len=len(BooClo2Seq0[BooCloLs0[0]])	
    C2CloLs={}
    BooClo2Seq={}	
#    print (len(BooCloLs0))	
    for B in BooCloLs0:
 #     print (B)	
      if B[1]=='G':
          if B.find('bc')==-1: C=1	  
          else: C=int(B.split('bc')[-1])
          C2CloLs[C]=C2CloLs.get(C,[])+[B]
          BooClo2Seq[B]=BooClo2Seq0[B]
    BooCloLs=list(BooClo2Seq.keys())		  
  #  Cls=list(C2CloLs.keys())
  #  Cls.sort()
   # print(Cls)
    Done=[]
    Clo2Iden={}	
    print (BooCloLs)	
    for B in BooCloLs:
     # if B=='>G39bc2': print (B) 
     # if B=='>G67bc1': print (B)	
      if Done.count(B)==0:
        Pri='n'	  
        Done.append(B)		
        SimLs=GetSimClo(BooClo2Seq[B],BooClo2Seq,Done,MaxDiff)	
        if B=='>G39bc2' or B=='>G67bc1' or SimLs.count('>G39bc2')!=0 or SimLs.count('>G67bc1')!=0: 
           print (B,SimLs)
           Pri='y'		   
        All=[]
        All+=SimLs		
        while SimLs!=[]:
            SimLsUp=[]		
            for B1 in SimLs:
                Done.append(B1)			
                SimLs0=	GetSimClo(BooClo2Seq[B1],BooClo2Seq,Done,MaxDiff)
                SimLsUp+=SimLs0	
            All+=SimLsUp	
            SimLs=list(set(SimLsUp))
        if Pri=='y': print (B,All)				
        Clo2Iden[B]=All	
    print (len(Clo2Iden),len(BooClo2Seq))
    print (Clo2Iden)	
    outFas='>Normal\n'+('A'*len(BooClo2Seq[B]))+'\n'
    out='Clone\tPosition\tAssign\tMutSupportProportion\tWild\tMut\tMiss\n'	
    out1='OriginalClone\tClone\n'	
    Old2NewCloID={}	
    NewCloLs=[]		
    for Clo in Clo2Iden:
       BC1=	BooClo2BooLs[Clo.replace('>','')]
      # if Clo.find('bc')==-1: Tot=1	
      # else: Tot=int(Clo.split('bc')[-1])
       IdenLs=Clo2Iden[Clo]
       #if Clo.find('bc')==-1:
        #  CloGeno=CountBaseAssign(BooClo2Seq[Clo],{},1)	            	   
       CloGeno=CountBaseAssign(BooClo2Seq[Clo],{},int(Clo.split('bc')[-1]))	  
       IdenLs=list(set(IdenLs))	   
       for I in IdenLs:
          BC1+=BooClo2BooLs[I.replace('>','')]	   
         # if Clo.find('bc')==-1: Tot+=1	   
         # else: Tot+=int(I.split('bc')[-1])
        #  if Clo.find('bc')==-1: CloGeno=CountBaseAssign(BooClo2Seq[I],CloGeno,1)			  
          CloGeno=CountBaseAssign(BooClo2Seq[I],CloGeno,int(I.split('bc')[-1]))
       BC1=list(set(BC1))
       Tot=len(BC1)	   
       Sup=1.0*Tot/TotBoo
      # if Sup>1: Sup=1	   
       if Sup>=MinPro:
          Cseq=''
          c=0	
          while c<Len:
             Info=CloGeno[c]
             C2NucLs=Invert(Info)
             Cls=list(C2NucLs.keys())
             Cls.sort()
          #   print (Cls)
            # Cnuc=C2NucLs[Cls[-1]]
            # if len(Cnuc)!=1: N='?'
            # else: N=Cnuc[0]
            # Cseq+=N	

             MutC=Info.get('T',0)
             WildC=Info.get('A',0)			 
             Tot=MutC+WildC+Info.get('?',0)			 
             MutPro=1.0*MutC/Tot	
             WildPro=1.0*WildC/Tot			 

             if MutPro>ConAssinMin: N='T'	
             elif WildPro>ConAssinMin: N='A'
             else: N='?'				
             out+=Clo.split('bc')[0]+'bc'+str(round(Sup*100))+'\t'+str(c+1)+'\t'+N+'\t'
             out+=str(MutPro)+'\t'+str(Info.get('A',0))+'\t'+str(Info.get('T',0))+'\t'+str(Info.get('?',0))+'\n'			 
             Cseq+=N			 
             c+=1			 
          outFas+=Clo.split('bc')[0]+'bc'+str(round(Sup*100))+'\n'+Cseq+'\n'
          out1+=Clo+'\t'+Clo.split('bc')[0]+'bc'+str(round(Sup*100))+'\n'
          Old2NewCloID[Clo]=Clo.split('bc')[0]+'bc'+str(round(Sup*100))	
          NewCloLs.append(Clo.split('bc')[0]+'bc'+str(round(Sup*100)))			  
          for I in IdenLs:
             out1+=I+'\t'+Clo.split('bc')[0]+'bc'+str(round(Sup*100))+'\n'
             Old2NewCloID[I]=Clo.split('bc')[0]+'bc'+str(round(Sup*100))				 
       else:
          out1+=Clo+'\tNone\n'
          Old2NewCloID[Clo]='None'		  
          for I in IdenLs:
              out1+=I+'\tNone\n'
              Old2NewCloID[I]='None'				  
    GetOut(Fas[:-6]+'_consen.fasta',outFas)		  
    GetOut(Fas[:-6]+'_consen.txt',out)			
    GetOut(Fas[:-6]+'_consen_CloneAnno.txt',out1)
    print (Old2NewCloID)
#    CloFreFile=open(CloFreFile,'r').readlines()[1:]	
#    Rep2Tu2NewCloFre={}

#    out='BooRep\tTumor\tOriginalCloneID\tBooCloneID\tCloneFreq\tNewBooCloneID\n'
#    for i0 in CloFreFile:
#        i=i0.split('\t')
#        Old=i[3]
#        New=Old2NewCloID['>'+Old]
#        Rep=i[0]
#        Tu=i[1]		
#        if Rep not in Rep2Tu2NewCloFre: Rep2Tu2NewCloFre[Rep]={}
#        if Tu not in Rep2Tu2NewCloFre[Rep]:Rep2Tu2NewCloFre[Rep][Tu]={}
#        Rep2Tu2NewCloFre[Rep][Tu][New]=Rep2Tu2NewCloFre[Rep][Tu].get(New,[])+[float(i[4])]
#        out+=i0[:-1]+'\t'+New.replace('>','')+'\n'
		
#    GetOut(Fas[:-6]+'CloFre.txt',out)	   		
#    print (Rep2Tu2NewCloFre)
#    RepNum=len(Rep2Tu2NewCloFre)
#    Tu2CloFre={}
#    for Rep in Rep2Tu2NewCloFre:
#        TuNewCloFre=Rep2Tu2NewCloFre[Rep]
#        for Tu in TuNewCloFre:
#            NewCloFre=TuNewCloFre[Tu]
#            if Tu not in Tu2CloFre: Tu2CloFre[Tu]={}
#            for NewClo in NewCloFre:
#               FreLs=NewCloFre[NewClo]
#               AveFre=1.0*sum(FreLs)/len(FreLs)
#               Tu2CloFre[Tu][NewClo]=Tu2CloFre[Tu].get(NewClo,[])+[AveFre]
#    out='Tumor\tClone\tAveFre\tSD\n'
#    out1='Tumor\t'+'\t'.join(NewCloLs)+'\n'
#    out1=out1.replace('>','')	
#    for Tu in Tu2CloFre:
#        CloFre=Tu2CloFre[Tu]
#        out1+=Tu		
#        for Clo in NewCloLs:
#           Fre=CloFre.get(Clo,[])
#           if Fre==[]:
#              Ave=0
#              SD=0	
#           else:
#              Zeros=RepNum-len(Fre)	
#              Fre=Fre+([0] * Zeros)			  
#              SD=numpy.std(Fre,ddof=1)	
#              Ave=1.0*sum(Fre)/len(Fre)			  
#           out+=Tu+'\t'+Clo.replace('>','')+'\t'+str(Ave)+'\t'+str(SD)+'\n'
#           out1+='\t'+str(Ave)	
#        out1+='\n'
#    GetOut(Fas[:-6]+'_consenCloFreSD.txt',out)
#    GetOut(Fas[:-6]+'_consenCloFre.txt',out1)
def ReadGV(Clones):
    Clones=open(Clones,'r').readlines()
    TipLs=[]
    Dec2Anc={}
    Anc2Dec={}
    for i in Clones:
      if i.find('-> ssnv')!=-1:

       i=i.strip().split('->')
       TipLs.append(i[0].strip())
      elif i.find('-> ')!=-1:
       i=i.strip().split('->')
       Dec=i[1].strip().replace(';','')
       Anc=i[0].strip()
       Dec2Anc[Dec]=Anc
       Anc2Dec[Anc]=Anc2Dec.get(Anc,[])+[Dec]
    for D in Dec2Anc:
       if D not in Anc2Dec: TipLs.append(D)
    TipLs=list(set(TipLs))
    return Dec2Anc,Anc2Dec,TipLs
def Ls2Tup(list):
    return tuple(i for i in list)	
def SumSCITEout(GV): #gv
    Dec2Anc,Anc2Dec,TipDec=ReadGV(GV)
    print (Dec2Anc,Anc2Dec,TipDec)
    NodeLs=[]
    for Anc in Anc2Dec:
        if len(Anc2Dec[Anc])>1: NodeLs.append(Anc)
    NodeLs+=TipDec
    NodeLs=list(set(NodeLs))
    print (NodeLs)
   # open('a','r').readlines()
    Dec2AncE={}
    Max0=9999999999999999999999999999999999999999999999999999999999999
    Node2MutLs={}
    for Node in NodeLs:
        if Node in Dec2Anc:
            MutLs=[Node]
            Anc=Dec2Anc[Node]
            while Anc in Dec2Anc:

                # print (Node,Anc)                
                 if NodeLs.count(Anc)!=0:
                 #    print (Anc)
                     Dec2AncE[Node]=Anc 
                     Node2MutLs[Node]=MutLs
                     Anc=Max0
                 else:
                    MutLs.append(Anc) 
                    Anc=Dec2Anc[Anc]
            if Anc!=Max0 :
                Dec2AncE[Node]=Anc
                Node2MutLs[Node]=MutLs        
    print(Dec2AncE,Node2MutLs)
    return (Dec2AncE,Node2MutLs)

               
        


def SumBooMutTree(oGV,cGV):
    cGV=open(cGV,'r').readlines()
    edge2sup={}
    for i in cGV:
       if i.find('->')!=-1: #75->41[label="2P20.0"];
           i=i.split('[')
           edge2sup[i[0].strip()]=float(i[1].split('\"')[1].split('P')[1])
    #print (edge2sup)
    Dec2AncE,Node2MutLs=SumSCITEout(oGV) 
    out='digraph {\n'
    for DecE in Dec2AncE:
      AncE=Dec2AncE[DecE]
      DecMutLs=Node2MutLs[DecE]
      AncMutLs=Node2MutLs.get(AncE,[AncE])
      In=[]
      for Dm in DecMutLs:
          Tot=0
          for Am in AncMutLs:
              E=Am+'->'+Dm
              S=edge2sup.get(E,0)
              Tot+=S
          for Dm1 in DecMutLs:
              if Dm!=Dm1:
                 E=Dm1+'->'+Dm
                 S=edge2sup.get(E,0)
                 Tot+=S                 			  
          In.append(Dm+': '+str(Tot))
      out+='m'+'m'.join(AncMutLs)+'->'+'m'+'m'.join(DecMutLs)+'[label=\"'+'\n'.join(In)+'\"];\n'
    out+='}\n'
    print (out)
    GetOut(oGV[:-3]+'_Boo.gv',out)	
def MakeConsenGeno_Li(Fas,MaxDiff,MinPro,TotBoo):	
    BooClo2BooLs={}
    File=Fas[:-6]+'.txt'
    CloFreFile='FullDataLichee.trees_All.txt'#Fas[:-6]+'CloFre.txt'	
    File=open(File,'r').readlines()[1:]
    for i in File:
       i=i.split('\t')
       BooClo=i[0]	   
       BooClo2BooLs[BooClo]=BooClo2BooLs.get(BooClo,[])+[i[1].strip()]	   
    BooCloLs0,BooClo2Seq0=ReadFasSeq(Fas)
    Len=len(BooClo2Seq0[BooCloLs0[0]])	
    C2CloLs={}
    BooClo2Seq={}	
#    print (len(BooCloLs0))	
    for B in BooCloLs0:
 #     print (B)	
      if B[1]=='G':
          if B.find('bc')==-1: C=1	  
          else: C=int(B.split('bc')[-1])
          C2CloLs[C]=C2CloLs.get(C,[])+[B]
          BooClo2Seq[B]=BooClo2Seq0[B]
    BooCloLs=list(BooClo2Seq.keys())		  
  #  Cls=list(C2CloLs.keys())
  #  Cls.sort()
   # print(Cls)
    Done=[]
    Clo2Iden={}	
    for B in BooCloLs:
      if Done.count(B)==0:	
        Done.append(B)		
        SimLs=GetSimClo(BooClo2Seq[B],BooClo2Seq,Done,MaxDiff)	
        Done+=SimLs		
        All=[]
        All+=SimLs		
        print ('\n',B,All)		
        while SimLs!=[]:
            SimLsUp=[]		
            for B1 in SimLs:
                Done.append(B1)			
                SimLs0=	GetSimClo(BooClo2Seq[B1],BooClo2Seq,Done,MaxDiff)
                print (B,B1,Done,SimLs0)				
                Done+=SimLs0				
                SimLsUp+=SimLs0	
            All+=SimLsUp	
            SimLs=list(set(SimLsUp))	
        Clo2Iden[B]=All	
    print (Clo2Iden)
  #  open('a','r').readlines()	
    print (len(Clo2Iden),len(BooClo2Seq))
    outFas='>Normal\n'+('A'*len(BooClo2Seq[B]))+'\n'
    out='Clone\tPosition\tAssign\tMutSupportProportion\tWild\tMut\tMiss\n'	
    out1='OriginalClone\tClone\n'	
    Old2NewCloID={}	
    NewCloLs=[]		
    for Clo in Clo2Iden:
       BC1=	BooClo2BooLs[Clo.replace('>','')]
      # if Clo.find('bc')==-1: Tot=1	
      # else: Tot=int(Clo.split('bc')[-1])
       IdenLs=Clo2Iden[Clo]
       #if Clo.find('bc')==-1:
        #  CloGeno=CountBaseAssign(BooClo2Seq[Clo],{},1)	            	   
       CloGeno=CountBaseAssign(BooClo2Seq[Clo],{},int(Clo.split('bc')[-1]))	  
       IdenLs=list(set(IdenLs))	   
       for I in IdenLs:
          BC1+=BooClo2BooLs[I.replace('>','')]	   
         # if Clo.find('bc')==-1: Tot+=1	   
         # else: Tot+=int(I.split('bc')[-1])
        #  if Clo.find('bc')==-1: CloGeno=CountBaseAssign(BooClo2Seq[I],CloGeno,1)			  
          CloGeno=CountBaseAssign(BooClo2Seq[I],CloGeno,int(I.split('bc')[-1]))
       BC1=list(set(BC1))
       Tot=len(BC1)	   
       Sup=1.0*Tot/TotBoo
      # if Sup>1: Sup=1	   
       if Sup>=MinPro:
          Cseq=''
          c=0	
          while c<Len:
             Info=CloGeno[c]
             C2NucLs=Invert(Info)
             Cls=list(C2NucLs.keys())
             Cls.sort()
          #   print (Cls)
             Cnuc=C2NucLs[Cls[-1]]
             if len(Cnuc)!=1: N='?'
             else: N=Cnuc[0]
             Cseq+=N	
             out+=Clo.split('bc')[0]+'bc'+str(round(Sup*100))+'\t'+str(c+1)+'\t'+N+'\t'
             MutC=Info.get('T',0)
             if MutC==0: MutPro=0
             else:			 
               MW=MutC+Info.get('A',0)			 
               MutPro=1.0*MutC/MW	
             out+=str(MutPro)+'\t'+str(Info.get('A',0))+'\t'+str(Info.get('T',0))+'\t'+str(Info.get('?',0))+'\n'	
             c+=1			 
          outFas+=Clo.split('bc')[0]+'bc'+str(round(Sup*100))+'\n'+BooClo2Seq[Clo]+'\n'
          out1+=Clo+'\t'+Clo.split('bc')[0]+'bc'+str(round(Sup*100))+'\n'
          Old2NewCloID[Clo]=Clo.split('bc')[0]+'bc'+str(round(Sup*100))	
          NewCloLs.append(Clo.split('bc')[0]+'bc'+str(round(Sup*100)))			  
          for I in IdenLs:
             out1+=I+'\t'+Clo.split('bc')[0]+'bc'+str(round(Sup*100))+'\n'
             Old2NewCloID[I]=Clo.split('bc')[0]+'bc'+str(round(Sup*100))				 
       else:
          out1+=Clo+'\tNone\n'
          Old2NewCloID[Clo]='None'		  
          for I in IdenLs:
              out1+=I+'\tNone\n'
              Old2NewCloID[I]='None'				  
    GetOut(Fas[:-6]+'_consen.fasta',outFas)		  
    GetOut(Fas[:-6]+'_consen.txt',out)			
    GetOut(Fas[:-6]+'_consen_CloneAnno.txt',out1)
	
	
def MakeConsenGeno_sc(Fas,MaxDiff,MinPro,TotBoo):	
    BooCloLs0,BooClo2Seq0=ReadFasSeq(Fas)
    Len=len(BooClo2Seq0[BooCloLs0[0]])	
    C2CloLs={}
    BooClo2Seq={}	
#    print (len(BooCloLs0))	
    for B in BooCloLs0:
 #     print (B)	
      if B[1]=='G' or B.find('Cell')!=-1:
          if B.find('bc')==-1: C=1	  
          else: C=int(B.split('bc')[-1])
          C2CloLs[C]=C2CloLs.get(C,[])+[B]
          BooClo2Seq[B]=BooClo2Seq0[B]
    BooCloLs=list(BooClo2Seq.keys())		  
    Cls=list(C2CloLs.keys())
    Cls.sort()
   # print(Cls)
    Done=[]
    Clo2Iden={}	
    for B in BooCloLs:
      if Done.count(B)==0:	
	
        SimLs=GetSimClo(BooClo2Seq[B],BooClo2Seq,Done,MaxDiff)	
        Done.append(B)
        Done+=SimLs		
        All=[]
        All+=SimLs		
      #  print ('\n',B,SimLs)		
        while SimLs!=[]:
            SimLsUp=[]		
            for B1 in SimLs:
		
                SimLs0=	GetSimClo(BooClo2Seq[B1],BooClo2Seq,Done,MaxDiff)
          #      print(B,B1,SimLs0)				
                Done+=SimLs0					
                SimLsUp+=SimLs0	
            All+=SimLsUp	
            SimLs=list(set(SimLsUp))	
        Clo2Iden[B]=All	
 #   print (Clo2Iden)
 #   open('a','r').readlines()	
    print (len(Clo2Iden),len(BooClo2Seq))
    outFas='>Normal\n'+('A'*len(BooClo2Seq[B]))+'\n'
    out='Clone\tPosition\tAssign\tMutSupportProportion\tWild\tMut\tMiss\n'	
    out1='Cell\tClone\n'	
    for Clo in Clo2Iden:
       if Clo.find('bc')==-1: Tot=1	
       else: Tot=int(Clo.split('bc')[-1])
       IdenLs=Clo2Iden[Clo]
       if Clo.find('bc')==-1:
          CloGeno=CountBaseAssign(BooClo2Seq[Clo],{},1)	            	   
       else: CloGeno=CountBaseAssign(BooClo2Seq[Clo],{},int(Clo.split('bc')[-1]))	  
       IdenLs=list(set(IdenLs))	   
       for I in IdenLs:
          if Clo.find('bc')==-1: Tot+=1	   
          else: Tot+=int(I.split('bc')[-1])
          if Clo.find('bc')==-1: CloGeno=CountBaseAssign(BooClo2Seq[I],CloGeno,1)			  
          else: CloGeno=CountBaseAssign(BooClo2Seq[I],CloGeno,int(I.split('bc')[-1]))			  
       Sup=1.0*Tot/TotBoo
       if Sup>1: Sup=1	   
       if Sup>=MinPro:
          Cseq=''
          c=0	
          while c<Len:
             Info=CloGeno[c]
             C2NucLs=Invert(Info)
             Cls=list(C2NucLs.keys())
             Cls.sort()
          #   print (Cls)
             Cnuc=C2NucLs[Cls[-1]]
             if len(Cnuc)!=1: N='?'
             else: N=Cnuc[0]
             Cseq+=N	
           #  out+=Clo.split('bc')[0]+'bc'+str(round(Sup*100))+'\t'+str(c+1)+'\t'+N+'\t'
             out+=Clo.replace('>','').replace('Cell','Clone')+'\t'+str(c+1)+'\t'+N+'\t'			 
             MutC=Info.get('T',0)
             MutPro=1.0*MutC/Tot	
             out+=str(MutPro)+'\t'+str(Info.get('A',0))+'\t'+str(Info.get('T',0))+'\t'+str(Info.get('?',0))+'\n'	
             c+=1			 
        #  outFas+=Clo.split('bc')[0]+'bc'+str(round(Sup*100))+'\n'+BooClo2Seq[Clo]+'\n'
          outFas+=Clo.replace('Cell','Clone')+'\n'+BooClo2Seq[Clo]+'\n'		  
          out1+=Clo.replace('>','')+'\t'+Clo.replace('>','').replace('Cell','Clone')+'\n'#split('bc')[0]+'bc'+str(round(Sup*100))+'\n'
          for I in IdenLs:
             out1+=I.replace('>','')+'\t'+Clo.replace('>','').replace('Cell','Clone')+'\n'#.split('bc')[0]+'bc'+str(round(Sup*100))+'\n'
       else:
          out1+=Clo.replace('>','')+'\tNone\n'
          for I in IdenLs:
              out1+=I.replace('>','')+'\tNone\n'		  
    GetOut(Fas[:-6]+'_consen.fasta',outFas)		  
    GetOut(Fas[:-6]+'_consen.txt',out)			
    GetOut(Fas[:-6]+'_consen_CellAnno.txt',out1)		
def GetOut(File,In):
    OutF=open(File,'w')
    OutF.write(In)
    OutF.close()	
