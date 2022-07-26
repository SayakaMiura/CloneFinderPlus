import sys
import os
import numpy as np

from alignments.MegaAlignment import MegaAlignment
class generate_replicate(object):
    def BooSum(self,OriTreFile,Rpath):
        Align=MegaAlignment()	
        cwd = os.getcwd()
        cwd=cwd.replace('\\','/')
        Rout=''
        Rout+='library(ape)\n'
        Rout+='true_tree <- ape::read.tree(\"'+cwd+'/'+OriTreFile+'\")\n'
        Rout+='lf <- list.files(\"'+cwd+'\", pattern = \"_map.nwk\", full.names = TRUE)\n'
        Rout+='x<- ape::rmtree(length(lf), true_tree$Nnode)  ## length(lf) = num of replicate trees\n'
        Rout+='for(j in 1:length(lf)){\n'
        Rout+='  x[j] <- list(ape::read.tree(lf[j]))\n'
        Rout+='}\n'
        Rout+='b <- phangorn::plotBS(true_tree, x, p =10,  \"phylogram\")\n'
        Rout+='true_boot_support <- as.numeric(b$node.label)\n'
        Rout+='ape::write.tree(b, file = \"'+cwd+'/Tree_bootstrap.nwk\")\n'
        Align.GetOut('RunR.r',Rout)
        os.system(Rpath+' RunR.r')
       # open('A','r').readlines()	
       # os.remove('RunR.r')	
		
	
    def map_clone_sc(self,Ori2Seq,OriCell2Clo,Boo2Seq,BooCell2Clo,Btre,Rep):
      Rep='Rep'+str(Rep)	
      Align=MegaAlignment()
      CloAnno=[]	  
      out='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n'
      Old2NewID={}
      for Ori in Ori2Seq:
         if Ori.replace('#','').replace('>','')!='Normal' and Ori.replace('#','').replace('>','')!='hg19':
             OSeq=Ori2Seq[Ori]
             BestC=9999
             BestSeq=''			 
             for Boo in Boo2Seq:
                 if Boo.replace('#','').replace('>','')!='Normal' and Boo.replace('#','').replace('>','')!='hg19':
                     BSeq=Boo2Seq[Boo]				 
                     DifC=Align.CountDifNum(OSeq,BSeq)
                     if DifC<BestC: 
                        BestC=DifC
                        BestSeq=Boo	
             Old2NewID[BestSeq.replace('#','')]=Old2NewID.get(BestSeq.replace('#',''),[])+[Ori.replace('#','')]
             out+='#'+Ori.replace('#','').replace('>','')+'\n'+Boo2Seq[BestSeq]+'\n'
             		 
      Align.GetOut(Rep+'_map.meg',out)			 
      print ('clone annotation',Old2NewID,len(OriCell2Clo))	
     # print (list(OriCell2Clo.keys())[0])	  
      for BCell in BooCell2Clo:
          BCell=BCell.replace('#','').replace('>','')	
          #print (OriCell2Clo['#'+BCell.replace('SA'+BCell.split('SA')[-1],'')])		  
          Oclo=OriCell2Clo['#'+BCell.replace('SA'+BCell.split('SA')[-1],'')]
          Bclo=	BooCell2Clo['#'+BCell]	  
          CloAnno.append(Rep+'\t'+BCell.replace('SA'+BCell.split('SA')[-1],'')+'\t'+Oclo+'\t'+';'.join(Old2NewID.get(Bclo,['NA']))+'\t'+Bclo)	  
 
      Extra=[]
      for Boo in Boo2Seq:
          Boo=Boo.replace('#','').replace('>','')	  
          if Old2NewID.get(Boo,'NA')=='NA' and Boo!='hg19' and Boo!='Normal': Extra.append(Boo)
      #if TreeEva=='y':		  
      BooTre=open(Btre,'r').readlines()[0]
      for Old in Old2NewID:
	  
          if BooTre.find('('+Old+':')==-1 and BooTre.find(','+Old+':')==-1:
              print ('tree does not have clone',Old,BooTre)
              open('A','r').readlines()
          else:
            NewIDls=Old2NewID[Old]
            if len(NewIDls)==1: In=	'TEMPORARYNAME'+NewIDls[0]		
            elif len(NewIDls)>1: 
                In='('+'TEMPORARYNAME'+NewIDls[0]+':0.0,'+'TEMPORARYNAME'+NewIDls[1]+':0.0)'
                if len(NewIDls)>2:				
                    Oths=NewIDls[2:]
                    for Oth in Oths:
                        In='('+In+':0.0,'+'TEMPORARYNAME'+Oth+')'
          if BooTre.find('('+Old+':')!=-1:
              BooTre=BooTre.replace('('+Old+':','('+In+':')	
              #print ('r','('+Old+':','('+In+':')			  
          elif BooTre.find(','+Old+':')!=-1:						
              BooTre=BooTre.replace(','+Old+':',','+In+':')
              #print ('r',','+Old+':',','+In+':')			  
          else: open('A','r').readlines()	
      print ('extra clone in boo',Extra)	
      Align.prune_tree(BooTre,Extra,Rep+'_map0.nwk')
      BooTre=open(Rep+'_map0.nwk','r').readlines()[0] 	   
      BooTre=BooTre.replace('TEMPORARYNAME','')	  
      Align.GetOut(Rep+'_map.nwk',BooTre)	  
      os.remove(Rep+'_map0.nwk')     
	  
	#  BooTre=BooTre.replace('TEMPORARYNAME','')	  
     # Align.prune_tree(BooTre,Extra,Rep+'_map.nwk') #Rep='Rep'+str(Rep)
    #  Align.GetOut(BooID+'_mapID.txt',str(Old2NewID))	  
      return CloAnno			
		
		
    def map_clone(self,OriMeg,BooID,TreeEva):
      Align=MegaAlignment()
      OriLs,Ori2Seq,out=Align.ReadMegSeq(OriMeg)
      BooLs,Boo2Seq,AA=Align.ReadMegSeq(BooID+'.meg')
      Old2NewID={}
      for Ori in OriLs:
         if Ori!='#Normal' and Ori!='#hg19':
             OSeq=Ori2Seq[Ori]
             BestC=9999
             BestSeq=''			 
             for Boo in BooLs:
                 if Boo!='#Normal' and Boo!='#hg19':
                     BSeq=Boo2Seq[Boo]				 
                     DifC=Align.CountDifNum(OSeq,BSeq)
                     if DifC<BestC: 
                        BestC=DifC
                        BestSeq=Boo	
             out+=Ori+'\n'+Boo2Seq[BestSeq]+'\n'
             Old2NewID[BestSeq.replace('#','')]=Old2NewID.get(BestSeq.replace('#',''),[])+[Ori.replace('#','')]
      print ('clone annotation',Old2NewID)			 
      Align.GetOut(BooID+'_map.meg',out)
      BooCloFre=BooID+'.txt'
      BooCloFre=open(BooCloFre,'r').readlines()
      Head=BooCloFre[0].strip().split('\t')
      CloneC=len(Head)	
#outCloFre=['Rep\tTumor\tRef clone\tRep Clone ID\trep clone frequency']	 
      CloFreLs=[] 
      out=''	  
      for Line in BooCloFre:	  
         Line=Line.strip().split('\t')	  
         out+=Line[0]
         Cc=1 
         while Cc<CloneC:		 
            OriClo =Head[Cc]
            if Line[0]=='Tumor': 
              if Old2NewID.get(OriClo,'NA')!='NA':			
               NewIDls=Old2NewID[OriClo]
               for NewID in NewIDls:		 
                   out+='\t'+NewID
            else:
              if Old2NewID.get(OriClo,'NA')!='NA':			
               NewIDls=Old2NewID[OriClo]
               Info=BooID.split('\\')[-1]+'\t'+Line[0]+'\t'			   
               for NewID in NewIDls:		 
                   out+='\t'+Line[Head.index(OriClo)]
                   CloFreLs.append(Info+NewID+'\t'+Line[Head.index(OriClo)])		##outCloFre=['Rep\tTumor\tRef clone\tRep Clone ID\trep clone frequency']			   
            Cc+=1				   
         out+='\n'
    #  BooCloFre=BooCloFre[1:]
    #  out+='\n'.join(BooCloFre)
	
      Align.GetOut(BooID+'_map.txt',out)	  
      BooTre=BooID+'.nwk'
      Extra=[]
      for Boo in BooLs:
          Boo=Boo.replace('#','')
          if Old2NewID.get(Boo,'NA')=='NA' and Boo!='hg19' and Boo!='Normal': Extra.append(Boo)
      if TreeEva=='y':		  
       BooTre=open(BooTre,'r').readlines()[0]
       for Old in Old2NewID:
	  
          if BooTre.find('('+Old+':')==-1 and BooTre.find(','+Old+':')==-1:
              print ('tree does not have clone',Old,BooTre)
              open('A','r').readlines()
          else:
            NewIDls=Old2NewID[Old]
            if len(NewIDls)==1: In=	'TEMPORARYNAME'+NewIDls[0]		
            elif len(NewIDls)>1: 
                In='('+'TEMPORARYNAME'+NewIDls[0]+':0.0,'+'TEMPORARYNAME'+NewIDls[1]+':0.0)'
                if len(NewIDls)>2:				
                    Oths=NewIDls[2:]
                    for Oth in Oths:
                        In='('+In+':0.0,'+'TEMPORARYNAME'+Oth+')'
          if BooTre.find('('+Old+':')!=-1:
              BooTre=BooTre.replace('('+Old+':','('+In+':')	
              #print ('r','('+Old+':','('+In+':')			  
          elif BooTre.find(','+Old+':')!=-1:						
              BooTre=BooTre.replace(','+Old+':',','+In+':')
              #print ('r',','+Old+':',','+In+':')			  
          else: open('A','r').readlines()	
       print ('extra clone in boo',Extra)
       Align.prune_tree(BooTre,Extra,BooID+'_map0.nwk')
       BooTre=open(BooID+'_map0.nwk','r').readlines()[0] 	   
       BooTre=BooTre.replace('TEMPORARYNAME','')	  
       Align.GetOut(BooID+'_map.nwk',BooTre)
       os.remove(BooID+'_map0.nwk')	   
      Align.GetOut(BooID+'_mapID.txt',str(Old2NewID))	  
      return CloFreLs	
    def GetHeadObsFreqTaHead(self, Head): 
      Head=Head.strip().split('\t')
      SampOrder=[]
      Name2Col={}
      c=0
      Len=len(Head)
      while c<Len:
        i=Head[c]
        if i.find(':')!=-1:
            Samp=i.split(':')[0]
            if Samp not in SampOrder:
                 SampOrder.append(Samp)
            Name2Col[i]=c
        c+=1
      return SampOrder,Name2Col
    def ObsFreqTaHead(self, Freq): 
      SNVOut=Freq.replace('.txt','_SNVfreq.txt')
      Freq=open(Freq,'r').readlines() 
      SampOrder,Name2Col=self.GetHeadObsFreqTaHead(Freq[0])
      Samp2FreqIn={}
      Samp2TotRead={}
      for Samp in SampOrder:
          Samp2FreqIn[Samp]=[]
          Samp2TotRead[Samp]=[]	  
      Freq=Freq[1:]
      outS=''
      SNVnum=0
      for i in Freq:
        i=i.strip().split('\t')
        TMP={}
        for Samp in SampOrder:
            if SNVnum==0: outS+=Samp+'\t'
            MutC=int(i[Name2Col[Samp+':alt']])
            RefC=int(i[Name2Col[Samp+':ref']])
            if (MutC+RefC)<1: MutFreq=0			
            else: MutFreq=1.0*MutC/(MutC+RefC)
            Tot=MutC+RefC
            Samp2FreqIn[Samp].append(MutFreq)
            Samp2TotRead[Samp].append(Tot)		
        SNVnum+=1
      outS=outS[:-1]+'\n'
      c=0
      while c<SNVnum:
        for Samp in SampOrder:
            outS+=str(Samp2FreqIn[Samp][c])+'\t'
        outS=outS[:-1]+'\n'
        c+=1
      #GetOut(SNVOut.split('\\')[-1],outS)
      return SampOrder, Samp2FreqIn	,SNVnum,Samp2TotRead    
    def sample_snvMEG(self, Meg,OutMeg):
        Align=	MegaAlignment()
        SeqLs,SeqDic,out=Align.ReadMegSeq(Meg)
        SNVnum=len(SeqDic[SeqLs[0]])
        #out='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n\n'		
        SNVposls=np.random.choice(SNVnum, SNVnum)
        for S in SeqLs:
           Seq=SeqDic[S]
           New=''
           for i in SNVposls:
             New+=Seq[i]
           out+=S+'\n'+New+'\n'
        self.GetOut(OutMeg,out)	
    def sample_read_snv(self, CFin,BooRep):
        SampOrder, Samp2FreqIn	,SNVnum,Samp2TotRead=self.ObsFreqTaHead(CFin)
        Head=''
        for Samp in SampOrder:
            Head+=Samp+':alt\t'+Samp+':ref\t'
        Head=Head[:-1]+'\n'	
        print (Head,SampOrder)
        Boo=1
        print (SNVnum)
        while Boo<=BooRep: 
            SNVposls=np.random.choice(SNVnum, SNVnum)
            print (len(SNVposls))
          #  open('A','r').readlines()	
            out=Head
            NewPos=0
            outPos='NewPos\tOriginal Pos\n'			
            #S=0
            for S in SNVposls:	
                outPos+=str(NewPos)+'\t'+str(S)+'\n'
                NewPos+=1				
            #while S<SNVnum:
                for Samp in SampOrder:
                     TotRead=Samp2TotRead[Samp][S]
                     Freq=Samp2FreqIn[Samp][S]			 
                     Samp=np.random.choice(2, TotRead, p=[1-Freq,Freq])
                     RefC=list(Samp).count(0)
                     AltC=list(Samp).count(1)
                     out+=str(AltC)+'\t'+str(RefC)+'\t'
                out=out[:-1]+'\n'
              #  S+=1
            self.GetOut(CFin[:-4]+'_Rep'+str(Boo)+'_Pos.txt',outPos)			  
            self.GetOut(CFin[:-4]+'_Rep'+str(Boo)+'.txt',out)		
            Boo+=1	
    def sample_cell(self, BEAMin,BooRep,Redun):
     Align=MegaAlignment()
     StLs, Cell2Seq,Head=Align.ReadMegSeq(BEAMin)
     Len=len(StLs)	 
     Boo=1
     while Boo<=BooRep: 
        out=Head	 
        Samp=list(np.random.choice(StLs,replace=True,size=Len))
        print ('sampling seq',Boo,Len,len(set(Samp)))
       
        SeqID=1	
        if Redun=='n':
            Samp=list(set(Samp))		
        for i in Samp:
                        out+=i+'SA'+str(SeqID)+'\n'+Cell2Seq[i]+'\n'
                        SeqID+=1	
        Align.GetOut('Rep'+str(Boo)+'.meg',out)
        Boo+=1					
    def sample_read(self, CFin,BooRep):

     SampOrder, Samp2FreqIn	,SNVnum,Samp2TotRead=self.ObsFreqTaHead(CFin)
     Head=''
     for Samp in SampOrder:
        Head+=Samp+':alt\t'+Samp+':ref\t'
     Head=Head[:-1]+'\n'	
     print (Head,SampOrder)
     Boo=1
     while Boo<=BooRep: 
        out=Head
        S=0
        while S<SNVnum:
            for Samp in SampOrder:
                 TotRead=Samp2TotRead[Samp][S]
                 Freq=Samp2FreqIn[Samp][S]			 
                 Samp=np.random.choice(2, TotRead, p=[1-Freq,Freq])
                 RefC=list(Samp).count(0)
                 AltC=list(Samp).count(1)
                 out+=str(AltC)+'\t'+str(RefC)+'\t'
            out=out[:-1]+'\n'
            S+=1
        self.GetOut(CFin[:-4]+'_Rep'+str(Boo)+'.txt',out)
        Boo+=1		
    def sample_read_1rep(self, CFin,Boo):

        SampOrder, Samp2FreqIn	,SNVnum,Samp2TotRead=self.ObsFreqTaHead(CFin)
        Head=''
        for Samp in SampOrder:
           Head+=Samp+':alt\t'+Samp+':ref\t'
        Head=Head[:-1]+'\n'	
        print (Head,SampOrder)
    # Boo=1
    # while Boo<=BooRep: 
        out=Head
        S=0
        while S<SNVnum:
            for Samp in SampOrder:
                 TotRead=Samp2TotRead[Samp][S]
                 Freq=Samp2FreqIn[Samp][S]			 
                 Samp=np.random.choice(2, TotRead, p=[1-Freq,Freq])
                 RefC=list(Samp).count(0)
                 AltC=list(Samp).count(1)
                 out+=str(AltC)+'\t'+str(RefC)+'\t'
            out=out[:-1]+'\n'
            S+=1
        self.GetOut(CFin[:-4]+'_Rep'+str(Boo)+'.txt',out)
      #  Boo+=1		
    def generate_CloSeq(self,Boo,MegFile,OutFName,Dir):
       Align=MegaAlignment()
       CloLs,Clo2Seq=Align.name2seq(open(MegFile,'r').readlines())
       NewSeqDic={}
       PosLsF=Dir+'\\FullData_Rep'+str(Boo)+'_Pos.txt'
       Pos=[]
       PosLsF=open(PosLsF,'r').readlines()[1:]
       Ori2New={}	   
       for i in PosLsF:
          Pos.append(int(i.split('\t')[1]))
          NewOri=i.split('\t')
          Ori=int(NewOri[1])
          New=int(NewOri[0])		  
          Ori2New[Ori]=Ori2New.get(Ori,[])+[]		  
       for Clo in CloLs:
           Seq=Clo2Seq[Clo]
           NewSeq=''
           for P in Pos:
              NewSeq+=Seq[P]
           NewSeqDic[Clo]=NewSeq 
       SeqIn=Align.UpMeg(NewSeqDic, CloLs)	
       Align.save_mega_alignment_to_file(OutFName, SeqIn)	   
    def make_consensus_cellSeq(self,BooC,ConsenCut):
        Cell2Pos2NucC={}
        Align=MegaAlignment()		
        B=1	
    #    TotB=0		
        while B<=BooC:		
           BEAM='Rep'+str(B)+'_BEAM.meg'
           if os.path.exists(BEAM)==True:	
          #  TotB+=1		   
            CellL0,CellS0,outMeg=Align.ReadMegSeq(BEAM)
            CellS={}
            for Cell0 in CellS0:
                Cell=Cell0.replace('SA'+Cell0.split('SA')[-1],'')
                CellS[Cell]=CellS0[Cell0]
            CellL=list(CellS.keys())				
            Len=len(CellS[CellL[0]])			
            for Cell in CellL:
               if Cell not in Cell2Pos2NucC: Cell2Pos2NucC[Cell]={}
               Seq=CellS[Cell]
			   
               c=0	
               while c<Len:
                   Base=Seq[c]			   
                   if c not in Cell2Pos2NucC[Cell]: Cell2Pos2NucC[Cell][c]={'A':0,'T':0,'?':0}
                   Cell2Pos2NucC[Cell][c][Base]=Cell2Pos2NucC[Cell][c].get(Base,0)+1	
                   c+=1	
           B+=1				   
 
        out='Cell\tPosition\tA\tT\tConsen\tBCV\n'
        #outMeg='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n'	
        for Cell in Cell2Pos2NucC:
            c=0
            NewSeq=''			
            while c<Len:
                NucDic=Cell2Pos2NucC[Cell][c]
                A=NucDic.get('A',0)
                T=NucDic.get('T',0)
                
                TotB=0	
                for Base in NucDic:
                   TotB+=NucDic[Base]				
                Apro=1.0*A/TotB	
                Tpro=1.0*T/TotB	
                out+=Cell.replace('#','')+'\t'+str(c+1)+'\t'+str(A)+'\t'+str(T)	+'\t'				
                if Apro>ConsenCut:
                    out+='A\t'+str(Apro)+'\n'
                    NewSeq+='A'					
                elif Tpro>ConsenCut:
                    out+='T\t'+str(Tpro)+'\n'
                    NewSeq+='T'	
                else:
                    NewSeq+='?'				
                    if Apro>Tpro: 
                        out+='A?\t'+str(Apro)+'\n'					
                    else: 	
                        out+='T?\t'+str(Tpro)+'\n'
                c+=1	
            outMeg+=Cell+'\n'+NewSeq+'\n'	
        self.GetOut('ConsensusCellSeq.txt',out)
        self.GetOut('ConsensusCellSeq.meg',outMeg)	
    def make_consensus_cellSeq_fas(self,BooC,ConsenCut):
        Cell2Pos2NucC={}
        Align=MegaAlignment()		
        B=1	
    #    TotB=0		
        while B<=BooC:		
           BEAM='Rep'+str(B)+'_BEAM.meg'
           if os.path.exists(BEAM)==True:	
          #  TotB+=1		   
            CellL0,CellS0,outMeg=Align.ReadMegSeq(BEAM)
            CellS={}
            for Cell0 in CellS0:
                Cell=Cell0.replace('SA'+Cell0.split('SA')[-1],'')
                CellS[Cell]=CellS0[Cell0]
            CellL=list(CellS.keys())				
            Len=len(CellS[CellL[0]])			
            for Cell in CellL:
               if Cell not in Cell2Pos2NucC: Cell2Pos2NucC[Cell]={}
               Seq=CellS[Cell]
			   
               c=0	
               while c<Len:
                   Base=Seq[c]			   
                   if c not in Cell2Pos2NucC[Cell]: Cell2Pos2NucC[Cell][c]={'A':0,'T':0,'?':0}
                   Cell2Pos2NucC[Cell][c][Base]=Cell2Pos2NucC[Cell][c].get(Base,0)+1	
                   c+=1	
           B+=1				   
 
        out='Cell\tPosition\tA\tT\tConsen\tBCV\n'
        outMeg=''##MEGA\n!Title SNVs;\n!Format datatype=dna;\n'	
        for Cell in Cell2Pos2NucC:
            c=0
            NewSeq=''			
            while c<Len:
                NucDic=Cell2Pos2NucC[Cell][c]
                A=NucDic.get('A',0)
                T=NucDic.get('T',0)
                
                TotB=0	
                for Base in NucDic:
                   TotB+=NucDic[Base]				
                Apro=1.0*A/TotB	
                Tpro=1.0*T/TotB	
                out+=Cell.replace('#','')+'\t'+str(c+1)+'\t'+str(A)+'\t'+str(T)	+'\t'				
                if Apro>ConsenCut:
                    out+='A\t'+str(Apro)+'\n'
                    NewSeq+='A'					
                elif Tpro>ConsenCut:
                    out+='T\t'+str(Tpro)+'\n'
                    NewSeq+='T'	
                else:
                    NewSeq+='?'				
                    if Apro>Tpro: 
                        out+='A?\t'+str(Apro)+'\n'					
                    else: 	
                        out+='T?\t'+str(Tpro)+'\n'
                c+=1	
            outMeg+=Cell.replace('#','>')+'\n'+NewSeq+'\n'	
        self.GetOut('ConsensusCellSeq.txt',out)
        self.GetOut('ConsensusCellSeq.fasta',outMeg)	
		
    def GetOut(self, OutFileName, In):
        OutF=open(OutFileName,'w')
        OutF.write(In)
        OutF.close()
 