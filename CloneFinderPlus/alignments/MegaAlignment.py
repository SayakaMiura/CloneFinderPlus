
class MegaAlignment():

    def name2seq(self, align_list):
        self.M = align_list	
        self.clone_order=[]
        self.clone_seq = {}
        Name=''		
        for Seq in self.M:
          Seq=Seq.strip()		
          if Seq!='':		
            if Seq[0]=='#' and Seq!='#MEGA':
                if Seq!='#hg19' and Seq!='#Normal':
                    self.clone_order.append(Seq)
                   				
                    self.clone_seq[Seq]=''
                    Name=Seq
                else: Name=''  					
            elif Name!='':
                self.clone_seq[Name] += Seq	
        return self.clone_order, self.clone_seq	

		
    def UpMeg(self, Name2Seq0, NameLs):
            if NameLs==[]:
                 for Name in Name2Seq0:
                     NameLs.append(Name)			 
            out=['#MEGA','!Title SNVs;','!Format datatype=dna;',' ']
            for Name in NameLs:
                if Name in Name2Seq0: Seq=Name2Seq0[Name]
                else: Seq=Name2Seq0['#'+Name]				
                if Name[0]!='#': Name='#'+Name
                				
                out+=[Name,Seq]#Name2Seq0[Name]]
            return out 
    def findLongExtBra(self,SeqLs):
        CloLs, Clo2Seq=self.name2seq(SeqLs)
      #  print (Clo2Seq)		
        Clo2UniPos={}
        Len=len(Clo2Seq[CloLs[0]])
        c=0
        while c<Len:	
            Tls=[]		
            for Clo in CloLs:
                Nuc=Clo2Seq[Clo][c]
                if Nuc=='T': Tls.append(Clo)
            if len(Tls)==1: Clo2UniPos[Tls[0]]=Clo2UniPos.get(Tls[0],[])+[c]
            c+=1
        return Clo2UniPos		
    def SumClo2UniMutLs(self,Clo2MutLs,T2C2F,Cut,MutCut):
        Clo2TuLs={}
        for T in T2C2F:
            C2F=T2C2F[T]
            for C in C2F:
               if C2F[C]>=Cut:
                   Clo2TuLs[C]=Clo2TuLs.get(C,[])+[T]
     #   print (Clo2TuLs)
        Tu2CloLs={}
        for Clo in Clo2MutLs:
	
            if len(Clo2MutLs[Clo])>=MutCut:
                  Clo=Clo.replace('#','')				
                  TuLs=Clo2TuLs.get(Clo,[])
                  if len(TuLs)==1: Tu2CloLs[TuLs[0]]=Tu2CloLs.get(TuLs[0],[])+[Clo]
     #   print (Tu2CloLs)				  
        return (Tu2CloLs)				  
    def GetUniMut(self,C2F,SeqLs,Cut):
       # HitCloLs=[]
        CloLs, Clo2Seq=self.name2seq(SeqLs)
        HitCloSeq={}		
        for C in C2F:
             if C2F[C]>=Cut: 
                 # HitCloLs.append(C)
                  HitCloSeq[C]=Clo2Seq['#'+C]
        HitSeqLs=self.UpMeg(HitCloSeq, list(HitCloSeq.keys()))				  
        Clo2Mut=self.findLongExtBra(HitSeqLs)
      #  print (Clo2Mut)
        return (Clo2Mut)		
    def MapSNVfre(self,Clo2MutLs,SNVfre):
        Clo2SNVfre={}	
        for Clo in Clo2MutLs:
            MutLs=Clo2MutLs[Clo]
            In=[]
            for i in MutLs:
                In.append(SNVfre[i])
            Clo2SNVfre[Clo]=In	
        #print (Clo2SNVfre)
        return Clo2SNVfre		
    def TestDecom(self,Clo2SNVfre,Clo2SNVpos,TarDecomClo,Clo2Fre,CloSeq):
        print ('cut long branch')
        CutSNVFre=Clo2Fre[TarDecomClo]/2
        TarSNVfre=Clo2SNVfre['#'+TarDecomClo]
        TarSNVpos=Clo2SNVpos['#'+TarDecomClo]		
        CutBra={'less':[],'more':[]}
        CutBraMutID={'less':[],'more':[]}		
        i=0	
        Len=len(TarSNVfre)
        while i<Len:
            SNV=TarSNVfre[i]		
        #for SNV in TarSNVfre:
            if SNV>=CutSNVFre:
                 CutBra['more'].append(SNV)
                 CutBraMutID['more'].append(TarSNVpos[i])				 
            else: 
                CutBra['less'].append(SNV)
                CutBraMutID['less'].append(TarSNVpos[i])				
            i+=1			
        print (CutBra)
        Min=9999
        Att=[]		
        if len(CutBra['more'])>=3 and len(CutBra['less'])>=3:		
         MoreAve=sum(CutBra['more'])/len(CutBra['more'])
         LessAve=sum(CutBra['less'])/len(CutBra['less'])
         MoreLessDif=(MoreAve-LessAve)*(MoreAve-LessAve)
         print ('more less diff',MoreLessDif)		 
         TotCloFre=0
         for C in Clo2Fre:
            TotCloFre+=Clo2Fre[C]
         if TotCloFre<1: TotCloFre=1			

         for Clo in Clo2Fre:
            if Clo!=TarDecomClo:
                   Aave=(Clo2Fre[Clo]/TotCloFre)/2
                   for i in CutBra:
                       snv=CutBra[i]
                       Ave=1.0*sum(snv)/len(snv)
                       Dif=(Aave-Ave)*(Aave-Ave)
                       if Dif<Min and Dif<MoreLessDif:
                                Min=Dif
                                Att=[Clo,i]
         print ('minimum',Att)
        if Att==[]: return []
        else:
            CloLs, Clo2Seq=self.name2seq(CloSeq)		
            MoveSNVID=CutBraMutID[Att[1]]
            FixCloSeq=Clo2Seq['#'+Att[0]]
            NewSeq=''
            OriTarSeq=Clo2Seq['#'+TarDecomClo]
            NewTarSeq=''			
            Len=len(FixCloSeq)
            c=0	
            while c<Len:
               if MoveSNVID.count(c)!=0: 
                   NewSeq+='T'
                   NewTarSeq+='A'				   
               else:
                   NewSeq+=FixCloSeq[c]
                   NewTarSeq+=OriTarSeq[c]				   
               c+=1
           # print (NewSeq,NewTarSeq)			   
            return [NewSeq,NewTarSeq]			   
	
    def FilterLostSNV(self,SNVpos,Tu2Freq):
        Fill=[]
        for i in SNVpos:
            Cou=0	
            for Tu in Tu2Freq:
                if Tu2Freq[Tu][i]>0: Cou+=1	
            if Cou==1: Fill.append(i)
        return Fill
    def AssignLostSNV(self,SNVpos,Tu2SNVfre,T2C2F,Cut,SeqLs):
        NewSeqLs=[]
        OldCloLs,OldClo2Seq=self.name2seq(SeqLs)
        Len=len(OldClo2Seq[OldCloLs[0]])		
        for T in Tu2SNVfre:
            Clo2AddM={}	
            Clo2Fre=T2C2F['T-'+T]	
            #Clo2Fre=T2C2F[T]			
            for SNV in SNVpos:
                Fre=Tu2SNVfre[T][SNV]*2			
                if Fre>0: 
                    Best=99999
                    Bclo=''
                    for Clo in Clo2Fre:
                        Cfre=Clo2Fre[Clo]
                        if Cfre>=Cut:
                            Dif=(Cfre-Fre)*(Cfre-Fre)
                            if Dif<Best:
                                Best=Dif	
                                Bclo=Clo
                    Clo2AddM[Bclo]=Clo2AddM.get(Bclo,[])+[SNV]
            if Clo2AddM!={}:
              for OldClo in Clo2AddM:			
                OldSeq=OldClo2Seq['#'+OldClo]
                AddM=Clo2AddM[OldClo]				
                NewSeq=''
                c=0	
                while c<Len:
                    if AddM.count(c)!=0: NewSeq+='T'
                    else: NewSeq+=OldSeq[c]
                    c+=1	
                print ('clone seq updated by adding lost SNV',T,OldClo,OldSeq,NewSeq)					
                NewSeqLs.append(NewSeq)					
        return NewSeqLs          				
             		
	
    def remove_redund_seqs(self, Meg):
        	
        print('removing redundant seqs...')
        	
        NameOrder, Name2Seq=self.name2seq(Meg)
        Name2IdenLs={}
        for Ref in NameOrder:		
      
            RefSeq=Name2Seq[Ref]
            Name2IdenLs[Ref]=self.find_redundant(RefSeq,Name2Seq)			

        Done=[]
        Nonredun_seqdic={}		
        for Name in NameOrder:
         
          if Done.count(Name)==0:
            IdenLs=Name2IdenLs[Name]
            Nonredun_seqdic[Name]=Name2Seq[Name]				
            Done+=IdenLs
        out2=self.UpMeg(Nonredun_seqdic,[])        
        return out2 
		
    def CountDifNum(self, Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif	
			
    def find_redundant(self,Seq,Seq_dic):
           Redun=[]
           for Name in Seq_dic:
                Seq1=Seq_dic[Name]
                DifNum=self.CountDifNum(Seq,Seq1)
                if DifNum==0: Redun.append(Name)
           return Redun				
			

		
    def ModSeq(self, CSeq0,ChangePosi,ChanNuc,Len):	
                      
                      c=0
                      CutCloSeq=''					  
                      while c<Len:
                          Code1=c in ChangePosi						  
                          if Code1==True: CutCloSeq+=ChanNuc
                          else: CutCloSeq+=CSeq0[c]
                          c+=1
                      return CutCloSeq	
					  
    def GetMutPos(self, Seq):
      TMP=[]
      Len=len(Seq)
      c=0
      while c<Len:
        if Seq[c]=='T': TMP.append(c)
        c+=1
      return TMP 
    def GetSharePosi1(self,Ori2Seq,ShareNuc):
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
	
    def get_mega_alignment_string(self, SeqLs):
        result = ''
        for item in SeqLs:#self._mega_seqs:
            result += item + "\n"
        return result
    
    def save_mega_alignment_to_file(self, filename, SeqLs):
        	
        destination = open(filename,'w')
        destination.write(self.get_mega_alignment_string(SeqLs))
        destination.close()  



    def AddNormal(self, seqs): 
        CellOrder, Cell2Seq = self.name2seq(seqs)	
        Len=len(Cell2Seq[CellOrder[0]])	
        seqs+=['#Normal\n',('A'*Len)+'\n']
        return seqs		

