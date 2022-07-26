from alignments.MegaAlignment import MegaAlignment

class CloneFrequencyAnalizer:
    def GetHeadCloneFreq(self,Head):
        Head=Head.strip().split('\t')
        Len=len(Head)
        c=1
        Name2Col={}
        while c<Len:
            Name2Col[Head[c]]=c
            c+=1
        return Name2Col
    def GetCloHitForTu(self,File,Cut):
        File=open(File,'r').readlines()
        Name2Col=self.GetHeadCloneFreq(File[0])
        File=File[1:]
        Tu2Clo={}
        T2C2F={}
        HitCloLs=[]		
        for i in File:
            i=i.strip().split('\t')	
            #Tu=i[0]
            Tu=i[0].replace('T-','')
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
	
    def remove_minor_clone(self,File,Cut,MegFile):
        Tu2Clo,CloLs,T2C2F=self.GetCloHitForTu(File,Cut)
        self.UpCloFreqTa1(T2C2F,CloLs,File[:-4]+'_'+str(Cut)+'.txt')
        Align=MegaAlignment()
        CloLsMeg,Clo2Seq=Align.name2seq(open(MegFile,'r').readlines())		
        NewSeqBuild=Align.UpMeg(Clo2Seq, CloLs)
        Align.save_mega_alignment_to_file(MegFile[:-4]+'_'+str(Cut)+'.meg', NewSeqBuild)		
    def UpCloFreqTa1(self,T2C2F,NameLs,OutN): #NameLs is Clone Ls
        Num=1
        Order=NameLs
        TMP='Tumor'	
        for Name in NameLs:
            Num+=1
            TMP+='\t'+Name
        TMP+='\n'
        TuLs=[]
        for T in T2C2F:
            TuLs.append(T)
        TuLs.sort()		
        for T in TuLs:
            TMP+=T
            for C in Order:
                Code=C in T2C2F[T]
                if Code!=True: In='\t0'			
                else:In='\t'+str(T2C2F[T][C])		
                TMP+=In
            TMP+='\n'			
        OutF=open(OutN,'w')
        OutF.write(TMP)
        OutF.close()	
    def Sort(self, Ls, C2F):
            
            FreLs=[]
            Fre2Clo={}
            for i in Ls:
               if i[0]=='#': i=i[1:]
               Code=i in C2F
               if Code==True:	   
                F=C2F[i]
                Code=F in FreLs
                if Code!=True:
                   FreLs.append(F)
                   Fre2Clo[F]=[]
                Fre2Clo[F].append(i)
            FreLs.sort(reverse=True) 
            TMP=[]
            for F in FreLs:
                TMP+=Fre2Clo[F]
            return TMP,Fre2Clo
			

    def save_frequency_table_to_file(self, file_name, T2C2F, NameLs):
            if NameLs == []:
                for T in T2C2F:
                    C2F=T2C2F[T]
                    for C in C2F:					
                        if C2F[C]>0: NameLs.append(C)
                NameLs=list(set(NameLs))
            TMP='Tumor'	
            for Name in NameLs: #Clone
              	TMP+='\t'+Name				
            TMP+='\n'
            TuLs=[]
            for T in T2C2F:
                TuLs.append(T)
            TuLs.sort()		
            for T in TuLs:
                TMP+=T
                for C in NameLs:
                    Code=C in T2C2F[T]
                    if Code!=True: In='\t0'			
                    else:In='\t'+str(T2C2F[T][C])		
                    TMP+=In
                TMP+='\n'
     
            destination = open(file_name,'w')
            destination.write(TMP)
            destination.close()    