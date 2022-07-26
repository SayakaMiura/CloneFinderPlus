from Bio import Phylo
#from Bio.Phylo.Consensus import *
from io import StringIO
import os
class MegaAlignment():
    def find_marker(self,MAO):
        Align = MegaAlignment()
        tree_builder = MegaML(MAO)
        tree_analyzer = TreeAnalizer() 		
        status = tree_builder.do_mega_ml(self.InMeg, 'Noresun')
        if status == True:
            tree1 = tree_builder.newick_trees
        else:
            print ('failed to run megaML')
        Tree_Rooted1 = tree_analyzer.RootTree_rootBottom(tree1, 'Normal')
	
        InferAncestor = MegaAncestor()
        InferAncestor.alignment_file = self.InMeg
        InferAncestor.input_tree_file = Tree_Rooted1	
     			
        self.ancestor_states, self.offspring2ancestor, cell2code, self.code2cell = InferAncestor.retrieve_ancestor_states()
        self.Marker=[] #count from 0	
        out='Count from 0. In subset\n'		
        c=0
        while c<self.SNVnum:
            MutBra=[]
            Back='n'
            for dec in self.offspring2ancestor:
              if dec!='Des1' and dec!='Des2':			
                anc=self.offspring2ancestor[dec]
                #print dec,anc,self.ancestor_states				
                decNuc=	self.ancestor_states['Node_'+dec][c].split('\t')[0]
                ancNuc=	self.ancestor_states['Node_'+anc][c].split('\t')[0]
               # print decNuc,ancNuc
                if decNuc=='T' and ancNuc=='A': MutBra.append(dec)
                if decNuc=='A' and ancNuc=='T': Back='y'
            if Back=='n' and len(MutBra)==1: 
                self.Marker.append(c)
                out+=str(c)+'\n'				
            c+=1
        Align.GetOut(self.ID[:-4]+'_markerPos.txt',out)			
       # print 'marker',self.Marker
    def make_consensus(self,CellLs,OriCell2Seq):
        Len=len(OriCell2Seq[CellLs[0]])
        Con=''
        c=0	
        while c<Len:
           Nuc2C={}
           for Cell in CellLs:
               Nuc=OriCell2Seq[Cell][c]
               Nuc2C[Nuc]=Nuc2C.get(Nuc,0)+1			   
           #    if Nuc not in Nuc2C.has_key(Nuc)!=True: Nuc2C[Nuc]=0	
            #   Nuc2C[Nuc]+=1	
           Best=0	
           BestNuc=[]
           for Nuc in Nuc2C:
               C=Nuc2C[Nuc]
               if C==Best: BestNuc.append(Nuc)
               elif C>Best:
                  Best=C	
                  BestNuc=[Nuc]
           if len(BestNuc)!=1: Con+='?'
           else: Con+=BestNuc[0]
           c+=1	
        return Con	
    def AnnotateClone(self,Meg):	
        CellLs, Cell2Seq,AA=self.ReadMegSeq(Meg)
        Seq2CellLs0={}
        for Cell in CellLs:
           Seq=Cell2Seq[Cell]
           Seq2CellLs0[Seq]=Seq2CellLs0.get(Seq,[])+[Cell]
        print ('clone C',len(Seq2CellLs0))#,len(CellLs))
        SeqLs=list(Seq2CellLs0.keys())
        Len=len(SeqLs)
        Seq2CellLs={}
        Done=[]
        c=0
        Max=Len-1	
        AllCell=0		
        while c<Len:
           if Done.count(c)==0:
               Seq0=SeqLs[c]
               NewCellLs=Seq2CellLs0[Seq0]			   
               c1=c+1
               while c1<Len:
                  if Done.count(c1)==0:
                      Seq1=SeqLs[c1]
                      DifC=self.CountDifNum_exMiss(Seq0,Seq1)
                      if DifC==0:
                          Done.append(c1)
                          NewCellLs+=Seq2CellLs0[Seq1]
                  c1+=1
				  
               Seq2CellLs[Seq0]=NewCellLs
               AllCell+=len(NewCellLs)			   
           c+=1			   
        print ('clone C merge',len(Seq2CellLs))#,AllCell)		
        ID=1
        Clo2Seq={}
        Clo2CellL={}
        Cell2Clo={}	
        Len=len(Seq)		
        out=''		
        for Seq in Seq2CellLs:
           CellLs=Seq2CellLs[Seq]
           CloID='Clone'+str(ID)
           out+='>'+CloID+'\n'+Seq+'\n'		   
           Clo2Seq[CloID]=Seq	
           Clo2CellL[CloID]=CellLs	
           for Cell in CellLs:
              	Cell2Clo[Cell]=CloID	
           ID+=1     				
		
        return Clo2Seq,Clo2CellL,Cell2Clo,out
    def annotate_cell1(self,OriSeq,OutFileID):
        self.ID=OutFileID	 
        self.find_marker()
        self.identify_lineage2()		
        OriCellLs, OriCell2Seq = Align.name2seq(OriSeq)
        #make consensus seq	
        out='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n'	
        outRef='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n'
        outRef_Simple='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n'	
        for Lin in self.Lin2CellLs:
          CellLs=self.Lin2CellLs[Lin]	
          Consen=Align.make_consensus(CellLs,OriCell2Seq)
          LinS=self.Lin2seq[Lin]		  
          for Cell in CellLs:
            outRef+=Cell+'_'+Lin.replace('#','')+'\n'+LinS+'\n'		  
            out+=Cell+'\n'+Consen+'\n'
            outRef_Simple+=Cell+'\n'+LinS+'\n'				
        Align.GetOut(OutFileID[:-4]+'_consenLin.meg',out)			
        Align.GetOut(OutFileID[:-4]+'_Marker_withLinID.meg',outRef)	
        Align.GetOut(OutFileID[:-4]+'_Marker.meg',outRef_Simple)
    def prune_tree(self,Tst,ExtraLs,OutF):
       #TreeLs=open(OriNwk,'r').readlines()
       #Cou=1 
       #for Tst in TreeLs:
        tree=Phylo.read(StringIO(Tst.strip()), "newick")   
      #  tree = Phylo.read(OriNwk, 'newick')
        for tip in ExtraLs:	
             tree.prune(tip)
        Phylo.write(tree, OutF,'newick')
        #Cou+=1
    def ReadFasSeq(self,Fas):
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

		
    def ReadMegSeq(self,Meg): 
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
    def CountDifNum_exMiss(self, Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c] and Seq0[c]!='?' and Seq1[c]!='?': Dif+=1
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
    def GetHead(self,Head):
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
    def GetCloHitForTu(self,File,Cut):
        File=open(File,'r').readlines()
        AA,Name2Col=self.GetHead(File[0])
        File=File[1:]
        Tu2Clo={}
        T2C2F={}	
        for i in File:
            i=i.strip().split('\t')	
            Tu=i[0]
            Tu2Clo[Tu]=[]
            T2C2F[Tu]={}		
            for Name in Name2Col:
              if Name!='Tumor':		
                Fre0=i[Name2Col[Name]]		
                Fre=float(Fre0)
                if Fre0.find('-')==-1 and Fre>Cut:
                     T2C2F[Tu][Name]=Fre
                     Tu2Clo[Tu].append(Name)				 
                else: T2C2F[Tu][Name]=0			
        return Tu2Clo,Name2Col,T2C2F
	
		
    def addGroup(self,CFMeg):
        CFfre=CFMeg[:-4]+'.txt'		
        CloL,Clo2Seq,outM=self.ReadMegSeq(CFMeg)
        Tu2Clo,Name2Col,T2C2F=self.GetCloHitForTu(CFfre,0)
      #  print (T2C2F)	
      #  open('A','r').readlines()		
        for Tu in Tu2Clo:
           CloLs=Tu2Clo[Tu]
           for Clo in CloLs:
              Seq=Clo2Seq['#'+Clo]
              outM+='#'+Clo+'TuSample_'+Tu+'{'+Tu+'}\n'+Seq+'\n'
           if len(CloLs)==1:
              outM+='#'+Clo+'TuSampleDup_'+Tu+'{'+Tu+'}\n'+Seq+'\n'		   
        self.GetOut(CFMeg[:-4]+'_group.meg',outM)
    def GetDiv(self,MegOut):	  
        DivFile=open(MegOut,'r').readlines()[-1].strip().split(' ')
        Div=[DivFile[0],DivFile[-1]]
        return Div		
    def ComputeDiv(self,CFMeg,Boo):
        self.addGroup(CFMeg)	#CFMeg[:-4]+'_group.meg'	
        Add=''	
        if Boo=='Y':
            Add='_boo'
        os.system('megacc -a distance_estimation_avg_in_sub_pops_nucleotide'+Add+'.mao -d '+CFMeg[:-4]+'_group.meg -o Sub.out')
        os.system('megacc -a distance_estimation_avg_overall_pops_nucleotide'+Add+'.mao -d '+CFMeg[:-4]+'_group.meg -o Pop.out')
        os.system('megacc -a distance_estimation_inter_pop_diversity_nucleotide'+Add+'.mao -d '+CFMeg[:-4]+'_group.meg -o Inter.out')
        os.system('megacc -a distance_estimation_proportion_inter_pop_diversity_nucleotide'+Add+'.mao -d '+CFMeg[:-4]+'_group.meg -o InterPro.out')		
        PiT=self.GetDiv('Pop.txt')	
        PiS=self.GetDiv('Sub.txt')			
        PiD=self.GetDiv('Inter.txt')	
        PiN=self.GetDiv('InterPro.txt')			
        os.remove('Pop.txt')
        os.remove('Sub.txt')
        os.remove('Inter.txt')
        os.remove('InterPro.txt')		
        os.remove('Pop_summary.txt')
        os.remove('Sub_summary.txt')
        os.remove('Inter_summary.txt')	
        os.remove('InterPro_summary.txt')			
       # print (PiT,PiS,PiD)
       # open('A','r').readlines()		
        return PiT,PiS,PiD,PiN	
    def Fas2Meg(self,Fas,Meg):
        SeqLs,SeqDic=self.ReadFasSeq(Fas)
        out='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n\n'
        for Seq in SeqLs:
            out+=Seq.replace('>','#')+'\n'+SeqDic[Seq]+'\n'
        self.GetOut(Meg,out)			
	
    def GetOut(self,OutFile,outIn):
     OutF=open(OutFile,'w')
     OutF.write(outIn)
     OutF.close() 	