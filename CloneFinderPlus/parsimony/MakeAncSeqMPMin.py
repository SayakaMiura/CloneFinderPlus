from alignments.MegaAlignment import MegaAlignment

class MakeAncSeqMPMin(object):
    
    def __init__(self):
        self._filenames = []
        self.ID = ''
        
    def GetRel(self, Ls):
        Anc2Seq={}	
        for Li in Ls:
           if Li!=[]:        	   
            Anc=Li[0]
            NodeN=Li[1]
            if NodeN!='hg19': 
             if NodeN=='-':
                NodeN=Anc
             SeqIn=''
             Nls=Li[4:]			
             for N in Nls:
                 if len(N)!=1: N='A'			 
                 SeqIn+=N
           #  print (NodeN,SeqIn)				 
             if SeqIn.find('T')!=-1:				 
                 Anc2Seq['#'+NodeN]=SeqIn            	
        return Anc2Seq  
        
    def get_best_alignment(self, files, ID, remove_redundant, tree_list): 
        Align = MegaAlignment()	
        AncFile = files[0]
        print('processing file: ' + AncFile)
        AncFile=open(AncFile,'r').readlines()
        Lines=[]
        Add='n'	
        for i in AncFile:
                if i.find('Index ')!=-1: Add='y'
                elif Add=='y': 
                    i=i.strip().split(' ')
                    In=[]			
                    for ii in i:
                        if ii!='': In.append(ii)
                    						
                    Lines.append(In)
        #print (Lines)					
        Anc2Seq=self.GetRel(Lines)   
       # print (Anc2Seq)		
        outseq=	Align.UpMeg(Anc2Seq,[])	
        outseq = Align.remove_redund_seqs(outseq)	
      #  print (outseq)
      #  open('a','r').readlines()		
        return outseq, tree_list[0]#, best_outset[1], best_outset[2], best_outset[4] #NadeMapInfo, mask_seq, Good_posi_info]
    
    
    def get_all_alignment(self, files, ID, remove_redundant, tree_list): 
      Align = MegaAlignment()
      SeqCol=[]	  
      for AncFile in files:	  
       # AncFile = files[0]
        print('processing file: ' + AncFile)
        AncFile=open(AncFile,'r').readlines()
        Lines=[]
        Add='n'	
        for i in AncFile:
                if i.find('Index ')!=-1: Add='y'
                elif Add=='y': 
                    i=i.strip().split(' ')
                    In=[]			
                    for ii in i:
                        if ii!='': In.append(ii)
                    						
                    Lines.append(In)
        #print (Lines)					
        Anc2Seq=self.GetRel(Lines)
        for A in Anc2Seq:
            Seq=Anc2Seq[A]
            SeqCol.append(Seq)			
       # print (Anc2Seq)
    #  print (len(SeqCol))	   
      SeqCol=list(set(SeqCol))
     # print (len(SeqCol))		  
      Anc2Seq={}
      ID=1
      for i in SeqCol:
         Anc2Seq[str(ID)]=i	
         ID+=1	
    #  print (Anc2Seq)	
    #  open('a','r').readlines()	  
      outseq=	Align.UpMeg(Anc2Seq,[])	
      outseq = Align.remove_redund_seqs(outseq)	
      #  print (outseq)
      #  open('a','r').readlines()		
      return outseq, tree_list[0]#, best_outset[1], best_outset[2], best_outset[4] #NadeMapInfo, mask_seq, Good_posi_info]
    
