import sys
import os
import Functions3
BooTot=sys.argv[1]
BooTot=int(BooTot)+1
OriMeg='FullDatasnv_CloneFinderPlus.meg'#'C:\\Users\\kumarlab\\Desktop\\Boostrap_clone0131\\input\\m8Mseed76\\m8Mseed76CFin_CloneFinder.meg'
BooMeg0='FullData_RepBOOIDsnv_CloneFinderPlus.meg'#'C:\\Users\\kumarlab\\Desktop\\Boostrap_clone0131\\input\\m8Mseed76\\boo\\FullData_RepBOOIDsnv_CloneFinder_map.meg'
if os.path.exists(OriMeg)==True: CloLs,Clo2Seq=Functions3.ReadMegSeq(OriMeg)
else:
   CloLs=[]
   Clo2Seq={} 
#SNVc=len(Clo2Seq[CloLs[0]])
RepLs=list(range(1,BooTot))
#Clo2MutPosOri=Functions3.Meg2PreAbs(OriMeg) #OriMeg[:-4]+'_PreAbs.txt' from 1
#Tu2CloOri,HitCloLsOri,T2C2FOri=Functions3.GetCloHitForTu(OriMeg[:-4]+'.txt',0)
Clo2MutPosLs={}
Tu2CloLs={}
BooC=0
Seq2BooLs={}
for Rep in RepLs:
  BooMeg=BooMeg0.replace('BOOID',str(Rep))
 # print (BooMeg)
  if os.path.exists(BooMeg)==True:
   BooC+=1
   BLs,Bseq=Functions3.ReadMegSeq(BooMeg)
   
   for i in Bseq:
      if i!='#hg19':
         bs=Bseq[i]	  
         Seq2BooLs[bs]=Seq2BooLs.get(bs,[])+[Rep]	  
#   Clo2MutPosBoo=Functions3.Meg2PreAbs(BooMeg[:-4]+'_map.meg') #OriMeg[:-4]+'_PreAbs.txt'
# #  print (Clo2MutPosBoo)
#   for Clo in Clo2MutPosBoo:
#        MutPosLs=Clo2MutPosBoo[Clo]
#        Clo2MutPosLs[Clo]=Clo2MutPosLs.get(Clo,[])+MutPosLs
		
#   Tu2CloBoo,HitCloLsBoo,T2C2FBoo=Functions3.GetCloHitForTu(BooMeg[:-4]+'.txt',0)	
#   for Tu in Tu2CloBoo:
#      CloBoo=Tu2CloBoo[Tu]
#      CloBoo=list(set(CloBoo))
#      Tu2CloLs[Tu]=Tu2CloLs.get(Tu,[])+CloBoo	  
print (BooC)
#print (Clo2MutPosLs)
#out='OriginalClone\tPosition\tAssign\tMutSupportProportion\n'	
#outPreAbs='OriginalClone\tTumor\tOriginalCloFre\tPreSupportProportion\n'  
#for Clo in CloLs:
#  if Clo!='#hg19':  	
#    c=0
#    BooMutPosLs=Clo2MutPosLs[Clo]	
#    while c<SNVc:
#        MutRep=BooMutPosLs.count(c+1)
#        Pro=1.0*MutRep/BooC		
#        out+=Clo.replace('#','')+'\t'+str(c+1)+'\t'+Clo2Seq[Clo][c]+'\t'+str(Pro)+'\n'
#        c+=1
#    for Tu in Tu2CloOri:
#        PreClo=Tu2CloLs[Tu].count(Clo.replace('#',''))
#        ProPre=1.0*PreClo/BooC		
#        outPreAbs+=Clo.replace('#','')+'\t'+Tu+'\t'+str(T2C2FOri[Tu].get(Clo.replace('#',''),0))+'\t'+str(ProPre)+'\n'	
#Functions3.GetOut(OriMeg[:-4]+'_genoBoo.txt',out)
#Functions3.GetOut(OriMeg[:-4]+'_CloPreBoo.txt',outPreAbs)

out='>Normal\n'+('A'*len(bs))+'\n'
g=1
out1='CloneID\tBooRepID\n'
Seq2BooCloID={}
for Seq in Seq2BooLs:
   BooLs=Seq2BooLs[Seq]
   bc=len(Seq2BooLs[Seq])
  # print (BooLs)
   out+='>G'+str(g)+'bc'+str(bc)+'\n'+Seq+'\n'
   Seq2BooCloID[Seq]='G'+str(g)+'bc'+str(bc)
   for i in BooLs:
      out1+= 'G'+str(g)+'bc'+str(bc)+'\t'+str(i)+'\n'  
   g+=1
for C in CloLs: 
   out+='>'+C.replace('#','')+'\n'+Clo2Seq[C]+'\n'  
Functions3.GetOut(OriMeg[:-4]+'_All.fasta',out)   
Functions3.GetOut(OriMeg[:-4]+'_All.txt',out1) 
print ('summarize clone frequency')
out='BooRep\tTumor\tOriginalCloneID\tBooCloneID\tCloneFreq\n'
for Rep in RepLs:
  BooMeg=BooMeg0.replace('BOOID',str(Rep))
  if os.path.exists(BooMeg)==True:
   BLs,Bseq=Functions3.ReadMegSeq(BooMeg)
   Tu2CloBoo,HitCloLsBoo,T2C2FBoo=Functions3.GetCloHitForTu(BooMeg[:-4]+'.txt',0)	
   for T in T2C2FBoo:
        C2FBoo= T2C2FBoo[T]
        for C in C2FBoo:
            F=C2FBoo[C]
            if float(F)>0:
                BooID=Seq2BooCloID[Bseq['#'+C]]
                out+=str(Rep)+'\t'+T+'\t'+C+'\t'+BooID+'\t'+str(F)+'\n'				
Functions3.GetOut(OriMeg[:-4]+'_AllCloFre.txt',out) 		