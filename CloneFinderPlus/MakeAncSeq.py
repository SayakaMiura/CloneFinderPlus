Ref='Normal'
import sys
Ta=sys.argv[1]
NodeMap=Ta.replace('.csv','_nodeMap.txt')
#OriMeg=Ta.replace('.csv','.meg')
Out=Ta.replace('.csv','WithAnc.fasta')
NodeMap=open(NodeMap,'r').readlines()
NodeMap=NodeMap[1:]
Node2Label={}
RefNodeID=''
Dec2Anc={}
Label2Count={}
Label2Seq={}

for i in NodeMap:
    i=i.strip().split(' ')
    ii=[]
    for Item in i:
        Item=Item.strip()
        if Item!='': ii.append(Item)
    Label=ii[0]
    Node=ii[1]
    Decs=[ii[2],ii[3]]
    if Label[:5]=='Node_': 
        for Dec in Decs:
            Dec2Anc[Dec]=Node
        Label2Count[Label]=0
    Label2Seq['Node_'+Node]=''
    if Label==Ref or Label=='Normal': RefNodeID=Node
    Node2Label[Node]=Label
#print (Label2Seq)
def GetBestNuc(Dic):
    TMP={}
    for Node in Dic:
        Nuc2Prob=Dic[Node]
        BestP=0
        BestNuc='A'
        for Nuc in Nuc2Prob:
            P=Nuc2Prob[Nuc]
            if P>BestP:
                BestNuc=Nuc
                BestP=P
        TMP[Node]=BestNuc
    return TMP
	
RefAnc=Dec2Anc[RefNodeID]

RmNode=Node2Label[RefAnc]
if RefAnc in Dec2Anc: print ('Error: Root is incorrect. Fix the nwk file and redo ancestor analysis', Ta)
else:
    AddNode=[]
    Ta=open(Ta,'r').readlines()
    Head=Ta[0].strip().split(',')
    Node2Col={}
    Node2Prob={}
    c=0
    Len=len(Head)
    while c<Len:
        i=Head[c].strip()
        Code=i in Label2Seq
        if Code==True: Node2Col[i]=c
        c+=1
    Ta=Ta[1:]

    PosiP=1
    for i in Ta:
        i=i.strip().split(',')
        PosiC=int(i[0])
        Nuc=i[1]
        if PosiC!=PosiP:
            PosiP=PosiC
            Node2BestNuc=GetBestNuc(Node2Prob)
            for Node in Node2BestNuc:
                Label2Seq[Node]+=Node2BestNuc[Node]
            Node2Prob={}
        for Node in Node2Col:
            Prob=i[Node2Col[Node]]
            Code=Node in Node2Prob
            if Code!=True: Node2Prob[Node]={}
            Node2Prob[Node][Nuc]=float(Prob)

Node2BestNuc=GetBestNuc(Node2Prob)
for Node in Node2BestNuc:
                Label2Seq[Node]+=Node2BestNuc[Node]
out='' 
for Label in Label2Seq:
 #  if Label!=RmNode:
      out+='>'+Label+'\n'+Label2Seq[Label]+'\n'
OutF=open(Out,'w')
OutF.write(out)
OutF.close()