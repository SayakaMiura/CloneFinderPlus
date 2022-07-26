import sys
import pydot

PF=sys.argv[1]
BooCut=0.5  #0.7 for boo 0.4 for consen
Out=PF[:-4]+'_BooCut'+str(BooCut)+'.gv'
Data=Out.split('\\')[-1][:-3]
out='digraph D {\nlabel = \"'+Data+'\";\n'
out+='labelloc = \"t\";\n'
PF=open(PF,'r').readlines()[1:]
for i in PF:
    i=i.split('\t')
    Edge=i[0].split('[')[0]
    Boo=float(i[-1])
    out+=Edge+' '
    B=round(100*Boo,0)
	
    if Boo>=BooCut:	
        out+='[label=\"'+str(B).split('.')[0]+'\",fontsize = 8]\n'
    else:
        out+='[color=grey, label=\"'+str(B).split('.')[0]+'\",fontcolor=grey,fontsize = 8]\n'	
out+='}\n'
OutF=open(Out,'w')
OutF.write(out)
OutF.close()

(graph,) = pydot.graph_from_dot_file(Out)
graph.write_png(Out[:-3]+'.png')