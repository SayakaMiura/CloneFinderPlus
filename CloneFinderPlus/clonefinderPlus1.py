# 'python CloneFinder.py snv [snv_input] [FastClone.meg]'
from parsers.DefaultTSPParser import DefaultTSPParser
from regression.CloneFrequencyComputer_cnv1 import CloneFrequencyComputer_cnv1
from decomposition.SNPGroupCombiner_cnv1 import SNPGroupCombiner_cnv1
from parsers.DefaultCNVarser import DefaultCNVarser
from config.ParamsLoader import ParamsLoader
from config.FormatInput import FormatInput
from alignments.FreqToMegaSeq import FreqToMegaSeq
from alignments.MegaAlignment import MegaAlignment
from output.CloneFrequencyAnalizer import CloneFrequencyAnalizer
from parsimony.MegaMP import MegaMP
from parsimony.TreeAnalizer import TreeAnalizer
from parsimony.MegaAncestor import MegaAncestor
from tsp_profiles.tsp_information import tsp_information
from significance_test.cluster_test import cluster_test
from output.OutputWrite import OutputWrite
import os
import sys
#import datetime

Significant_cutoff=0.05

#startTime = datetime.datetime.now()
#print (startTime)
dir = os.getcwd()
tree_builder = MegaMP()
tree_builder.mao_file = dir + '/infer_MP_nucleotide.mao'
try: 
    loader = ParamsLoader()
    loader.params_file = 'options.ini' #options_file
    params = loader.load_params()
    total_read_cut = params.total_read_cut
    mutant_read_cut = params.mutant_read_cut	
    summary_file= "input data file: "+params.input_data_file+'\n'+"total read count cutoff: "+str(params.total_read_cut)+'\n'"mutant read count cutoff: "+str(params.mutant_read_cut)+'\n'"clone frequency cutoff: "+str(params.freq_cutoff)+'\n\n'
except:
    print ('Errors in options.ini')	
parser = DefaultTSPParser()
parser.input_data_file = params.input_data_file
parser_cnv = DefaultCNVarser() 
AnalyzeTree = TreeAnalizer()
OutFile = OutputWrite()
Align = MegaAlignment()
CloFreAna = CloneFrequencyAnalizer()
Format = FormatInput()
#path = "Data/"+params.input_id

if parser.parse() == False:
        print (parser.messages)
else:
       tsp_list = parser.get_tumor_sample_profile_list() 			 
       CNV_information = parser_cnv.get_tumor_cnv_profile(params.cnv_data_file)		
       original_align_builder = FreqToMegaSeq()
       original_align_builder.initialize(tsp_list, True) #True to remove duplicates; False to keep duplicates
       num_sites = tsp_list.num_read_counts()
       all_tsp = tsp_information(tsp_list)    				
       v_obs = all_tsp.tumor2alt_frequency()	   
       total_read, alt_read, ReadCount_table = all_tsp.tumor2read_count()				
       tumor_seqs = original_align_builder.get_mega_allalignment()
       initial_seq0 = original_align_builder.get_mega_alignment()
       initial_seq, A1, A2  = OutFile.ReNameCloFreMeg(initial_seq0, {}, 'number') 
       if len(sys.argv)==4:
           IniCloSeq=sys.argv[3]
           initial_seq=open(IniCloSeq,'r').readlines()		   
     #  print (initial_seq)
     #  open('a','r').readlines()	   
       CNV_information_test = Format.add_low_quality_SNV_info(CNV_information,total_read, alt_read,total_read_cut,mutant_read_cut)          	
       print ('add initial ancestor')	
    #   print (tumor_seqs)
       TuLs,Tu2Seq=Align.name2seq(tumor_seqs)	  
     #  print (TuLs,Tu2Seq)	   
     #  open('a','r').readlines()	   
       id = 'mega_alignment' # id will be used internally for file names
       status = tree_builder.do_mega_mp(initial_seq, id)
       	   
       if status != True:
                print ('CloneFinder cannot make inferences, because samples do not have evolutionary structure.')		
       else:
        seqs_with_ancestor, A1 = tree_builder.alignment_least_back_parallel_muts()  	
        Repeat='y'	
        ItC=0		
        while Repeat=='y' and ItC<=5:	
            print ('compute clone frequencies with ancestor')	
           		
            clone_frequency_cnv = CloneFrequencyComputer_cnv1(seqs_with_ancestor, v_obs, CNV_information_test, params.freq_cutoff)
            clone_frequency_cnv.regress_cnv()	
          #  print (clone_frequency_cnv.Tumor2Clone_frequency)	
          #  CloFreAna.save_frequency_table_to_file(IniCloSeq[:-4]+'WithAnc.txt', clone_frequency_cnv.Tumor2Clone_frequency, [])
           # OutFile.GetOut(IniCloSeq[:-4]+'WithAnc.meg','\n'.join(seqs_with_ancestor))	

            print ('assign lost mutations')
            seqs_addedMut = CloneFrequencyComputer_cnv1(clone_frequency_cnv.hitclone_seq_builder, v_obs, CNV_information_test, params.freq_cutoff)			
            seqs_addedMut=seqs_addedMut.AddMut(clone_frequency_cnv.Tumor2Clone_frequency,Tu2Seq,0.01)			
          #  print (seqs_addedMut)
            print ('recompute clone frequency')			
            clone_frequency_cnv = CloneFrequencyComputer_cnv1(seqs_addedMut, v_obs, CNV_information_test, params.freq_cutoff)
            clone_frequency_cnv.regress_cnv()	
          #  print (clone_frequency_cnv.Tumor2Clone_frequency)	
           # CloFreAna.save_frequency_table_to_file(IniCloSeq[:-4]+'AddMiss.txt', clone_frequency_cnv.Tumor2Clone_frequency, [])
           # OutFile.GetOut(IniCloSeq[:-4]+'AddMiss.meg','\n'.join(clone_frequency_cnv.hitclone_seq_builder))	
			
          #  open('a','r').readlines()

		   
            print ('decompose incorrect candidate clone genotypes')
           # print (clone_frequency_cnv.hitclone_seq_builder,clone_frequency_cnv.Tumor2Clone_frequency)
            	
            print ('find long external branch (single tumor specific)')
            Clo2UniMutLs=Align.findLongExtBra(clone_frequency_cnv.hitclone_seq_builder) #count from 0
          #  print (Clo2UniMutLs)
            Tu2LongExtBraCloLs_SingleTu=Align.SumClo2UniMutLs(Clo2UniMutLs,clone_frequency_cnv.Tumor2Clone_frequency, params.freq_cutoff,10) #Long branch:>=10 mut
          #  print (Tu2LongExtBraCloLs_SingleTu)
            RmCloLs=[]
            AddClo={}
            NewCloId=1			
            for Tu in Tu2LongExtBraCloLs_SingleTu:
                   DecomCloLs=Tu2LongExtBraCloLs_SingleTu[Tu]
                   if len(DecomCloLs)==1: #target clone [0]
                       Clo2UniMut=Align.GetUniMut(clone_frequency_cnv.Tumor2Clone_frequency[Tu],clone_frequency_cnv.hitclone_seq_builder, params.freq_cutoff) #count from 0 					   
                     #  print (Clo2UniMut)			
                       if len(Clo2UniMut)>1:
                          print ('try to decompose long branch:',Tu)
                          Clo2SNVfre=Align.MapSNVfre(Clo2UniMut,v_obs[Tu[2:]])
                        #  print (Clo2SNVfre)
                          DecomSeqLs=Align.TestDecom(Clo2SNVfre,Clo2UniMut,DecomCloLs[0],clone_frequency_cnv.Tumor2Clone_frequency[Tu],clone_frequency_cnv.hitclone_seq_builder)
                          if DecomSeqLs!=[]:
                                RmCloLs.append(DecomCloLs[0])
                                AddClo['#N'+str(NewCloId)]=DecomSeqLs[0]
                                AddClo['#N'+str(NewCloId+1)]=DecomSeqLs[1]
                                NewCloId+=2	
                                summary_file+=Tu+' was decomposed\n'								
						  
          #  print (RmCloLs,AddClo)
          #  open('a','r').readlines()		   
           # decomposition = SNPGroupCombiner_cnv1(clone_frequency_cnv.hitclone_seq_builder, v_obs, clone_frequency_cnv.Tumor2Clone_frequency, CNV_information, params.freq_cutoff, tumor_seqs)			
           # hit_seq_dic, Message = decomposition.get_decomposed_seq()			
            #print (Message) 

           			
            if RmCloLs==[]:#Message[:10]!='decomposed' or DecomTuNum>(1.0*len(v_obs)/2):  
                        Repeat='n'                      
                        final_seq = clone_frequency_cnv.hitclone_seq_builder
                        final_clofre = clone_frequency_cnv.Tumor2Clone_frequency
                       # print (seqs_with_ancestor,clone_frequency_cnv.Tumor2Clone_frequency)	
                       # open('a','r').readlines()						
                        print ('no more new clones')#,clone_frequency_cnv.Tumor2Clone_frequency)
				
            else: 
                       OriCloLs,OriCloSeqDic=Align.name2seq(clone_frequency_cnv.hitclone_seq_builder)
                       for Oclo in OriCloLs:
                          if RmCloLs.count(Oclo.replace('#',''))==0: AddClo[Oclo]=OriCloSeqDic[Oclo]					   
                       AddClo['#hg19']='A'*num_sites
					   
					   
					   
                       hit_seq_build = Align.UpMeg(AddClo,list(AddClo.keys()))		
                    #   print (hit_seq_build)
                     #  open('a','r').readlines()					   
                       id = 'dec_mega_alignment' 				   
                       status = tree_builder.do_mega_mp(hit_seq_build, id)
                       if status == True :
                             decomseqs_with_ancestor, A1 = tree_builder.alignment_least_back_parallel_muts() # True will execute code to remove redundant seqs (True is default)
							# self._cleanup_temp_files()	
                             clone_frequency_cnv2 = CloneFrequencyComputer_cnv1(decomseqs_with_ancestor, v_obs, CNV_information_test, params.freq_cutoff)
                             clone_frequency_cnv2.regress_cnv()	
                             decomseqs_with_ancestor_clone_freq=clone_frequency_cnv2.Tumor2Clone_frequency		
                             decomseqs_with_ancestor=	clone_frequency_cnv2.hitclone_seq_builder	
                             seqs_with_ancestor, A1, A2  = OutFile.ReNameCloFreMeg(decomseqs_with_ancestor, decomseqs_with_ancestor_clone_freq, 'number')
                             ItC+=1	
                             
							 
                             print ('clones were decomposed, so repeat the decomposition')									 
                       else:
                             print ('decomposed clones do not have evolutionary structure.. So, we will not decompose them.')  
                             summary_file+=	'decomposed clone genotypes are not good (no evolutionary structure), so clones were not decomposed.\n'								 					 
                             Repeat='n'
                             final_seq = clone_frequency_cnv.hitclone_seq_builder
                             final_clofre = clone_frequency_cnv.Tumor2Clone_frequency					 				 
					                  							
        print ('clone decomposition is complete!')
        print ('assign lost SNVs')
        FinalCloLs,FinalClo2Seq=Align.name2seq(final_seq)		
        LostSNVpos=Align.GetSharePosi1(FinalClo2Seq,'A')
      #  print (LostSNVpos)
        SingleTuLostSNVpos=Align.FilterLostSNV(LostSNVpos,v_obs)	
        AddCloSeq=Align.AssignLostSNV(SingleTuLostSNVpos,v_obs,final_clofre,params.freq_cutoff,final_seq)	
        if AddCloSeq!=[]:		
         print ('lost SNV was assigned so compute clone frequency')		
         OriCloLs,OriClo2Seq=Align.name2seq(final_seq)
         CloID=1	
         for Aclo in AddCloSeq:
             OriClo2Seq['#A'+str(CloID)]=Aclo	
             CloID+=1	
			 
         OriClo2Seq['#hg19']='A'*num_sites
					   
					   
					   
         up_seq_build = Align.UpMeg(OriClo2Seq,list(OriClo2Seq.keys()))		
	
         clone_frequency_cnv2 = CloneFrequencyComputer_cnv1(up_seq_build, v_obs, CNV_information_test, params.freq_cutoff)
         clone_frequency_cnv2.regress_cnv()	
         final_clofre=clone_frequency_cnv2.Tumor2Clone_frequency		
         final_seq=	clone_frequency_cnv2.hitclone_seq_builder
      #  final_seq1, final_clone_frequency1, final_clone_order1  = OutFile.ReNameCloFreMeg(final_seq,final_clofre,   'number') ###	

     #   print ('test clone hit and remove insignificant clones')
     #   significant_clone=cluster_test()		
     #   final_seq, final_clofre= significant_clone.remove_insignificant_clones_add(v_obs, final_clone_frequency1,  final_seq1, CNV_information_test, Significant_cutoff)
        print ('making output files')	

        final_seq10, final_clone_frequency1, final_clone_order1  = OutFile.ReNameCloFreMeg(final_seq,final_clofre,   'number') 		
        final_seq1=[]
        for i in final_seq10:
           	final_seq1.append(i.replace('?','A'))	
        Align.save_mega_alignment_to_file(params.input_id + '_CloneFinderPlus.meg', final_seq1)
        CloFreAna.save_frequency_table_to_file(params.input_id + '_CloneFinderPlus.txt',  final_clone_frequency1,  [])	

        id = 'final'		
        status = tree_builder.do_mega_mp(final_seq1, id)
        if status == True:
                A1, tree = tree_builder.alignment_least_back_parallel_muts() 
                Rooted=AnalyzeTree.RootTree(tree)
                InferAncestor = MegaAncestor()
                InferAncestor.alignment_file = final_seq1
                InferAncestor.input_tree_file = Rooted
     			
                ancestor_states, offspring2ancestor, cell2code, code2cell = InferAncestor.retrieve_ancestor_states()
                RescaledTree=InferAncestor.get_scaledNWK()					
                OutFile.GetOut(params.input_id+ '_CloneFinderPlus.nwk', RescaledTree.replace('hg19:','Normal:'))				                     
	
		

		
	
#######################
os.remove(params.input_id + '.txt')
os.remove(params.input_id + '-CNV.txt')

#timeFile = params.input_id + '_summary.txt'
#endTime = datetime.datetime.now()
#print (endTime)
#totalTime = (endTime - startTime)
#print (totalTime)
#summary_file += 'Run time: ' + str(totalTime) + '\n'
#OutTime=open(timeFile,'w')
#OutTime.write(summary_file)
#OutTime.close()
#######################			
