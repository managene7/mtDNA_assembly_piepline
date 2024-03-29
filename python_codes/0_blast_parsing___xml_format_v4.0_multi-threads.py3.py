#Parsing of XML file that contains BLAST results
def writeAlignSeq(q_name, q_len, q_from, q_to, h_name, h_len, h_from, h_to, score, evalue, idty, q_seq, midline, h_seq, align_len, seq_out, align_out):
    q_len=str(q_len)
    q_from=str(q_from)
    q_to=str(q_to)
    h_len=str(h_len)
    h_from=str(h_from)
    h_to=str(h_to)
    score=str(score)
    evalue=str(evalue)
    idty=str(idty)
    
    q_seq_name=">"+q_name+"___q"+q_len+"_"+q_from+"-"+q_to
    h_seq_name=">"+h_name+"___h"+h_len+"_"+h_from+"-"+h_to
    # write sequence
    seq_out.write(q_seq_name+"\n")
    seq_out.write(q_seq.replace('-','')+"\n")
    seq_out.write(h_seq_name+"\n")
    seq_out.write(h_seq.replace('-','')+"\n")

    align_out.write("<div style='width:1090px'>\n")
    align_out.write("<p><br><br></p>\n")
    # write align
    q_name_write="<p style='font-size: medium; font-weight: bold;'>Query-HSP name: %s</p> " % q_name
    q_len_write="<p style='font-size: medium; font-weight: bold;'>Query-HSP length: %s </p>\n" % q_len   
    q_posi_write="<p style='font-size: medium; font-weight: bold;'>From: %s To: %s</p>\n" % (q_from, q_to)    
    align_out.write(q_name_write)
    align_out.write(q_len_write)
    align_out.write(q_posi_write)
    align_out.write("<p><br></p>\n")

    h_name_write="<p style='font-size: medium; font-weight: bold;'>Hit-HSP name: %s </p>\n" % h_name  
    h_len_write="<p style='font-size: medium; font-weight: bold;'>Hit-HSP length: %s</p>\n" % h_len 
    h_posi_write="<p style='font-size: medium; font-weight: bold;'>From: %s To: %s</p>\n" % (h_from, h_to)  
    align_out.write(h_name_write)
    align_out.write(h_len_write)
    align_out.write(h_posi_write)

    align_info="<p style='font-size: medium; font-weight: bold;'>Score: %s Evalue: %s Identity: %s</p>\n" % (score, evalue,str(idty))
    align_out.write(align_info)
    align_out.write("<p><br></p>\n")

    align_out.write("<p>Alignment</p>\n")
    show_align_len=int(align_len)
    show_unit_num=int(show_align_len/50)
    show_unit_rest=int(show_align_len%50)


    if show_unit_num>=1:
        init=0
        scale_list=[]
        for i in range(show_unit_num):
            if init==0:
                scale_list.append(str(1)+"-"*48)
                init=2
            elif init==2:
                scale_list.append(str(i*50)+"-"*48)
                init=3
            else:
                scale_list.append(str(i*50)+"-"*47)
        if show_unit_rest != 0:
            if 50*show_unit_num+show_unit_rest >= 100:
                scale_list.append(str(50*show_unit_num)+"-"*(show_unit_rest-3)+str(50*show_unit_num+show_unit_rest))
            else:
                scale_list.append(str(50*show_unit_num)+"-"*(show_unit_rest-2)+str(50*show_unit_num+show_unit_rest))
        else:
            scale_list.append(str(50*show_unit_num))


        scale="".join(scale_list)
    else:
        scale="1"+"-"*(show_unit_rest-2)+str("show_unit_rest")
    
    

    align_out.write("<p>%s</p>\n" % scale)

    aligned_q_seq=q_seq
    aligned_midline=midline.replace(" ", "-")
    aligned_h_seq=h_seq
    
    repeat=int(len(aligned_q_seq)/show_align_len)+1
    for i in range(repeat):
        l_posi=i*show_align_len
        r_posi=(i+1)*show_align_len

        q_seq=aligned_q_seq[l_posi:r_posi]
        align_out.write("<p>%s</p>\n" % q_seq)
        
        midline=aligned_midline[l_posi:r_posi]
        align_out.write("<p>%s</p>\n" % midline)
        
        h_seq=aligned_h_seq[l_posi:r_posi]
        align_out.write("<p>%s</p>\n" % h_seq)
        align_out.write("<p><br></p>\n")
    align_out.write("</div>")

    
#________________ option parse _______________________________
import sys 

args = sys.argv[1:]

option_dict={'-lim_evalue':10, '-ov_len':1, '-seq_parse':'0', '-lim_identity':0.6,'-hsp_identity':0.6,'-lim_score':50,'-hit_function':'1','-num_hsp':10000000,'-num_hit':10000000,'-out':"",'-num_threads':'8',"-mem_per_thread":"auto"}
help_description="""
________________________________________________________________________________


Usage;

If default value exists, corresponding value can be omitted.

-xml            file name of blast result (xml format)
-num_threads    number of threads (default is 8)
-mem_per_thread RAM size per threads (ex: 20000000 (bytes) default: auto)
-out            output file name without extension that you want to (ex.. test)
-seq_parse      write the aligned sequence in separated files.
                int number (ex: 150): write aligned sequence and set the align length to be shown in a line 
                0: don't write. (Default is 0)
                150 or 200 or 250: write
-lim_evalue     evalue threshold of BLAST hit. Default is e-10
-lim_identity   identity threshold of BLAST hit. Default is 0.60
-hsp_identity   identity threshold of BLAST hsp. Default is 0.60
-num_hit        number of hit to be parsed
-num_hsp        number of hsp to be parsed
-lim_score      score threshold of hsps in a hit. Default is 50.
-hit_function   to filter out 'unknown' or 'hypothetical' in hit description.
                1: not filtering 2: filtering. Default is 1
-ov_len         length threshold of result overlap result. Default is 1.

________________________________________________________________________________
"""
if args==[]:
    print (help_description)
    quit()
    
for i in range(len(args)):
    if args[i].startswith("-"):
        try:
            option_dict[args[i]]=args[i+1]
        except:
            if args[0]=="-help":
                print(help_description)
                quit()
#        else:
#            print "There maybe an option error!"
#            quit()
import csv

def blast_xml_parsing(job_id, return_dict, xml, xml_sub_cont, num_threads, seq_parse, lim_evalue, lim_identity, hsp_identity, num_hit, num_hsp, lim_score, hit_function, ov_len):
    option_dict['-xml']=xml
    option_dict['-seq_parse']=seq_parse
    option_dict['-lim_evalue']=lim_evalue
    option_dict['-lim_identity']=lim_identity
    option_dict['-hsp_identity']=hsp_identity
    option_dict['-num_hit']=num_hit
    option_dict['-num_hsp']=num_hsp
    option_dict['-lim_score']=lim_score
    option_dict['-hit_function']=hit_function
    option_dict['-ov_len']=ov_len
    option_dict['-num_threads']=num_threads

    #__file open to write aligned sequence__
    if option_dict['-seq_parse']!="0":
        in_name=option_dict['-xml']
        seq_out=open(in_name[:-4]+".fasta",'w')
        align_out=open(in_name[:-4]+".html",'w')
        align_out.write("<div style='display: container; width:auto; font-family: consolas; font-size: small; line-height: 0.3em'>\n")
    #____________________________________________


    input_file=open(option_dict['-xml'], 'r')

    #print "\n\nBLAST xml parsing in progression......" ##################################################################

    
    csv.field_size_limit(100000000)
    import shelve
    
    csv_file_cont=[]

    query_def="temporary_query_def"
    #This is MODULE 1 to read parameter of a blastp result by dictionary form
    #structure {query_def:[query_len,[hit_num, hit_id, hit_len, total_score, E-value, mean_identity, query_coverage, subject_coverage]
    #Pickle out_file will be returned


    #________________________________________query circuit_____________________
    parameter_dic={}
    while 1:
        if xml_sub_cont==[]:
            break
        line=xml_sub_cont[0].strip();del xml_sub_cont[0]

        if not line:
            #print "\n\n\nBLAST xml parsing completed. Open the csv file by MicroSoft EXCEL..." ############################3
            break 

        
        if line.startswith("<Iteration>"):
            while 1:
                if xml_sub_cont==[]:
                    break
                line=xml_sub_cont[0].strip();del xml_sub_cont[0]
                if line.startswith("<Iteration_query-def>"):
                    query_def=line[21:-22]
                if line.startswith("<Iteration_query-len>"):
                    query_len=line[21:-22]
                    parameter_dic[query_def]=[int(query_len)]
                    parameter_dic[query_def].append({})
                if line.startswith("<Iteration_hits></Iteration_hits>"):
                    csv_file_cont.append([query_def,parameter_dic[query_def][0], "-", "-", "-", "-", "-", "-", "-","-", "-"]) 
                    break
                if line.startswith("</Iteration>"):
                    
                    
    #__________________garbage hit deletion______________________________



    #This is MODULE 2 to analysis
        #Select maximum score hit
            #csv export

                    if query_def != "temporary_query_def":


                #___________________________________hit deletion containing nameless___________________  

                        if str(option_dict['-hit_function'])=="2":
                            hit_id_list3=list(parameter_dic[query_def][1].keys())
            
                            hit_e_value_list=[]
                            for h_id in hit_id_list3:
                                hit_e_value_list.append([parameter_dic[query_def][1][h_id]['e-value'], (1000/parameter_dic[query_def][1][h_id]['total_score']), parameter_dic[query_def][1][h_id]['mean_identity'], h_id, parameter_dic[query_def][1][h_id]['hit_range_merged'], parameter_dic[query_def][1][h_id]['hit_def']])######### To delete Vitis vinifera hit, hit_def added
                            hit_e_value_list.sort()
            
            
                            hit_garbage_list=[]
                        
                            for key in hit_e_value_list:
                            
                                if 'hypothetical protein' in key[5]: ######### To delete nameless protein hit######
                                    hit_garbage_list.append(key[3])######### To delete nameless protein hit######
                                if 'predicted protein' in key[5]: ######### To delete nameless protein hit######
                                    hit_garbage_list.append(key[3])######### To delete nameless protein hit######
                                if 'uncharacterized' in key[5]:######### To delete nameless protein hit######
                                    hit_garbage_list.append(key[3])######### To delete nameless protein hit######
                                if 'uncharacterized protein' in key[5]:######### To delete nameless protein hit######
                                    hit_garbage_list.append(key[3])######### To delete nameless protein hit######
                                    
                                else:
                                    break
            
                            if len(hit_garbage_list) >= 1:
                                for key in hit_garbage_list:
                                    del parameter_dic[query_def][1][key]  #hit grabage deletion by e-value

                            hit_garbage_list=[] #To clear memory

                
                        #__________________garbage hit deletion by overlap______________________________
        
                        hit_id_list=list(parameter_dic[query_def][1].keys())
                        if len(hit_id_list) > 1:
                            hit_e_value_list=[]
                            for h_id in hit_id_list:
                                hit_e_value_list.append([                                
                                    (1000/parameter_dic[query_def][1][h_id]['max_score']),
                                    (100/parameter_dic[query_def][1][h_id]['identity_of_max_score']),
                                    parameter_dic[query_def][1][h_id]['e-value'],
                                    h_id,
                                    parameter_dic[query_def][1][h_id]['hit_range_merged']])
                            hit_e_value_list.sort()

                            n=0
                            hit_garbage_list=[]
                
                            for key in hit_e_value_list:
                                n=n+1
                                if n>int(option_dict['-num_hit']):
                                    hit_garbage_list.append(key[3])
                            
                                if n == len(hit_e_value_list):
                                    break
                                for key_2 in hit_e_value_list[n:]:
                                    if key[3] != key_2[3]:  # to avoid muplicated hit deletion
                                        if len(key_2[4].intersection(key[4])) > int(option_dict["-ov_len"]):
                                            hit_garbage_list.append(key_2[3])
                                    else: pass


                            temp=set(hit_garbage_list)
                            hit_garbage_new=list(temp)

                            if len(hit_garbage_new) >= 1:
                                for key in hit_garbage_new:
                                    del parameter_dic[query_def][1][key]  #hit grabage deletion by range comparison

                            hit_garbage_new=[] #To clear memory
                            temp=[] #To clear memory
                            hit_e_value_list=[] #To clear memory


                #___________________________________hit grabage deletion by e-value and identity threshold___________________  

                        hit_id_list2=list(parameter_dic[query_def][1].keys())
                
            
                        hit_e_value_list=[]
                        for h_id in hit_id_list2:
                            hit_e_value_list.append([
                                parameter_dic[query_def][1][h_id]['e-value'],
                                parameter_dic[query_def][1][h_id]['max_score'],
                                parameter_dic[query_def][1][h_id]['identity_of_max_score'],
                                h_id,
                                parameter_dic[query_def][1][h_id]['hit_range_merged']])
                        hit_e_value_list.sort()
            
            
                        hit_garbage_list=[]
                        set_e_value='1e-%s' % str(option_dict["-lim_evalue"])
                        for key in hit_e_value_list:
                            if float(key[0]) > float(set_e_value):
                                hit_garbage_list.append(key[3])
                            if float(key[2]) < float(option_dict['-lim_identity']):
                                hit_garbage_list.append(key[3])
                            if float(key[1]) <float(option_dict['-lim_score']):
                                hit_garbage_list.append(key[3])
                                
                            
                        new_hit_garbage_list=list(set(hit_garbage_list))
                        if len(hit_garbage_list) >= 1:
                            for key in new_hit_garbage_list:
                                del parameter_dic[query_def][1][key]  #hit grabage deletion by e-value
                        new_hit_garbage_list=[] # to clear memory
                        hit_garbage_list=[] # to clear memory
                        hit_e_value_list=[] # to clear memory


    #_________________________csv, sequence_dic writing_________________________________________________    


                    #parameter_dic={query_def:'[query_len,{hit_id: {'Hsp_all':hsp_dic,hit_range_merged, hit_def, hit_accession, hit_length,total_score, e-value, mean_identity, query_coverage, subject_coverage}]}'}'}


                        No_result_num=0
        
                        hit_id=list(parameter_dic[query_def][1].keys())
                        evalue_hit_list=[]
                        if len(list(parameter_dic[query_def][1].keys()))==0:
                            No_result_num=No_result_num+1
                            csv_file_cont.append([query_def,parameter_dic[query_def][0], "-", "-", "-", "-", "-", "-", "-","-", "-"])

                        temp_dic={}
                        #seq_file=shelve.open(option_dict["-out"]+"_sequence_db.dic") # DB construction for matched query region 
                            
                        if len(list(parameter_dic[query_def][1].keys()))==1:
                            for sub_name in hit_id:
                                evalue_hit_list.append([parameter_dic[query_def][1][sub_name]['e-value'],sub_name])
                            min_hit_list=min(evalue_hit_list)
                            min_hit_id=min_hit_list[1]


                            
                            if len(parameter_dic[query_def][1][min_hit_id]['Hsp_all'])==1:

                                sub_key=list(parameter_dic[query_def][1][min_hit_id]['Hsp_all'].keys())[0]

                                csv_value_list=[
                                    query_def,
                                    parameter_dic[query_def][0],
                                    parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_query-from'],
                                    parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_query-to'],
                                    parameter_dic[query_def][1][min_hit_id]['hit_def'],
                                    parameter_dic[query_def][1][min_hit_id]['hit_length'],
                                    parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_hit-from'],
                                    parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_hit-to'],
                                    parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_bit-score'],
                                    parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_evalue'],
                                    round(float(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_identity'])/float(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_align-len']),4)
                                    ]
                                if int(option_dict['-seq_parse'])!=0:
                                    csv_value_list.append(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_qseq'])
                                    csv_value_list.append(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_midline'])
                                    csv_value_list.append(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][sub_key]['Hsp_hseq'])
                                    
                                    csv_value_list=csv_value_list+[option_dict['-seq_parse'], seq_out, align_out]
                                    writeAlignSeq(csv_value_list[0], csv_value_list[1], csv_value_list[2], csv_value_list[3], csv_value_list[4], csv_value_list[5], csv_value_list[6], csv_value_list[7], 
                                                csv_value_list[8], csv_value_list[9], csv_value_list[10], csv_value_list[11], csv_value_list[12], csv_value_list[13], 
                                                csv_value_list[14], csv_value_list[15], csv_value_list[16])
                                    csv_value_list=csv_value_list[0:-6]
                                csv_file_cont.append(csv_value_list)



                            if len(parameter_dic[query_def][1][min_hit_id]['Hsp_all']) > 1:


                                Hsp_id_temp=[]    # To sort hsps by position-start
                                
                                Hsp_keys=list(parameter_dic[query_def][1][min_hit_id]['Hsp_all'].keys())
                                
                                for key in Hsp_keys:
                                    Hsp_id_temp.append((parameter_dic[query_def][1][min_hit_id]['Hsp_all'][key]['Hsp_query-from'],key))

                                Hsp_id_temp.sort()

                                Hsp_id=[]    
                                for id in Hsp_id_temp:   
                                    Hsp_id.append(id[1])  # To sort hsps by position-end

                                Hsp_id_temp=[]    #To clear memory
                                Hsp_keys=[]       #To clear memory
                                
                                m=0
                                for hsp in Hsp_id:
                                    m=m+1
                                    csv_value_list=[
                                        query_def,
                                        parameter_dic[query_def][0],
                                        parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_query-from'],
                                        parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_query-to'],
                                        parameter_dic[query_def][1][min_hit_id]['hit_def']+'_continued-'+str(m),
                                        parameter_dic[query_def][1][min_hit_id]['hit_length'],
                                        parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_hit-from'],
                                        parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_hit-to'],
                                        parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_bit-score'],
                                        parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_evalue'],
                                        round(float(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_identity'])/float(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_align-len']),4)
                                        ]

                                    if int(option_dict['-seq_parse'])!=0:
                                        csv_value_list.append(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_qseq'])
                                        csv_value_list.append(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_midline'])
                                        csv_value_list.append(parameter_dic[query_def][1][min_hit_id]['Hsp_all'][hsp]['Hsp_hseq'])

                                        csv_value_list=csv_value_list+[option_dict['-seq_parse'], seq_out, align_out]
                                        writeAlignSeq(csv_value_list[0], csv_value_list[1], csv_value_list[2], csv_value_list[3], csv_value_list[4], csv_value_list[5], csv_value_list[6], csv_value_list[7], 
                                                    csv_value_list[8], csv_value_list[9], csv_value_list[10], csv_value_list[11], csv_value_list[12], csv_value_list[13], 
                                                    csv_value_list[14], csv_value_list[15], csv_value_list[16])
                                        csv_value_list=csv_value_list[0:-6]

                                    csv_file_cont.append(csv_value_list)

                                
                            

                            


                            #temp_dic[">"+query_def]=parameter_dic[query_def][1][min_hit_id]['Hsp_qseq']# DB writing for matched query region (left)
                            
                            #seq_name_list.write(">"+query_def) # Matched sequence name writing
                            #seq_name_list.write("\n")
                            
                            #temp_dic[">"+parameter_dic[query_def][1][min_hit_id]['hit_def']]=parameter_dic[query_def][1][min_hit_id]['Hsp_hseq'] # DB writing for matched query region (right)
                            
                            #seq_name_list.write(">"+parameter_dic[query_def][1][min_hit_id]['hit_def']) # Matched sequence name writing
                            #seq_name_list.write("\n")

                            #seq_file[">"+query_def]=temp_dic # DB writing for matched query region by dictionary
                            



                        if len(list(parameter_dic[query_def][1].keys())) >= 2:
                            for sub_name in hit_id:
                                evalue_hit_list.append([
                                    parameter_dic[query_def][1][sub_name]['e-value'],(1000/parameter_dic[query_def][1][sub_name]['total_score']),parameter_dic[query_def][1][sub_name]['mean_identity'],sub_name])
                            evalue_hit_list.sort()


                            if len(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'])==1:

                                sub_key=list(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'].keys())[0]
                                
                                csv_value_list=[
                                    query_def+'_continued-1',
                                    parameter_dic[query_def][0],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_query-from'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_query-to'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['hit_def'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['hit_length'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_hit-from'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_hit-to'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_bit-score'],
                                    parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_evalue'],
                                    round(float(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_identity'])/float(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_align-len']),4)
                                    ]
                                if int(option_dict['-seq_parse'])!=0:
                                    csv_value_list.append(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_qseq'])
                                    csv_value_list.append(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_midline'])
                                    csv_value_list.append(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][sub_key]['Hsp_hseq'])

                                    csv_value_list=csv_value_list+[option_dict['-seq_parse'], seq_out, align_out]
                                    writeAlignSeq(csv_value_list[0], csv_value_list[1], csv_value_list[2], csv_value_list[3], csv_value_list[4], csv_value_list[5], csv_value_list[6], csv_value_list[7], 
                                                csv_value_list[8], csv_value_list[9], csv_value_list[10], csv_value_list[11], csv_value_list[12], csv_value_list[13], 
                                                csv_value_list[14], csv_value_list[15], csv_value_list[16])
                                    csv_value_list=csv_value_list[0:-6]

                                csv_file_cont.append(csv_value_list)



                            if len(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all']) > 1:

                                Hsp_id_temp=[]    # To sort hsps by position-start
                                Hsp_keys=list(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'].keys())    
                                for key in Hsp_keys:
                                    Hsp_id_temp.append((parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][key]['Hsp_query-from'],key))

                                Hsp_id_temp.sort()
                                                                                                                        
                                Hsp_id=[]    
                                for id in Hsp_id_temp:   
                                    Hsp_id.append(id[1])  # To sort hsps by position-end

                                Hsp_id_temp=[] # To clear memory
                                Hsp_keys=[]    # To clear memory

                                m=0
                                for hsp in Hsp_id:
                                    m=m+1
                                    csv_value_list=[
                                        query_def+'_continued-1',
                                        parameter_dic[query_def][0],
                                        parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_query-from'],
                                        parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_query-to'],
                                        parameter_dic[query_def][1][evalue_hit_list[0][3]]['hit_def']+'_continued-'+str(m),
                                        parameter_dic[query_def][1][evalue_hit_list[0][3]]['hit_length'],
                                        parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_hit-from'],
                                        parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_hit-to'],
                                        parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_bit-score'],
                                        parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_evalue'],
                                        round(float(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_identity'])/float(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_align-len']),4)
                                        ]

                                    if int(option_dict['-seq_parse'])!=0:
                                        csv_value_list.append(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_qseq'])
                                        csv_value_list.append(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_midline'])
                                        csv_value_list.append(parameter_dic[query_def][1][evalue_hit_list[0][3]]['Hsp_all'][hsp]['Hsp_hseq'])

                                        csv_value_list=csv_value_list+[option_dict['-seq_parse'], seq_out, align_out]
                                        writeAlignSeq(csv_value_list[0], csv_value_list[1], csv_value_list[2], csv_value_list[3], csv_value_list[4], csv_value_list[5], csv_value_list[6], csv_value_list[7], 
                                                    csv_value_list[8], csv_value_list[9], csv_value_list[10], csv_value_list[11], csv_value_list[12], csv_value_list[13], 
                                                    csv_value_list[14], csv_value_list[15], csv_value_list[16])
                                        csv_value_list=csv_value_list[0:-6]

                                    csv_file_cont.append(csv_value_list)


                            l=0
                            for key in evalue_hit_list[1:]:
                                l=l+1


                                if len(parameter_dic[query_def][1][key[3]]['Hsp_all'])==1:

                                    sub_key=list(parameter_dic[query_def][1][key[3]]['Hsp_all'].keys())[0]
                                
                                    csv_value_list=[
                                        '%s_continued-%s' % (query_def,str(l+1)),
                                        parameter_dic[query_def][0],
                                        parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_query-from'],
                                        parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_query-to'],
                                        parameter_dic[query_def][1][key[3]]['hit_def'],
                                        parameter_dic[query_def][1][key[3]]['hit_length'],
                                        parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_hit-from'],
                                        parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_hit-to'],
                                        parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_bit-score'],
                                        parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_evalue'],
                                        round(float(parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_identity'])/float(parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_align-len']),4)
                                        ]
                                    if int(option_dict['-seq_parse'])!=0:
                                        csv_value_list.append(parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_qseq'])
                                        csv_value_list.append(parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_midline'])
                                        csv_value_list.append(parameter_dic[query_def][1][key[3]]['Hsp_all'][sub_key]['Hsp_hseq'])

                                        csv_value_list=csv_value_list+[option_dict['-seq_parse'], seq_out, align_out]
                                        writeAlignSeq(csv_value_list[0], csv_value_list[1], csv_value_list[2], csv_value_list[3], csv_value_list[4], csv_value_list[5], csv_value_list[6], csv_value_list[7], 
                                                    csv_value_list[8], csv_value_list[9], csv_value_list[10], csv_value_list[11], csv_value_list[12], csv_value_list[13], 
                                                    csv_value_list[14], csv_value_list[15], csv_value_list[16])
                                        csv_value_list=csv_value_list[0:-6]

                                    csv_file_cont.append(csv_value_list)



                                if len(parameter_dic[query_def][1][key[3]]['Hsp_all']) > 1:

                                    Hsp_id_temp=[]    # To sort hsps by position-start
                                    Hsp_keys=list(parameter_dic[query_def][1][key[3]]['Hsp_all'].keys())
                                    for key0 in Hsp_keys:
                                        Hsp_id_temp.append((parameter_dic[query_def][1][key[3]]['Hsp_all'][key0]['Hsp_query-from'],key0))
                                    Hsp_id_temp.sort()

                                    Hsp_id=[]    
                                    for id in Hsp_id_temp:   
                                        Hsp_id.append(id[1])  # To sort hsps by position-end

                                    Hsp_id_temp=[] #To clear memory
                                    Hsp_keys=[]    #To clear memory

                                    
                                    m=0
                                    for hsp in Hsp_id:
                                        m=m+1
                                        csv_value_list=[
                                            '%s_continued-%s' % (query_def,str(l+1)),
                                            parameter_dic[query_def][0],
                                            parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_query-from'],
                                            parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_query-to'],
                                            parameter_dic[query_def][1][key[3]]['hit_def']+'_continued-'+str(m),
                                            parameter_dic[query_def][1][key[3]]['hit_length'],
                                            parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_hit-from'],
                                            parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_hit-to'],
                                            parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_bit-score'],
                                            parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_evalue'],
                                            round(float(parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_identity'])/float(parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_align-len']),4)
                                            ]

                                        if int(option_dict['-seq_parse'])!=0:
                                            csv_value_list.append(parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_qseq'])
                                            csv_value_list.append(parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_midline'])
                                            csv_value_list.append(parameter_dic[query_def][1][key[3]]['Hsp_all'][hsp]['Hsp_hseq'])

                                            csv_value_list=csv_value_list+[option_dict['-seq_parse'], seq_out, align_out]
                                            writeAlignSeq(csv_value_list[0], csv_value_list[1], csv_value_list[2], csv_value_list[3], csv_value_list[4], csv_value_list[5], csv_value_list[6], csv_value_list[7], 
                                                        csv_value_list[8], csv_value_list[9], csv_value_list[10], csv_value_list[11], csv_value_list[12], csv_value_list[13], 
                                                        csv_value_list[14], csv_value_list[15], csv_value_list[16])
                                            csv_value_list=csv_value_list[0:-6]

                                        csv_file_cont.append(csv_value_list)


                    parameter_dic={} # To clear memory
                    hsp_dic={} # To clear memory
                    hit_parameter_dic={} # To clear memory
                    hsp_sub_dic={} # To clear memory
                    break









    #____________________________Hit Circuit____________________________________
                
                if line.startswith("<Hit>"):

                    
                    hsp_dic={}
                    
                    while 1:
                        if xml_sub_cont==[]:
                            break
                    
                        line=xml_sub_cont[0].strip();del xml_sub_cont[0]
                                                
                        if line.startswith("<Hit_id>"):
                            hit_id=line[8:-9]
                        if line.startswith("<Hit_def>"):
                            hit_def=line[9:-10]
                        if line.startswith("<Hit_accession>"):
                            hit_accession=line[15:-16]
                            
                        if line.startswith("<Hit_len>"):
                            hit_length=int(line[9:-10])#hit_length

                                

    #________________________________Garbage hsp deletion________________________________________________

                                
                        if line.startswith("</Hit_hsps>"):
                            hsp_dic_keys=list(hsp_dic.keys())
                                            
                            

                            if len(hsp_dic_keys) ==1:
                                set_value=list(range(hsp_dic[hsp_dic_keys[0]]['Hsp_query-from'], hsp_dic[hsp_dic_keys[0]]['Hsp_query-to']+1)) 
                                hsp_hit_merge=set(list(set_value))
                                
                            
                            if len(hsp_dic_keys) > 1:
                                for k in hsp_dic_keys:
                                    hsp_dic[k]['hit_range'] = set(list(range(hsp_dic[k]['Hsp_query-from'], hsp_dic[k]['Hsp_query-to']+1))) # 'hit_range' entered to hsp value
                                    
                                hsp_score_sort=[]
                                for key in hsp_dic_keys:
                                    hsp_score_sort.append([hsp_dic[key]['Hsp_bit-score'], key]) #To sort by score
                                
                                hsp_score_sort.sort()
                                hsp_score_sort.reverse()                                        #To sort by score

                                n=0                                                         #To sort by score
                                hsp_garbage_list=[]                                         #To sort by score
                                for key in hsp_score_sort:                                  #To sort by score
                                    n=n+1
                                    if n == len(hsp_score_sort):                                #To sort by score
                                        break
                                    for key_2 in hsp_score_sort[n:]:
                                        if len(hsp_dic[key[1]]['hit_range'].intersection(hsp_dic[key_2[1]]['hit_range'])) >= int(option_dict['-ov_len']):
                                            hsp_garbage_list.append(key_2[1])

                                temp1=set(hsp_garbage_list)
                                hsp_garbage_new=list(temp1)
                                
                                            
                                if len(hsp_garbage_new) >= 1:
                                    for key in hsp_garbage_new:
                                        del hsp_dic[key]  #Hsp garbage deletion by overlap

                                hsp_garbage_list=[] #To clear memory
                                hsp_garbage_new=[] #To clear memory
                                temp1=[] #To clear memory
                                hsp_score_sort=[] #To clear memory

                                ##___________garbage hsp deletion by score_start______________________

                                hsp_dic_keys2=list(hsp_dic.keys())
                
                                if len(hsp_dic_keys2)>1:
                                    hsp_score_sort=[]
                                    for key in hsp_dic_keys2:
                                        hsp_score_sort.append([hsp_dic[key]['Hsp_bit-score'], float(hsp_dic[key]['Hsp_identity'])/float(hsp_dic[key]['Hsp_align-len']), key])
                                    hsp_score_sort.sort()
                                    hsp_score_sort.reverse()
                                        
                                    hsp_sort_list=hsp_score_sort[1:]
                                    hsp_garbage_list=[]
                                    hsp_num=1
                                    for key in hsp_sort_list:
                                        hsp_num=hsp_num+1
                                        if hsp_num > int(option_dict['-num_hsp']):
                                            hsp_garbage_list.append(key[2])
                                        else:
                                            if float(key[0]) < float(option_dict['-lim_score']):
                                                hsp_garbage_list.append(key[2])
        
                                            if float(key[1]) < float(option_dict['-hsp_identity']):##############################################################
                                                hsp_garbage_list.append(key[2])

                                    hsp_garbage_list=list(set(hsp_garbage_list))
                                        
                                    if len(hsp_garbage_list) >= 1:
                                        for key in hsp_garbage_list:
                                            del hsp_dic[key]
                                    hsp_score_sort=[] #To clear memory
                                    hsp_garbage_list=[] #To clear memory
                                    hsp_sort_list=[] #To clear memory
                                ##___________garbage hsp deletion by score_end______________________
                                        



                                new_key=list(hsp_dic.keys())

                                if len(new_key) > 1:
                                    hsp_hit_merge=set([])
                                    for key in new_key:
                                        hsp_hit_merge=set(hsp_dic[key]['hit_range']).union(hsp_hit_merge)
                                    
                                if len(new_key) == 1:
                                    hsp_hit_merge=hsp_dic[new_key[0]]['hit_range']   #Hsp range merge
                                    
                                    
                                if len(new_key) == 0:
                                    hsp_hit_merge=set([])
                                    
                            
                            
                            
                                #hsp_dic={'1':
                                    #{'hit_range': "234,2345,34,5346,34,6...",
                                    #'Hsp_bit-score': 234235,
                                    #'Hsp_query-to':8140,
                                    #'Hsp_align-len':3019,
                                    #'Hsp_identity':3019,
                                    #'Hsp_query-from':5122,
                                    #'Hsp_evalue':0.0, 
                                    #'Hsp_hit-from':9212,
                                    #'Hsp_hseq':'atatgc...',
                                    #'Hsp_midline':'||||...',
                                    #'Hsp_qseq':'atgc...',
                                    #'Hsp_hit-to':12230}}

    #____________________________________Value calculate and add___________________________________                        

                            
                            
                            hit_parameter_dic={}
                            hit_parameter_dic['hit_range_merged']=hsp_hit_merge  #merged hit range
                            hit_parameter_dic['hit_def']=hit_def #hit_def value
                            hit_parameter_dic['hit_accession']=hit_accession #hit_accession value
                            hit_parameter_dic['hit_length']=hit_length #hit_length value


                            hsp_keys=list(hsp_dic.keys())
                            
                            for k in hsp_keys:
                                hit_parameter_dic['query_from']=hsp_dic[k]['Hsp_query-from']#########
                                hit_parameter_dic['query_to']=hsp_dic[k]['Hsp_query-to']############
                                hit_parameter_dic['hit_from']=hsp_dic[k]['Hsp_hit-from']############
                                hit_parameter_dic['hit_to']=hsp_dic[k]['Hsp_hit-to']############

                            
                            evalue_list=[]
                            for k in hsp_keys:
                                evalue_list.append(hsp_dic[k]['Hsp_evalue'])
                            hit_parameter_dic['e-value']=min(evalue_list)#e-value
                                
                            
                            hsp_identity=0
                            hsp_length=0
                            hsp_identity_list=[]
                            for k in hsp_keys:
                                hsp_identity=hsp_identity + hsp_dic[k]['Hsp_identity']############
                                hsp_length=hsp_length + hsp_dic[k]['Hsp_align-len']############
                                individual_identity=float(hsp_dic[k]['Hsp_identity'])/float(hsp_dic[k]['Hsp_align-len'])
                                hsp_identity_list.append(individual_identity)
                                
                            hit_parameter_dic['mean_identity']=float(hsp_identity) / float(hsp_length) #mean_identity value############
                            hit_parameter_dic['max_identity']=max(hsp_identity_list) #max_identity value############
                            hit_parameter_dic['min_identity']=min(hsp_identity_list) #min_identity value############

                            total_score=0
                            max_score_list=[]
                            for k in hsp_keys:
                                total_score=total_score + hsp_dic[k]['Hsp_bit-score']
                                max_score_list.append([float(hsp_dic[k]['Hsp_bit-score']),float(hsp_dic[k]['Hsp_identity'])/float(hsp_dic[k]['Hsp_align-len'])])
                            hit_parameter_dic['total_score']=total_score #total_score value
                            hit_parameter_dic['max_score']=max(max_score_list)[0]#max_score value
                            hit_parameter_dic['identity_of_max_score']=max(max_score_list)[1]#identity of max_score value


                            red_len_query=0
                            split_length=[]
                            for k in hsp_keys:
                                hit_covered_len_query = hsp_dic[k]['Hsp_query-to'] - hsp_dic[k]['Hsp_query-from'] + 1
                                
                                hit_parameter_dic['query_coverage']=int(hit_covered_len_query)  #query_coverage value

                            hit_covered_len_subject=0
                            for k in hsp_keys:
                                hit_covered_len_subject = hit_covered_len_subject + (hsp_dic[k]['Hsp_hit-to'] - hsp_dic[k]['Hsp_hit-from'] + 1)
                            hit_parameter_dic['subject_coverage'] = abs(int(hit_covered_len_subject)) #subject_coverage value

                            hit_parameter_dic['Hsp_qseq']=hsp_dic[k]['Hsp_qseq'] # Aligned query sequence
                            hit_parameter_dic['Hsp_midline']=hsp_dic[k]['Hsp_midline']
                            hit_parameter_dic['Hsp_hseq']=hsp_dic[k]['Hsp_hseq'] # Aligned hit sequence
                            hit_parameter_dic['Hsp_all']=hsp_dic

        
                            parameter_dic[query_def][1][hit_id]=hit_parameter_dic #{query_def:'[query_len,{hit_id: {'Hsp_all':hsp_dic,hit_range_merged, hit_def, hit_accession, hit_length,total_score, e-value, mean_identity, query_coverage, subject_coverage}]}'}'}


                            break
                    

    #_____________________________________Hsp Circuit__________________________

                        if line.startswith("<Hsp>"):
                            hsp_sub_dic={} #{hit_range, hsp_bit-score,hsp_evalue,hsp_query-from,hsp_query-to,hsp_hit-from,hsp_hit-to,hsp_identity, hsp_align-len}
                            while 1:
                                if xml_sub_cont==[]:
                                    break
                                line=xml_sub_cont[0].strip();del xml_sub_cont[0]
                                if line.startswith("<Hsp_num>"):
                                    hsp_num=line[9:-10]
                                if line.startswith("<Hsp_bit-score>"):
                                    hsp_sub_dic['Hsp_bit-score']=float(line[15:-16])
                                if line.startswith("<Hsp_evalue>"):
                                    hsp_sub_dic['Hsp_evalue']=float(line[12:-13])
                                if line.startswith("<Hsp_query-from>"):
                                    hsp_sub_dic['Hsp_query-from']=int(line[16:-17])
                                if line.startswith("<Hsp_query-to>"):
                                    hsp_sub_dic['Hsp_query-to']=int(line[14:-15])
                                if line.startswith("<Hsp_hit-from>"):
                                    hsp_sub_dic['Hsp_hit-from']=int(line[14:-15])
                                if line.startswith("<Hsp_hit-to>"):
                                    hsp_sub_dic['Hsp_hit-to']=int(line[12:-13])
                                if line.startswith("<Hsp_identity>"):
                                    hsp_sub_dic['Hsp_identity']=int(line[14:-15])
                                if line.startswith("<Hsp_align-len>"):
                                    hsp_sub_dic['Hsp_align-len']=int(line[15:-16])
                                if line.startswith("<Hsp_qseq>"):
                                    Hsp_qseq_list=[]
                                    while 1:
                                        Hsp_qseq_list.append(line)
                                        line=xml_sub_cont[0].strip();del xml_sub_cont[0]
                                        if line.startswith("<Hsp_hseq>"):
                                            Hsp_qseq="".join(Hsp_qseq_list)[10:-11]#.replace("-","") # Conserve the alignment format
                                            hsp_sub_dic['Hsp_qseq']=Hsp_qseq
                                            break
                                if line.startswith("<Hsp_hseq>"):
                                    Hsp_hseq_list=[]
                                    while 1:
                                        Hsp_hseq_list.append(line)
                                        line=xml_sub_cont[0].strip();del xml_sub_cont[0]
                                        if line.startswith("<Hsp_midline>"):
                                            Hsp_hseq="".join(Hsp_hseq_list)[10:-11]#.replace("-","") # Conserve the alignment format
                                            hsp_sub_dic['Hsp_hseq']=Hsp_hseq
                                            break
                                if line.startswith("<Hsp_midline>"):
                                    Hsp_midline_list=[]
                                    while 1:
                                        Hsp_midline_list.append(line)
                                        line=xml_sub_cont[0].strip();del xml_sub_cont[0]
                                        if line.startswith("</Hsp>"):
                                            Hsp_midline="".join(Hsp_midline_list)[13:-14]#.replace("-","") # Conserve the alignment format
                                            hsp_sub_dic['Hsp_midline']=Hsp_midline
                                            break


                                if line.startswith("</Hsp>"):
                                    hsp_dic[hsp_num]=hsp_sub_dic
                                    break
                            

    if option_dict['-seq_parse']!="0":
        align_out.write("</div>")

    if option_dict['-seq_parse']!="0":
        seq_out.close()
        align_out.close()
    return_dict[job_id]=csv_file_cont
    #return csv_file_cont

def main():
    import multiprocessing
    from os.path import getsize
    xml_size=getsize(option_dict['-xml'])
    in_xml=open(option_dict['-xml'],'r')
    init=0
    size=0
    tot_size=0
    manager=multiprocessing.Manager()
    return_dict=manager.dict()
    
    job_id=0
    work_list=[]
    sub_list=[]
    num_threads=int(option_dict['-num_threads'])
    if option_dict['-mem_per_thread']=="auto":
        option_dict['-mem_per_thread']=int(xml_size/(float(option_dict['-num_threads'])))
    print ("xml file size: ", xml_size,"byte")
    print ("Used RAM per thread: ", option_dict['-mem_per_thread'])
    while 1:
        line=in_xml.readline()
        tot_size+=len(line)
        line=line.strip()
        if "<Iteration>" in line:
            init=1
        if init==1:
            if "</Iteration>" in line:
                size+=len(line)
                
                sub_list.append(line)    
                if size >= int(option_dict['-mem_per_thread']):
                    job_id=job_id+1
                    work_list.append((job_id, return_dict, option_dict['-xml'], sub_list, option_dict['-num_threads'], option_dict['-seq_parse'], option_dict['-lim_evalue'], option_dict['-lim_identity'], option_dict['-hsp_identity'], option_dict['-num_hit'], option_dict['-num_hsp'], option_dict['-lim_score'], option_dict['-hit_function'], option_dict['-ov_len']))
                    sub_list=[]
                    size=0
                    if len(work_list)==num_threads:
                        jobs=[]
                        for arg_list in work_list:
                            p=multiprocessing.Process(target=blast_xml_parsing, args=arg_list)
                            jobs.append(p)
                            p.start()
                            
                        for proc in jobs:
                            proc.join()
                        progress=round((tot_size/xml_size)*100, 2)
                        print (str(progress)+"% has processed..\n")
                        work_list=[]
            elif line=="":
                if sub_list!=[]:
                    job_id=job_id+1
                    work_list.append((job_id, return_dict, option_dict['-xml'], sub_list, option_dict['-num_threads'], option_dict['-seq_parse'], option_dict['-lim_evalue'], option_dict['-lim_identity'], option_dict['-hsp_identity'], option_dict['-num_hit'], option_dict['-num_hsp'], option_dict['-lim_score'], option_dict['-hit_function'], option_dict['-ov_len']))
                jobs=[]
                for arg_list in work_list:
                    #print (arg_list)
                    p=multiprocessing.Process(target=blast_xml_parsing, args=arg_list)
                    jobs.append(p)
                    p.start()
                for proc in jobs:
                    proc.join()
                
                progress=100
                print ("Blast xml parsing: "+str(progress)+"% has processed..")
                
                
                out_csv_list=[]
                for id, cont in return_dict.items():
                    if cont!=[]:
                        for sub_cont in cont:
                            out_csv_list.append(sub_cont)
                out_csv_list.sort()
                
                break
            else:
                size+=len(line)
                sub_list.append(line)
            

    if option_dict["-out"]=="":
        option_dict["-out"]=option_dict["-xml"][:-4]+".csv"

    csv_new_open=open(option_dict["-out"],"w", newline="")
    csv_file=csv.writer(csv_new_open)

    title_line=["Query_name","Query_length", "Query_from", "Query_to", "Hit_description", "Hit_length", "Hit_from", "Hit_to", "Score", "e-value", "Identity"]
    csv_file.writerow(title_line)
    for row in out_csv_list:
        csv_file.writerow(row)
    print ("Blast xml parsing has completed!\n")
if __name__=="__main__":
    main()

                                 



