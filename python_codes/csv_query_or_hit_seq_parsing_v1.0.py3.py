#________________ option parse _______________________________
import sys 

args = sys.argv[1:]

option_dict={'-out':"", '-target':'2', '-seq_type':"1", "-method":"1", '-idty':'0.6'}
help_description="""
________________________________________________________________________________

Usage;

If default value exists, corresponding value can be omitted.
-csv        csv file name
-seq        sequence file name
-seq_type   1: fastq (default), 2: fasta
-target     1: query, 2: hit (default is 2)
-method     1: select hits (default), 2: filter out hits
-idty       identity threshold to select (default: 0.6)
-out        Enter output name (option)
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
                print (help_description)
                quit()

def csv_target_parsing(csv_file, target):
    print ("Parsing blast query or hit seq name..\n")
    import csv
    incsv=csv.reader(open(csv_file,'r'))
    out_list=[]
    init=0
    for row in incsv:
        if init==0:
            init=1
            pass
        else:
            if option_dict['-method']=="1":
                if row[4]!="-":
                    if float(row[10])>=float(option_dict['-idty']):
                        if target=="2":
                            out_list.append(row[4].split("_cont")[0])
                        elif target=="1":
                            out_list.append(row[0].split("_cont")[0])
            else:
                if row[4]=="-":
                    if target=="2":
                        out_list.append(row[4].split("_cont")[0])
                    elif target=="1":
                        out_list.append(row[0].split("_cont")[0])
    final=list(set(out_list))
    return final

def fasta_to_dic(fasta):
    print ("Converting fasta to dic..\n\n")
    infile=open(fasta,'r')
    out_dic={}
    init=0
    while 1:
        line=infile.readline().strip()
        if line=="":
            seqs="".join(seq_list)
            out_dic[name]=seqs
            break
        elif ">" in line:
            if init==0:
                name=line[1:]
                seq_list=[]
                init=1
            else:
                seqs="".join(seq_list)
                out_dic[name]=seqs
                name=line[1:]
                seq_list=[]
        else:
            seq_list.append(line)
    return out_dic

def fastq_parsing(fastq, name_list, out):
    print ("Parsing fastq seq..\n\n")
    infile=open(fastq,'r')
    out_file=open(out,'w')
    while 1:
        line1=infile.readline().strip()
        line2=infile.readline().strip()
        line3=infile.readline().strip()
        line4=infile.readline().strip()
        if line1=="":
            break
        elif line1[1:] in name_list:
            out_file.write(line1+"\n")
            out_file.write(line2+"\n")
            out_file.write(line3+"\n")
            out_file.write(line4+"\n")
    return


def main():
    seq_name_list=csv_target_parsing(option_dict['-csv'], option_dict['-target'])
    if option_dict['-seq_type']=="2":
        fasta_dic=fasta_to_dic(option_dict['-seq'])
        print ("Parsing sequences..\n")
        if option_dict['-out']=="":
            if option_dict['-target']=="2":
                option_dict['-out']=option_dict['-csv'][:-4]+"_hit_seq.fasta"
            elif option_dict['-target']=="1":
                option_dict['-out']=option_dict['-csv'][:-4]+"_query_seq.fasta"
        out_file=open(option_dict['-out'],'w')
        for name in seq_name_list:
            if name in fasta_dic:
                out_file.write(">"+name+"\n")
                out_file.write(fasta_dic[name]+"\n")
            else:
                print ("No sequence identified with name: ", name)
        print ("Completed!")
    elif option_dict['-seq_type']=="1":
        if option_dict['-out']=="":
            if option_dict['-target']=="2":
                option_dict['-out']=option_dict['-csv'][:-4]+"_hit_seq.fastq"
            elif option_dict['-target']=="1":
                option_dict['-out']=option_dict['-csv'][:-4]+"_query_seq.fastq"
        fastq_parsing(option_dict['-seq'],seq_name_list, option_dict['-out'])

if __name__=="__main__":
    main()