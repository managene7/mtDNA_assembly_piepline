#________________ option parse _______________________________
import sys 

args = sys.argv[1:]

option_dict={'-out':"", '-chop_size':"1000"}
help_description="""
________________________________________________________________________________


Usage;

If default value exists, corresponding value can be omitted.

-fasta      (required)  fasta file name
-chop_size  (option)    size for chopping (default: 1000 bp)
-out        (option)    output file name

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

def fasta_to_dic(fasta):
    infile=open(fasta,'r')
    out_dic={}
    init=True
    while 1:
        line=infile.readline().strip()
        if line=="":
            seq="".join(seqs)
            out_dic[name]=seq
            break
        elif init:
            if line.startswith(">"):
                name=line
                seqs=[]
                init=False
            else:
                print ("-Error!- The file is not a fasta format.")
                quit()
        else:
            if ">" not in line:
                seqs.append(line)
            else:
                seq="".join(seqs)
                out_dic[name]=seq
                name=line
                seqs=[]
    return out_dic

def chop_dic_fasta(dic_fasta, size, out_file):
    out=open(out_file,'w')
    
    for name,seq in dic_fasta.items():
        left=0
        serial=0
        iter=int(len(seq)/int(size))
        for i in range(iter):
            serial+=1
            right=left+int(size)
            out.write(name+"-"+str(serial)+"\n")
            if i==iter-1:
                out.write(seq[left:]+"\n")
            else:
                out.write(seq[left:right]+"\n")
            left=right

def main():
    if option_dict['-out']=="":
        option_dict['-out']=option_dict['-fasta']+"_%s.fasta" % option_dict['-chop_size']
    fasta_dic=fasta_to_dic(option_dict['-fasta'])
    chop_dic_fasta(fasta_dic, option_dict['-chop_size'], option_dict['-out'])
    print ("\nChopping sequences in %s file into %s bp has been completed!\n" % (option_dict['-fasta'], option_dict['-chop_size']))

if __name__=="__main__":
    main()
