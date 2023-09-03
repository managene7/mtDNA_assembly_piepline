#________________ option parse _______________________________
import sys 

args = sys.argv[1:]

option_dict={'-out':""}
help_description="""
________________________________________________________________________________


Usage;

If default value exists, corresponding value can be omitted.

-total_size     (Required)  total extraction size of the sequence 
-fasta_name     (Required)  fasta file name
-out            (option)    output file name
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

def convert_fasta(fasta, total_size, out):
    infile=open(fasta,'r')
    init=True
    out_list=[]
    while 1:
        line=infile.readline().strip()
        if line=="":
            seq="".join(seqs)
            seq_len=len(seq)
            out_list.append([seq_len, name, seq])
            break
        elif init:
            if ">" in line:
                name=line
                seqs=[]
                init=False
            else:
                print ("The first line is not a fasta name format.")
                quit()
        elif line[0]!=">":
            seqs.append(line)
        else:
            seq="".join(seqs)
            seq_len=len(seq)
            out_list.append([seq_len, name, seq])
            name=line
            seqs=[]
    out_list.sort()
    out_list.reverse()
    
    out_file=open(out, 'w')
    tot_size=0
    for cont in out_list:
        tot_size+=cont[0]
        if tot_size>int(total_size):
            out_file.write(cont[1]+"\n")
            out_file.write(cont[2]+"\n")
            return "%s sequence has parsed.\n" % str(total_size)
        else:
            out_file.write(cont[1]+"\n")
            out_file.write(cont[2]+"\n")
    return "-Alert!- Sequence is not enough: "+str(tot_size)+"/%s has parsed.\n" % str(total_size)
def main():
    if option_dict['-out']=="":
        option_dict['-out']=option_dict['-fasta_name']+"_top"+option_dict['-total_size']+"bp.fasta"
    result=convert_fasta(option_dict['-fasta_name'], option_dict['-total_size'], option_dict['-out'])
    print (result)

if __name__=="__main__":
    main()

        