#________________ option parse _______________________________
import sys 
args = sys.argv[1:]
option_dict={'-out':"", '-min_len':'20000', '-out_type':"1"}
for i in range(len(args)):
    if args[i].startswith("-"):
        try:
            option_dict[args[i]]=args[i+1]
        except:
            if args[0]=="-help" or args[0]=="":
                print ("""
________________________________________________________________________________
Usage;
If default value exists, corresponding value can be omitted.
-fastq      fastq sequence file name
-out_type   1: fastq (default), 2: fasta, 3: both
-min_len    20000 (bp) (default)
-out        Enter output name (option)
________________________________________________________________________________
""")
                quit()

def fastq_parsing(fastq, min_len, out_type, out):
    print ("Size selection of fastq reads over %s bp..\n" % option_dict['-min_len'])
    infile=open(fastq,'r')
    total_length=0
    if out_type=="1" or out_type=="2":
        out_file=open(out,'w')
    if out_type=="3":
        out_file1=open(out+".fastq",'w')
        out_file2=open(out+".fasta",'w')
    while 1:
        line1=infile.readline().strip()
        line2=infile.readline().strip()
        line3=infile.readline().strip()
        line4=infile.readline().strip()
        
        if line1=="":
            break
        elif len(line2)>=int(min_len):
            total_length += len(line2)
            if out_type=="1":
                out_file.write(line1.split()[0]+"\n")
                out_file.write(line2+"\n")
                out_file.write(line3+"\n")
                out_file.write(line4+"\n")
            elif out_type=="2":
                out_file.write(">"+line1.split()[0][1:]+"\n")
                out_file.write(line2+"\n")
            elif out_type=="3":
                out_file1.write(line1.split()[0]+"\n")
                out_file1.write(line2+"\n")
                out_file1.write(line3+"\n")
                out_file1.write(line4+"\n")

                out_file2.write(">"+line1.split()[0][1:]+"\n")
                out_file2.write(line2+"\n")
    return total_length

def main():
    if option_dict['-out']=="":
        if option_dict['-out_type']=="1":
            option_dict['-out']=option_dict['-fastq']+"_min"+option_dict['-min_len']+".fastq"
        elif option_dict['-out_type']=="2":
            option_dict['-out']=option_dict['-fastq']+"_min"+option_dict['-min_len']+".fasta"
        elif option_dict['-out_type']=="3":
            option_dict['-out']=option_dict['-fastq']+"_min"+option_dict['-min_len']
        else:
            print ("Please choose '-out_type' option. ")
            quit()
    total_len=fastq_parsing(option_dict['-fastq'], option_dict['-min_len'], option_dict['-out_type'], option_dict['-out'])
    total_len=round(total_len/1000000,1)
    print ("Total %s Mbp reads are isolated.\n" % str(total_len))

if __name__=="__main__":
    main()
