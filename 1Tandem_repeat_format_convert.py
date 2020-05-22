#usage: python xx.py input_name out_name
import sys
import re
pattern = re.compile("[0-9]{1,}")

outfile = open(sys.argv[2],'w')
with open(sys.argv[1]) as fh:
    block = {}
    for row in fh:
        if row.strip().startswith("Sequence:"):
            chr_num = row.strip().split(" ")[1]
            block[chr_num]=[]
        else:
            repeat_start = row.strip().split(" ")[0]
            if re.match(pattern, repeat_start):
                block[chr_num].append(row.strip().split(" "))
                outfile.write(chr_num+"\t"+"\t".join(row.strip().split(" "))+"\n")
print 'ok'
            
