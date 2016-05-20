import sys
import gzip

def index_headers(bank_file_name):
    headers= []
    headers.append("NULL") # reads are 1-indexed we add this dummy value
    if "gz" in bank_file_name:
        sequencefile=gzip.open(bank_file_name,"r")
    else: 
        sequencefile=open(bank_file_name,"r")
        
    # i=1
    for line in sequencefile.readlines():
        if line[0]=='>':
            # print line[1:],
            headers.append(line.rstrip()[1:])            #
            # print i
            # print line[1:],
            # print headers[i]
            # i+=1
    return headers
    


def sequence_sizes(bank_file_name):
    sizes= []
    sizes.append(0) # reads are 1-indexed we add this dummy value
    
    if "gz" in bank_file_name:
        sequencefile=gzip.open(bank_file_name,"r")
    else: 
        sequencefile=open(bank_file_name,"r")
        
    for line in sequencefile.readlines():
        if line[0]!='>':
            # print line[1:],
            sizes.append(len(line))
    return sizes
            

def convert_SRC_linker_output(headers, sizes, k, SRC_linker_output_file_name):
    print "qseqid\tsseqid\tevalue\tpident"
    if "gz" in SRC_linker_output_file_name:
        sequencefile=gzip.open(SRC_linker_output_file_name,"r")
    else: 
        sequencefile=open(SRC_linker_output_file_name,"r")
        
    #31:1246-3 479-3 1043-3 820-3 
    
    for line in sequencefile.readlines():
        if line[0]=='#': #header
            continue
        line=line.rstrip()
        query_read_id=int(line.split(':')[0])
        targets=line.split(':')[1].split(' ')
        for target in targets:
            target_read_id=int(target.split('-')[0])
            if target_read_id == query_read_id: continue
            coverage = 200*k*int(target.split('-')[1])/float(sizes[query_read_id]+sizes[target_read_id])
            print headers[query_read_id]+"\t"+headers[target_read_id]+"\t0.0\t%.4g"%(coverage)
        
        
headers=index_headers(sys.argv[1])
sizes=sequence_sizes(sys.argv[1])
k=int(sys.argv[3])
convert_SRC_linker_output(headers,sizes,k, sys.argv[2])
