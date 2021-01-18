import sys
import numpy as np

# z increment and rotation for A-RNA
z_increment = 2.81
rot = 0.5707226654

# read database
def readline(line):

    res_type = line[17:20].strip()

    res_num = 1
    atom_type = line[12:16].strip()
    atom_num = int(line[6:11].strip())

    X = float(line[30:38])
    Y = float(line[38:46])
    Z = float(line[46:54])
    vv = [atom_num,atom_type,res_type,res_num,X,Y,Z]
    return vv


# write pdb to stirng
def pdb_string(dline,atom_shift,residue_shift):
    
    insert = ""

    anum = dline[0] + atom_shift
    atype = dline[1]
    if("H" in atype):
        bfac,occ  = 0,0
    else:
        bfac,occ = 1,1
        
    res_type = dline[2]
    res_num = dline[3] + residue_shift
    x = dline[4]
    y = dline[5]
    z = dline[6] 
    chain = " "
    
    if(len(atype)==4):
        string = "%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.0f%6.0f\n" \
            % ("ATOM",anum,atype," ",res_type,chain,res_num,insert,x,y,z,occ,bfac)
    else:
        string = "%-6s%5d  %-3s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.0f%6.0f\n" \
            % ("ATOM",anum,atype," ",res_type,chain,res_num,insert,x,y,z,occ,bfac)
    return string

# translate and rotate template 
def tr_rot(rr,idx,sign):
    for k in range(len(rr)):
        rr[k][6] *= sign
        rr[k][6] += idx*z_increment
        x = rr[k][4]
        y = rr[k][5]
        rr[k][4] = x*np.cos(sign*idx*rot) - y*np.sin(sign*idx*rot)
        rr[k][5] = sign*(x*np.sin(sign*idx*rot) + y*np.cos(sign*idx*rot))
        
    

def create_ARNA(s1,s2,atom_shift1=0,residue_shift1=0,atom_shift2=0,residue_shift2=0,term5a=True,term3a=True,term5b=True,term3b=True):

    # check that sequences are the same
    #assert(len(s1)==len(s2))
    print("# lenght of strand 1:%d" % len(s1))
    print("# lenght of strand 2:%d" % len(s2))
    ll1 = len(s1)
    string = ""
    # first strand
    for i in range(ll1):
        el1 = s1[i]
        # read 3' or 5' end
        if(term3a==True):
            if(i==ll1-1):el1 += "3"
        if(term5a==True):
            if(i==0):el1 += "5"
        # read residue in database
        rr1 = [readline(line) for line in open("pdbs/%s.pdb" % el1)]
        # translate and rotate 
        tr_rot(rr1,i,+1)
        # add to string 
        for el in rr1:
            string += pdb_string(el,atom_shift1,residue_shift1)
            
        atom_shift1 += len(rr1)
        residue_shift1 += 1

    #string += "TER\n"
    if(atom_shift2==0):
        atom_shift2 = atom_shift1
    if(residue_shift2==0):
        residue_shift2 = residue_shift1

    ll2 = len(s2)
    
    # second strand
    for i in range(ll2):
        el1 = s2[i]
        if(term3b==True):
            if(i==ll2-1): el1 += "3"
        if(term5b==True):
            if(i==0):el1 += "5"
        rr1 = [readline(line) for line in open("pdbs/%s.pdb" % el1)]
        tr_rot(rr1,ll1-1-i,-1.)
        for el in rr1:
            string += pdb_string(el,atom_shift2,residue_shift2)
        atom_shift2 += len(rr1)
        residue_shift2 += 1
    string += "TER\n"

    return string

def main():
    
    s1 = sys.argv[1]
    if(len(sys.argv)==2):
        print(create_ARNA(s1,""))
    else:
        s2 = sys.argv[2]
        print(create_ARNA(s1,s2))


    
if __name__ == '__main__':
   main()
