import sys,os
import numpy as np
import build_aform as create

# read sequence and SS
fh = open(sys.argv[1])
outname=sys.argv[2]
info = fh.readlines()
fh.close()
seq = info[0][:-1]

# create single stranded extended structure
ss = create.create_ARNA(seq,"")
fhw = open("%s_extended.pdb" % (outname),"w")
fhw.write(ss)
fhw.close()
print("# A-form single stranded structure is %s_extended.pdb" % (outname))

# create double helical templates
for j in range(1,len(info)):
    #print(info[j])
    sec_struc = info[j][:-1]
    # read ss and perform checks here and there
    assert(len(sec_struc)==len(seq))
    myseq1 = ""
    myseq2 = ""
    is_started1=False
    is_started2=False
    for k,el in enumerate(sec_struc):
        if(el=="("):
            myseq1 += seq[k]
            if(is_started1==False):
                assert(len(sec_struc[k:].split(".("))==1)
                # find starting residue and atom
                os.system("awk '$5 == %d' %s_extended.pdb | head -1 | awk '{print $2}' > tmp1" % (k+1,outname))
                residue_shift1 = (k)
                atom_shift1 = np.loadtxt("tmp1",dtype=int)-1
                is_started1=True
                
        if(el==")"):
            myseq2 += seq[k] 
            if(is_started2==False):
                assert(len(sec_struc[k:].split(".)"))==1)
                os.system("awk '$5 == %d' %s_extended.pdb | head -1 | awk '{print $2}' > tmp2" % (k+1,outname))
                residue_shift2= (k)
                atom_shift2 = np.loadtxt("tmp2",dtype=int)-1

                is_started2=True
    print(atom_shift1,atom_shift2,residue_shift1,residue_shift2)
    #cmd="python create_arna4.py %s %s %d %d %d %d > %s_template_%d.pdb" % (myseq1,myseq2,residue_shift1, atom_shift1, outname,(j))
    assert(len(myseq1)==len(myseq2))
    print("# found helix of lenght %d (%s %s)" % (len(myseq1),myseq1,myseq2))
    term5a,term3a,term5b,term3b= False, False,False, False
    if(residue_shift1==0):
        term5a=True
    if(residue_shift2+len(myseq2)==len(seq)):
        term3b = True
        
    print(term3b,residue_shift2,len(myseq2),len(seq))
    ss = create.create_ARNA(myseq1,myseq2,\
                       atom_shift1=atom_shift1,residue_shift1=residue_shift1,\
                            atom_shift2=atom_shift2,residue_shift2=residue_shift2,\
                            term5a=term5a,term3a=term3a,term5b=term5b,term3b=term3b)
    fhw = open("%s_template_%d.pdb" % (outname,j),"w")
    fhw.write(ss)
    fhw.close()
    print("# A-form double stranded template %s_template_%d.pdb" % (outname,j))
