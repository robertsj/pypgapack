"""
Crude testing using the examples and with reference output.
"""

import os, re
from math import fabs

def soft_equiv(x_ref, x) :
    if (fabs(x_ref - x) <= fabs(x_ref) * 1e-8) :
        return True
    elif (x_ref == 0.0) :
        return (fabs(x) <= 1e-8)
    else :
        return False

def test(index) :
    """
    Run example exampleINDEX.py
    """
    inp = 'example' + index + '.py'
    out = 'output/example' + index
    ref = 'output/example' + index + '_ref'
    # run the example
    print " testing ", inp
    os.system('rm ' + out) # remove any old run
    os.system('python ' + inp + ' >> ' + out)
    # read the output and reference and build arrays of lines
    f_out = open(out, 'r')
    f_ref = open(ref, 'r')
    lines_out = f_out.readlines()
    lines_ref = f_ref.readlines()
    if len(lines_ref) != len(lines_out) :
        print "Error: " + out + " line count mismatch."
        return False
    n = 0
    while n < len(lines_ref) :
        words_ref = lines_ref[n].split()
        words_out = lines_out[n].split()
        if len(words_ref) == 3 and words_ref[1] == 'Best' and words_ref[0] != 'The' :
            if not soft_equiv(float(words_ref[2]), float(words_out[2])) :
                print "Error: " + out + " in Iter #" + words_ref[0]
                return False
        elif len(words_ref) == 4 and words_ref[2] == 'Evaluation:' :
            if not soft_equiv(float(words_ref[3].strip('.')), float(words_out[3].strip('.'))) :
                print "Error: " + out + " in Best Evaluation" 
                return False
        elif len(words_ref) == 3 and words_ref[2] == 'String:' :
            n = n + 1
            while n < len(lines_ref) :
                s_r = lines_ref[n]
                s_o = lines_out[n]
                if len(s_r) == 0 :
                    if len(s_o) > 0 :
                        print "Error: " + out + " in line #" + str(n)
                        return False
                    return True
                elif s_r == '***Destroying PGA context***\n' :
                    if s_o != s_r :
                        print "Error: " + out + " in line #" + str(n)
                        return False
                    return True 
                s_r = re.sub(r"([\[\],#])", " ", s_r )
                s_r = re.sub(r"[\d][\:]", " ", s_r )
                s_r = s_r.split()
                s_o = re.sub(r"([\[\],#])", " ", s_o )
                s_o = re.sub(r"[\d][\:]", " ", s_o )    
                s_o = s_o.split()              
                for i in range (0, len(s_r)) :
                    if not soft_equiv(float(s_r[i]), float(s_o[i])) :
                        print "Error: " + out + " in Iter #" + words_ref[0]
                        return False   
                n += 1
        n += 1

n = 3 # number of examples so far

for i in range(1, n + 1) :
    if i < 10 :
        index = "0"+str(i)
    else :
        index = str(i)
    if (test(index)) :
        print "Test " + index + " passed."
    else :
        print "Test " + index + " failed."
    

