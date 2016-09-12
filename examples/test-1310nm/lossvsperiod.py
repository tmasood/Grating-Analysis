#!/usr/bin/env python3

import re
import os

def main():
    # Calculate Loss .vs. Grating Period for the structure 
    print('Calculate Loss .vs. Period')
    gdatfile = open("lossvsprd.dat", "w")

    periodAint = 1800
    while (periodAint < 4200):
        cmd = './../../gratanal test-1310nm-ex1 2>&1 | tee tmp-file'
        os.system(cmd)

        # Open file
        infile = open("test-1310nm-ex1", "r")
        outfile = open("test-1310nm-ex1-tmp", "w")
        tmpfile = open("tmp-file", "r")

        lines1 = tmpfile.readlines()
        tmpfile.close()

        # look for loss
        for line in lines1:
            stregexp = re.compile(r'.+loss\s=\s(-?\d+\.\d+)')
            match = stregexp.match(line)
            if match:
                loss = match.group(1)
                print loss
    
        lines = infile.readlines()
        infile.close()

        # look for PRD
        for line in lines:
            stregexp = re.compile(r'PRD\s=\s(\d+\.\d+)')
            match = stregexp.match(line)
            if match:
                period = match.group(1)
                period = float(period)+0.0001
                periodA = period * 10000
                gpdata = str(periodAint) + "    " + str(loss) + "\n"
                gdatfile.write(gpdata)
                periodAint = periodA
                print periodAint                
                line = "PRD = " + str(period) + "\n"

            outfile.write(line)
            
        outfile.close()
        os.rename("test-1310nm-ex1-tmp","test-1310nm-ex1")
    
    gdatfile.close()
    
main()
