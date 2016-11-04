#!/bin/python
from __future__ import print_function
from math import *
import sys, re

GAP='-'

# Assumes only two fasta entries
def read_fasta(fasta):
    a=''
    a_name=''
    b=''
    b_name=''
    state=0
    for line in open(fasta, 'r'):
        line = line.rstrip('\n')
        line = line.rstrip('\r')
        if 0 == state:
            if None == re.search('^>', line):
                print('error: not valid fasta format')
            a_name = line[1:]
            state = state + 1
        elif 1 == state:
            if None == re.search('^>', line):
                a = a + line
            else:
                state = state + 1
                b_name = line[1:]
        else:
            b = b + line
    return (a,a_name,b,b_name)

def read_scoring(scoring):
    state = 0
    header = []
    sm = {}
    for line in open(scoring, 'r'):
        line = line.rstrip('\n')
        line = line.rstrip('\r')
        if None == re.search('^#', line):
            a = line.split(' ')
            if 0 == state:
                header = a
                state = state + 1
            else:
                if len(a) != 1 + len(header):
                    print('error: the number of entries did not match the header')
                    sys.exit(2)
                for i in range(len(header)):
                    sm[a[0] + header[i]] = int(a[i+1])
   
    return sm

def globalAlign(a,b,sm):
    n=len(a) # rows
    m=len(b) # cols
    d=-5 # gap penalty

    # Initialize grid of size len(a)+1 x len(b)+1 with 0
    grid = [[0 for j in range(m+1)] for i in range(n+1)]

    # Initialize first column and row of grid by multiplying distance by gap penalities
    for i in range(n+1):
        grid[i][0] = i * d
    for j in range(m+1):
        grid[0][j] = j * d
        
    # Create scoring matrix
    for i in range(1,n+1):
        for j in range(1,m+1):
            amino_acid = a[i-1] + b[j-1]
            match = grid[i-1][j-1] + sm[amino_acid]
            delete = grid[i-1][j] + d
            insert = grid[i][j-1] + d
            grid[i][j] = max(match, insert, delete)

    x = ""  # align sequence a
    y = ""  # align sequence b

    # Start with the last amino acid in the sequence
    i = n 
    j = m 

    # create the alignment seqeunces by choosing optimal path with direction being diagonal, up, or left
    while (i > 0 and j > 0):
        amino_acid = a[i-1] + b[j-1]
        path = grid[i][j]
        pathDiag = grid[i - 1][j - 1]
        pathUp = grid[i][j - 1]
        pathLeft = grid[i - 1][j]

        # choose optimal direction based on highest score
        if (path == pathDiag + sm[amino_acid]):
            # choose diagonal path
            x = a[i-1] + x
            y = b[j-1] + y
            # update path
            i -= 1
            j -= 1
        elif (path == pathLeft + d):
            # choose left path
            x = a[i-1] + x
            y = GAP + y
            # update path
            i -= 1
        elif (path == pathUp + d):
            # choose up path
            x = GAP + x
            y = b[j-1] + y
            # update path
            j -= 1

    # fill in beginning of sequence a and b with gaps if needed
    while (i > 0):
        x = a[i-1] + x
        y = GAP + y
        i -= 1

    while (j > 0):
        x = GAP + x
        y = b[j-1] + y
        j -= 1

    # calculate global score
    score = 0
    for i in range(len(x)):
        amino_acid = x[i] + y[i]

        # add score when amino acids match 
        if x[i] == y[i]:
            score += sm[amino_acid]

        # apply gap penalty if there is a gap
        elif x[i] == GAP or y[i] == GAP:          
            score += d
    
        # get score if amino acids do not match and do not contain gaps
        elif x[i] != GAP and y[i] != GAP and x[i] != y[i]: 
            score += sm[amino_acid]

    return (x,y,score)

def localAlign(a,b,sm):
    n=len(a) # rows
    m=len(b) # cols
    d=-5 # gap penalty

    # Initialize grid of size len(a)+1 x len(b)+1 with 0
    grid = [[0 for j in range(m+1)] for i in range(n+1)]

    # Need tracer matrix in order to store the optimal path
    tracer = [[0 for j in range(m+1)] for i in range(n+1)]

    score = 0  
    # Create scoring matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            amino_acid = a[i-1] + b[j-1]
            match = grid[i-1][j-1] + sm[amino_acid]
            delete = grid[i][j-1] + d   
            insert = grid[i-1][j] + d     
            grid[i][j] = max(0, insert, delete, match)

            # reached end of path
            if grid[i][j] == 0:
                tracer[i][j] = 0

            if grid[i][j] == insert:
                tracer[i][j] = 1

            if grid[i][j] == delete:
                tracer[i][j] = 2

            if grid[i][j] == match:
                tracer[i][j] = 3

            if grid[i][j] >= score:
                score = grid[i][j];
                iNew = i
                jNew = j
                
    
    x = ""  # align sequence a
    y = ""  # align sequence b

    # Start with the highest score in the grid
    i = iNew 
    j = jNew 
    
    # create the alignment seqeunces by choosing optimal path with direction being diagonal, up, or left
    while tracer[i][j] != 0:
        if tracer[i][j] == 3:
            x = a[i-1] + x
            y = b[j-1] + y
            # update path
            i -= 1
            j -= 1

        elif tracer[i][j] == 2:
            x = GAP + x
            y = b[j-1] + y
            # update path
            j -= 1

        elif tracer[i][j] == 1:
            x = a[i-1] + x
            y = GAP + y
            # update path
            i -= 1

    # calculate local score
    score = 0
    for i in range(len(x)):
        amino_acid = x[i] + y[i]

        # get score when amino acids match 
        if x[i] == y[i]:
            score += sm[amino_acid]

        # apply gap penalty if there is a gap
        elif x[i] == GAP or y[i] == GAP:          
            score += d
    
        # get score if amino acids do not match and do not contain gaps
        elif x[i] != GAP and y[i] != GAP and x[i] != y[i]: 
            score += sm[amino_acid]

    return (x,y,score)

def main():
    # get arguments
    if 3 != len(sys.argv):
        print('error: two input files required')
        sys.exit(2)
    fasta = sys.argv[1]
    scoring = sys.argv[2]

    # read in input files
    (a,a_name,b,b_name)=read_fasta(fasta)
    (sm)=read_scoring(scoring)

    (ag,bg,sg)=globalAlign(a,b,sm)
    (al,bl,sl)=localAlign(a,b,sm)

    # print output
    print("#GLOBAL SCORE=", sg, sep=' ')
    print(a_name,ag)
    print(b_name,bg)
    print("#LOCAL SCORE=", sl, sep=' ')
    print(a_name,al)
    print(b_name,bl)

if __name__ == '__main__':
    main()
