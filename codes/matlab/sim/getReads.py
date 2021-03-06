#!/usr/bin/python

import sys
from numpy import *
import random
import numpy.random as nrd

from optparse import OptionParser
parser = OptionParser(usage="-r REF expressionFile1.exp [expressionFiles2.exp]\n\n        Program generates reads from fasta file based on read counts provided in the expression files files (generated by setExpression2.py). The reads are saved into expressionFile1.fastq")
#parser.add_option("-f", "--fold", dest="fold", default=2.0, type="float")
#parser.add_option("-N", "--molN", dest="N", default=1000000, type="int")
#parser.add_option("-s", "--skip", dest="skip", default=1, type="int")
#parser.add_option("-v", "--verbose", default=False, dest="verbose", action="store_true", help="Print out separate histograms")
#parser.add_option("-e", "--empty", default=False, dest="empty", action="store_true", help="")
#parser.add_option("-l", "--logged", default=True, dest="logged", action="store_false")
parser.add_option("-r","--ref", dest="ref", type="string", help="reference fasta file with transcript names matching the names in .exp files")

(options, args) = parser.parse_args()

# quality probs {{{
# starting at ascii 30 data from SRR039631.fastq
qualProb = [0.0, 0.0, 0.0, 991.0, 66931.0, 0.0, 16354.0, 41889.0, 0.0, 78957.0, 32727.0, 22735.0, 74530.0, 82952.0, 42842.0, 5700.0, 58968.0, 139482.0, 50500.0, 76961.0, 31504.0, 171013.0, 223970.0, 16481.0, 199051.0, 257191.0, 509415.0, 77964.0, 508007.0, 1127517.0, 5918205.0, 166107.0, 829.0, 227.0, 0.0];
qPSum=sum(qualProb);
qP = zeros(len(qualProb));
for i in xrange(1,len(qualProb)):
   qP[i] = qP[i-1] + qualProb[i]/qPSum;
#}}}
def getRead(rd,fl): #{{{
   if fl:
      rd = rd[::-1]; #  REALLYY?
      rd2=""
      for i in xrange(len(rd)):
         if rd[i] == 'A' or rd[i] == 'a': rd2+='T';
         elif rd[i] == 'T' or rd[i] == 't': rd2+='A';
         elif rd[i] == 'G' or rd[i] == 'g': rd2+='C';
         elif rd[i] == 'C' or rd[i] == 'c': rd2+='G';
         else: rd2+='N'
      rd = rd2;
   qual="";
   read=""
   for i in xrange(len(rd)):
      x = nrd.random();
      q = 30;
      while qP[q-1]>x:q-=1;
      while qP[q]<x:q+=1;
      qch = q+30;
      qual += chr(qch);
      q = qch - 33;
      pm = 10**(q/-10)
      x = nrd.random();
      if x<pm:
         errCh = random.sample(['A','T','G','C','N'],1)[0]
         while errCh == rd[i]:
            errCh = random.sample(['A','T','G','C','N'],1)[0]
         read += errCh;
      else:
         read += rd[i];

   return (read,qual);
#}}}

#fastaN = "/localhome/work/refHG/ensembl/ensemblGenes.fasta"
if options.ref:
   fastaN=options.ref;
else:
   sys.exit("Please provide reference file");
#fastaN = "ensemblGenes.fasta"
rLen = 50;

for exp in args:
   print "Processing file: ",exp;
   if exp[-5:] == ".exps": exp=exp[:-5];
   
   exps = []
   mapp = {};
   inF=open(exp+".exps","r");
   i = 0 ;
   norm = 0.;
   for line in inF:
      if line[0] == '#':continue;
      lA=line.split();
      ct = int(lA[0]);
      #d1 = float(lA[1]);
      #d2 = float(lA[2]);
      #l = float(lA[4]);
      name = lA[3];
      exps.append([ct,name]);
      mapp[name]=i;
      i+=1;
   inF.close();
   M = i;

   chrF = open(fastaN,"r");
   outF = open(exp+".fastq","w");

   reads = [];
   trN = ""
   trS = ""
   strand = 0;
   readT = 0
   counter=729;
   for line in chrF:
      if readT>counter:
         print "  ",readT;
         counter*=3;
      if len(reads)>1000:
         samp = len(reads)/2;
         for j in xrange(samp):
            i = random.randint(0,len(reads)-1);
            fl = int(random.getrandbits(1));
            outF.write("@"+reads[i][1]+"-"+str(fl)+"-"+str(readT)+"\n");
            readT+=1;
            (read,qual) = getRead(reads[i][0],fl);
            outF.write( read +"\n+\n"+qual+"\n");
            reads.pop(i);

      if line[0] == '>':
         if trN!="":
            trL = len(trS);
            if strand: trS += "A"*rLen;
            else: trS ="A"*rLen + trS;
            ct = exps[mapp[trN]][0];
            for i in xrange(ct):
               st = random.randint(0,trL-1);
               reads.append([trS[st:st+rLen],trN+"-"+str(st)]);

         trN = line.split()[0][1:];
         if line.split()[2][-2] == '-':strand = 0;
         else: strand = 1;
         trS = "";
      else: trS+=line.rstrip();
   
   for i in xrange(len(reads)):
      fl = int(random.getrandbits(1));
      outF.write("@"+reads[i][1]+"-"+str(fl)+"-"+str(readT)+"\n");
      readT+=1;
      (read,qual) = getRead(reads[i][0],fl);
      outF.write( read +"\n+\n"+qual+"\n");
   
   outF.close();
   chrF.close();
   print "Reads written: ",readT;

