import os,sys
newpath='/Users/uqmkel13/Documents/workspace/MainRepository'
sys.path.append(newpath)

from DataAnalysisTools import genotypeFileTools
import sqlite3


from DataAnalysisTools import datatools
from DataAnalysisTools import genotypeFileTools

conn = sqlite3.connect('/Users/uqmkel13/Documents/Projects/ReproCSIRO/CRCReproData/CRCGenotypes.db')
cur=conn.cursor()
genotypes=genotypeFileTools.genomeFileTools(conn)
NumberAnimals=10
NumberSNPs=10

#select SNPs
#select animals
snpNames="select SNPName from SNPLocations limit %s"%NumberSNPs
cur.execute(snpNames)
x=cur.fetchall()
snpNames=[i[0] for i in x ]


genotypes=genotypeFileTools.genomeFileTools(conn)
genotypes.createSelectedMarkerList(["%s"%i for i in snpNames])

Ans="select barcode from genotypes limit %s "%NumberAnimals

cur.execute(Ans)
x=cur.fetchall()
Animals=[i[0] for i in x]
print Animals
for i in x:
    x=genotypes.getGenotype(i)
    print x 
