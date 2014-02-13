import os,csv,sys
from random import choice
from string import split
import fileinput
from datatools import makeTable
from zlib import compress,decompress


#from zlib import compressobj
import zlib
#genotype file tools
#
# Create Genable file
# 
# Recode file
# select Rows 
# Select columns
#
def recodeWombat(t):
#    recodeDict={"AA":'0',"AB":'1',"BB":'2'}
#    tmp=t.split(",")
#    tmp=[recodeDict[i] for i in tmp]
#    tmp="".join(tmp)
#    return "".join(tmp)
    print t 
    return t



#set up a table containing the positions of the markers to select
def selectMarkers():
    
    return

class genomeFileTools:
    def __init__(self,connection=None):
        self.con=connection
        if not connection==None:
            self.cur=self.con.cursor()
            self.con.create_function("recodeWombat", 1, recodeWombat)
      #  self.con.text_factory = str
            self.con.create_function("pullSelectedMarkers", 1, self.pullSelectedMarkers)
            self.con.create_function("countHetSelectedMarkers", 1, self.countHetSelectedMarkers)

        pass
    def DBtoCRCformat(self,outfile,snpmap):
#        1 read SNPMAP
#        2 create markerList in order
#        3 push to 
        
        return
    
    def calcAF(self):
        """Calculate Allele frequency for a data set """
        pass
    def pullSelectedMarkers(self,t):
        return "".join([t[i] for i in self.selectedMarkers])

    def countHetSelectedMarkers(self,t):  
        return sum([1 for i in self.selectedMarkers if t[i]=='1'])


    def runQuery(self,query):
        self.cur.execute(query)
        return self.cur.fetchone()

    def createSelectedMarkerList(self,markerList):
        sql='select Position,SNPName from GenotypesFileMap where SNPName in (%s) and GTFileNumber=1'
        SNPdict={}
        for i in markerList:
            if not SNPdict.has_key(i):
                SNPdict[i]=0
            SNPdict[i]+=1
            
        
        markerListT=["'%s'"%i for i in markerList]
       
        self.cur.execute(sql%",".join(markerListT))
        x=self.cur.fetchall()
        SNPdict={}
        for i in x:
            SNPdict[i[1]]=i[0]
        
        self.selectedMarkers=[SNPdict[i] for i in markerList]
        print len(x),len(markerList),'If these match all is good'
        if not len(x)==len(markerList):
            print 'marker list size wrong either duplicate markers or markers not in database'
            sys.exit()
        
    
    def selectEvenlySpacedMarkers(self,NMarkersToCheck):
        allmrks="""SELECT Position from GenotypesFileMap """    
        self.cur.execute(allmrks)
        allMarkers=[i[0] for i in self.cur.fetchall()]
        stepSize=(len(allMarkers)/NMarkersToCheck)
        selmarkers= allMarkers[::stepSize]
        self.selectedMarkers=selmarkers[0:NMarkersToCheck]    
    #The join in TMP may sometimes need some string manipulation not sure why it currently does not
    def getGenotype(self,barcodes):
        tmp=",".join(["'%s'"%i for i in barcodes])
        phensql="""SELECT barcode,pullSelectedMarkers(genotype) from Genotypes where barcode in (%s) """%tmp
      #  print phensql
       
        self.cur.execute(phensql) 
        x=self.cur.fetchall() 
        return x
    def loadGenotypes(self,genotypesFile,cur):
        self.recodeDict={"AA":'0',"AB":'1',"BB":'2','--':'9'}
        
        recodeDict=self.recodeDict
        print "Loading map"
        f=open(genotypesFile,'rb')
        genotypesData=csv.reader(f)
        #add file to List of genotypeFiles
        cur.execute("INSERT INTO GenotypesFiles (GenotypeFileName) VALUES (?)",(genotypesFile,))
        
        #Get RowNumber for this file
        cur.execute("select ROWID from GenotypesFiles where GenotypeFileName='%s'"%genotypesFile)
        x=cur.fetchone()
        GTfileNumber= x[0] 

        #read and load snps into GenotypesFileMap
        header=genotypesData.next()
        #need to remove the ID column for mapping
        del header[0]
        SNPtupForInsert=[(GTfileNumber,header[i],i) for i in range(len(header))]
        #Insert  into map
        SQL="insert into  GenotypesFileMap (GTfileNumber,SNPName,Position) VALUES (?,?,?)"
        cur.executemany(SQL,SNPtupForInsert)
        print "Loading data"
        counter=0
        #load genotypes
        for i in genotypesData:
            anID=i[0]
            del i[0]
            SQL="insert into  Genotypes (barcode,GTfileNumber,genotype) VALUES (?,?,?)"
            tmp="".join([recodeDict[j] for j in i])
            cur.execute(SQL,(anID,GTfileNumber,tmp))
            counter+=1
            if counter%100==0:
                print 'up to row',counter
                
        f.close()
        self.con.commit()
        pass
    
   
    def calcGEBV(self,a,b):
         return sum([a[i]*float(b[i]) for i in range(len(a))])
     #takes a genotype file csv delimited makes tables and stores in sqllite
    def createDataStructureAndLoad(self,genotypesFile,cur):
        #GenotypesFiles (FileName)
        #GenotypesFileMap (GenotypeFileRowID,SNPname,SNPpostion)
        #Genotypes (ID,GenotypeFileRowID,genotypes)
        print "Setting up tables"
        makeTable('GenotypesFiles', ["GenotypeFileName"],["TEXT"], cur)
        makeTable('GenotypesFileMap',['GTfileNumber','SNPName','Position'], ["INTEGER","TEXT","INTEGER"], cur)
        makeTable('Genotypes', ["barcode","GTfileNumber","genotype"],["TEXT","INTEGER","TEXT"], cur)
        self.loadGenotypes(genotypesFile,cur)
        pass

    def writePLINKFromDB(self,barcodes,markers,fileName):
        recodeDict={'0':"A A",'1':"A B",'2':'B B','9':'0 0'}
        mrks=self.createSelectedMarkerList(markers)
        
        fout=open('%s.ped'%fileName,'wb')
        fout2=open('Phen%s'%fileName,'wb')
       # fout3=open('gen%s.txt','wb')
        fout2.write("id\tsex\tt1\n")
        for an in barcodes:
            i=self.getGenotype([an])[0]
     #       print i,len(i[0])
            tmp=[]
            tmp.extend(['1',i[0],'0','0','1','%s'%choice([1,2])])
            tmp1=[recodeDict[j] for j in i[1]]
            tmp.extend(tmp1)
#MK LATE 16/10
#            for j in i[1]:
#                gty=recodeDict[j]
#                tmp.append(gty)

     #       print i[0],len(i[1])
            fout2.write("%s\t%s\t%s\n"%(tmp[1],tmp[4],tmp[5]))
            tmp=" ".join(tmp)
            fout.write("%s\n"%tmp)
            
        fout.close()    
        fout=open('%s.map'%fileName,'wb')
 #       fout.write("chrom name position \n")
        
        sqlMap="""Select chromosome, SNPName , '0' as CMposition,BPposition from SNPlocations where SNPName='%s'"""
        for i in markers:
            self.cur.execute(sqlMap%i)
            xx=self.cur.fetchone()
            print xx
            fout.write("%s %s %s %s \n"%(xx[0],xx[1],xx[2],xx[3]))


    def writeGenableFromDB(self,barcodes,markers,fileName):
        recodeDict={'0':"A/A",'1':"A/B",'2':'B/B','9':'0/0'}
        mrks=self.createSelectedMarkerList(markers)
        anns=self.getGenotype(barcodes)
        fout=open(fileName,'wb')
        fout2=open('Phen%s'%fileName,'wb')
        fout2.write("id\tsex\tt1\n")
        for i in anns:
            tmp=[]
            tmp.extend(['1',i[0],'0','0','1','%s'%choice([1,2])])
            for j in i[1]:
                gty=recodeDict[j]
                tmp.append(gty)
     #       print i[0],len(i[1])
            fout2.write("%s\t%s\t%s\n"%(tmp[1],tmp[4],tmp[5]))
            tmp=" ".join(tmp)
            fout.write("%s\n"%tmp)
            
        fout.close()    
        fout=open('SNPmap.txt','wb')
#        fout.write("chrom name position \n")
        
        sqlMap="""Select chromosome, SNPName , BPposition from SNPlocations where SNPName='%s'"""
        for i in markers:
            self.cur.execute(sqlMap%i)
            xx=self.cur.fetchone()
            print xx
            fout.write("%s %s %s %s \n"%(xx[0],xx[1],xx[2]))
#        for i in f:
#            tmp=i.split()     
#            fout.write("%s %s %s \n"%(tmp[2],tmp[0],tmp[3]))

        #NEED TO Print an accompanying map file for use in genable
        
#makeTable(self,tableName,fieldTypes,c)    

# PUSH to db
# create table for markers
# create index of locations of each marker
# create table for genotypes 
# Two columns id and genotype

#
    
    
def convertToGenable(fileNameIn,fileNameOut,recodeDict,missvalue):
    f=open(fileNameIn,'rU')
    fcsv=csv.reader(f)
    fout=open(fileNameOut,'wb')
    
    counter=0
    header=fcsv.next()
 #   foutcsv.writerow(header)
    idlist=[]
    print 'Number of fields in header',len(header)
    for i in fcsv:
        counter+=1
        tmp=[]
        tmp.extend(['1',i[0],'0','0','1','%s'%choice([1,2])])
        idlist.append(i[0])
        for j in i[1:]:
            #print j
            gty=missvalue
            if recodeDict.has_key(j):
                gty=recodeDict[j]
            tmp.append(gty)
        if counter==1: print 'Number of fields in line',len(tmp)
     #   sys.exit()
        tmp=" ".join(tmp)
        fout.write("%s\n"%tmp)
     #   if counter>10:
     #       sys.exit()
    
    f.close()
    fout.close()
    return idlist

def convertLongtoCSV(genotypeFile,mapfile,outputfile,recodeDict,postions=[2,3]):
    skip=10
    fout=open(outputfile,'wb')
    foutcsv=csv.writer(fout)   
    snpHeader=[]
    snpDict={}
    f=open(mapfile,'rb') 
    counter=0

    for line in f:
        counter+=1
        tmp=line.split()
        snpHeader.append(tmp[0])
        snpDict[tmp[0]]=0

    f.close()

    snpHeader.insert(0,'ID')

    print 'Markers in MAP',len(snpHeader)
    foutcsv.writerow(snpHeader)
    f=open(genotypeFile,'rb') 
    for i in range(skip):
        tmp=f.readline()
        print tmp[:-1]
    AnimalID=''
    counter=0
    from copy import deepcopy
    baseGenotypes={}
    for i in snpHeader:
        baseGenotypes[i]='--'
        
    genotypeDict=deepcopy(baseGenotypes)
    #genotype
    
    for line in f:
        tmp=line[:-1].split()
      #  print tmp
        if counter>0 and AnimalID!=tmp[1]:
            genotypeDict['ID']=AnimalID
#            print 'Write animal out'
            thisLine=[genotypeDict.has_key(i) for i in snpHeader]
 #           print 'MatchingReadyToGo',sum(thisLine)
  #          print 'Not Matching',sum([not genotypeDict.has_key(i) for i in snpHeader])
  #          print 'Number of snps ignored',sum([not snpDict.has_key(i) for i in genotypeDict.keys() ])
            
   #         print thisLine[0:10]

            thisLine=[genotypeDict[i] for i in snpHeader]
            
            foutcsv.writerow(thisLine)
            AnimalID=tmp[1]
            genotypeDict=deepcopy(baseGenotypes)
        genotypeDict[tmp[0]]='%s%s'%(tmp[postions[0]],tmp[postions[1]])
        
        AnimalID=tmp[1]
        
        counter+=1       
 
     #   if counter%500==0: print 'upto',counter
        tmp=split(line)
#        recodeDictTmp=dict(recodeDict)
#        recodeDictTmp[tmp[0]]=tmp[0]
#        foutcsv.writerow([recodeDictTmp[i] for i in tmp])        
    genotypeDict['ID']=AnimalID
    print 'Write animal out'
    thisLine=[genotypeDict[i] for i in snpHeader]
    foutcsv.writerow(thisLine)
    f.close()
 
    
    return 
def convertBEAGLEtoCSV(genotypeFiles,mapfile,outputfile,recodeDict):
    fout=open(outputfile,'wb')
    foutcsv=csv.writer(fout)
    
    snpHeader=['ID']
    f=open(mapfile,'rb') 
    counter=0
    for line in f:
        counter+=1
        tmp=split(line)
        snpHeader.append(tmp[0])
    f.close()
    foutcsv.writerow(snpHeader)
    for genotypeFile in genotypeFiles:
        f=open(genotypeFile,'rb') 
        counter=0
        f.readline()
        header=f.readline().split()
        print header[0:10]
       
        while 1:
            counter+=1
            if counter%500==0: print 'upto',counter
            line1=f.readline()
            if not line1:
                break
            line2=f.readline()
            
            tmp1=split(line1)
            tmp2=split(line2)
            thisGenotype=[(header[i],'%s'%recodeDict['%s%s'%(tmp1[i],tmp2[i])]) for i in range(1,len(header))]
            thisGenotype=dict(thisGenotype)
            thisGenotype['ID']=tmp1[0]
            if not thisGenotype.has_key('BovineHD0100000026'):
                thisGenotype['BovineHD0100000026']='--'
            foutcsv.writerow([thisGenotype[j] for j in snpHeader])   
        
        f.close()


    fout.close()
    return
def convertCRCtoCSV(genotypeFile,mapfile,outputfile,recodeDict):
    fout=open(outputfile,'wb')
    foutcsv=csv.writer(fout)
    
    snpHeader=['ID']
    f=open(mapfile,'rb') 
    counter=0
    for line in f:
        counter+=1
        tmp=split(line)
        snpHeader.append(tmp[0])
    f.close()
    foutcsv.writerow(snpHeader)
    
    f=open(genotypeFile,'rb') 
    counter=0
    for line in fileinput.input([genotypeFile]):
        counter+=1
        if counter%500==0: print 'upto',counter
        tmp=split(line)
        recodeDictTmp=dict(recodeDict)
        recodeDictTmp[tmp[0]]=tmp[0]
        foutcsv.writerow([recodeDictTmp[i] for i in tmp])        

    f.close()
#for line in fileinput.input([genotypeFile]):        
    return
def    createGenebleSNPmapFromCRC(snpMapin,snpMapout):
    fout=open(snpMapout,'wb')
    f=open(snpMapin,'rb')
    fout.write("chrom name position \n")
    for i in f:
        tmp=i.split()     
        fout.write("%s %s %s \n"%(tmp[2],tmp[0],tmp[3]))
    
    
if __name__ == "__main__":
    """
    os.chdir('/Users/uqmkel13/Documents/Projects/PIC/FileS1')
    genotypeFile="genotypes.txt"
    genotypesOut='genableGTY.txt'
    idlist=convertToGenable(genotypeFile,genotypesOut,{'0.00':'A/A','1.00':'A/B','2.00':'B/B'},'0/0')
    """
    os.chdir('/Users/uqmkel13/Documents/Projects/FattyAcidComposition/Testing')
    genotypeFile="testGTYPes.csv"
    import sqlite3
    os.system("rm test.db")
    conn = sqlite3.connect('test.db')
    cur=conn.cursor()
    gty=genomeFileTools(conn)
    gty.createDataStructureAndLoad(genotypeFile, cur)
    conn.create_function("recodeWombat", 1, recodeWombat)
    cur.execute("select recodeWombat(genotype) from genotypes")
    #cur.execute("select recodeWombat(genotype) from genotypes")
    x=cur.fetchall()
    for i in x:
        print i 
    
