#!env python

from BioDist import *

def tacheName(target,tacheNum):
    return target+"-job%09d"%(tacheNum)

def rstName(target,tacheNum):
    return target+"-rst%09d"%(tacheNum)

def slaveClaissic():
    subprocess.check_output(["rpyc_classic"])

class exprLevelJob(Master):
    def __init__(self,slaveList,uname,queryFile,tacheSize,dbName,dbVolPath,doDivision=False, doInitDB=False, dbDataFile=None, dbTitle=None):
        super(exprLevelJob,self).__init__(slaveList,uname)
        self.target = queryFile
        self.tacheSize= tacheSize
        self.dbName = dbName
        # self.slaveList = Array(c_char_p, self.slaveList)
        self.doDivision = doDivision
        self.doInitDB = doInitDB
        self.dataFile = dbDataFile
        self.dbTitle  = dbTitle
        # self.slaveID = Value('i',0,lock=False)
        # self.port   = Value('i',portNum,lock=False)
        self.slaveID  = 0
        self.port     = portNum
        self.lock   = threading.Lock()
        self.slaveCMDPath = '/users/nfs/Etu9/3404759/Workspace/Semestre02/BimProjet/Part2/Script/slave.py'
        # self.dbVolumePath = Array('c',dbVolPath)
        self.dbVolumePath = dbVolPath

    def initBlastDB(self):
        '''
        makeblastdb [-h] [-help] [-in input_file] [-input_type type]
        -dbtype molecule_type [-title database_title] [-parse_seqids]
        [-hash_index] [-mask_data mask_data_files] [-mask_id mask_algo_ids]
        [-mask_desc mask_algo_descriptions] [-gi_mask]
        [-gi_mask_name gi_based_mask_names] [-out database_name]
        [-max_file_sz number_of_bytes] [-logfile File_Name] [-taxid TaxID]
        [-taxid_map TaxIDMapFile] [-version]
        makeblastdb -in ERR420359_62_MERGED_FASTQ_readsWithpCDS.fasta -title RNA_Marker113 -out RNA_Marker113 -dbtype nucl
        '''
        return subprocess.check_output(["makeblastdb", "-in", self.dataFile, "-title", self.dbTitle, "-out", os.path.join(os.path.dirname(self.dataFile),self.dbName), "-dbtype", "nucl"])
                
    def divideQueryFile(self):
        iCount    = 0
        iFlag     = False
        tgrHandle = open(self.target,'r')
        jobHandle = open(tacheName(self.target,0),'w')
        for line in tgrHandle:
            if line[0] == '>':
                iCount += 1
                if not iFlag:
                    iFlag = True
            if iCount % self.tacheSize == 0 and iFlag:
                jobHandle.close()
                jobHandle = open(tacheName(self.target,iCount/self.tacheSize),'w')
                iFlag = False
            jobHandle.write(line)
        jobHandle.close()
        return int(np.ceil(1.0*iCount/self.tacheSize))

    def getTacheNum(self):
        iCount = -1
        while os.path.exists(tacheName(self.target,iCount+1)):
            iCount += 1
        if iCount < 0:
            raise ValueError("No Tache Exists.")
        return iCount + 1
    
    def preJob(self):
        if self.doDivision:
            print "Dividing Query File..."
            self.tacheNum = self.divideQueryFile()
        else:
            self.tacheNum = self.getTacheNum()
        if self.doInitDB:
            print "init Blast DB..."
            self.initBlastDB()
        self.currTaNum= 0
        self.tacheLock= Lock()

    def getTacheID(self):
        self.tacheLock.acquire()
        tacheID = self.currTaNum.value
        self.currTaNum.value += 1
        self.tacheLock.release()
        return tacheID
        
    def blastHits(self):
        self.lock.acquire()
        slaveID = self.slaveID.value
        self.slaveID.value += 1
        self.lock.release()
        conn = self.connect_slave(slaveID)
        STOP = False
        while True:
            tacheID = self.getTacheID()
            if tacheID >= self.tacheNum.value:
                break
            print "tacheID: ", tacheID
            if 0 <> conn.root.blastHits(self.dbName,tacheName(self.target,tacheID),rstName(self.target,tacheID),self.dbVolumePath.value):
                print "blast tache [%d] is failed. "%tacheID
        self.disconnect_slave(conn)

    def loadHits(self,conn):
        STOP = False
        while True:
            tacheID = self.getTacheID()
            if tacheID >= self.tacheNum.value:
                break
            print "tacheID: ", tacheID
            if 0 <> conn.root.hitsToNivExpr(rstName(self.target,tacheID)):
                print "load Hits tache [%d] is failed. "%tacheID

    def sortHits(self,conn):
        conn.root.sortNivExprDict()
                
    def job2(self):
        self.lock.acquire()
        slaveID = self.slaveID.value
        self.slaveID.value += 1
        self.lock.release()
        conn = self.connect_slave(slaveID)
        self.loadHits(conn)
        self.sortHits(conn)
        self.disconnect_slave(conn)
        
    def job(self):
        self.blastHits()
        
    def persistResultat(self):
        pass

slaveList = [machine_address%(5,j,i) for i in range(2,17) for j in [2,3] ]
uname     = '3404759'
tacheSize = 10000

# queryFile = '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS828_DNA_N3/ERR420353_54_MERGED_FASTQ_readsWithpCDS.fasta'
# dbDataFile= '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS828_RNA_N3/ERR420349_52_MERGED_FASTQ_readsWithpCDS.fasta'
# dbVolPath = '/Vrac/WJX3404759/BioDB/Axial_FS828_RNA_N3/'

# queryFile = '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS821_DNA_Marshmallow/ERR420344_48_MERGED_FASTQ_readsWithpCDS.fasta'
# dbDataFile= '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS821_RNA_Marshmallow/ERR420342_43_MERGED_FASTQ_readsWithpCDS.fasta'
# dbVolPath = '/Vrac/WJX3404759/BioDB/Axial_FS821_RNA_Marshmallow/'

queryFile = '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS861_DNA_Marker113/ERR420363_68_MERGED_FASTQ_readsWithpCDS.fasta'
dbDataFile= '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS861_RNA_Marker113/ERR420359_62_MERGED_FASTQ_readsWithpCDS.fasta'
dbVolPath = '/Vrac/WJX3404759/BioDB/Axial_FS861_RNA_Marker113'

dbName    = 'FS861_RNA_Marker113' # "FS821_RNA_Marshmallow" # "FS828_RNA_N3" # "FS857_RNA_Anemone"
dbTitle   = 'FS861_RNA_Marker113' # "FS821_RNA_Marshmallow" # "FS828_RNA_N3" # "FS857_RNA_Anemone"

doDivision= True
doInitDB  = True

job = exprLevelJob(slaveList,
                   uname,
                   queryFile,
                   tacheSize,
                   dbName,
                   dbVolPath,
                   doDivision,
                   doInitDB,
                   dbDataFile,
                   dbTitle)

# job.reset()
job.run()
# job.initBlastDB()

# "/users/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS857_DNA_Anemone/ERR420357_58_MERGED_FASTQ_readsWithpCDS.fasta"
# '/Vrac/WJX3404759/BioDB/Axial_FS857_RNA_Anemone/'
    
