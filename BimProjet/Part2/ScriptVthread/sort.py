#!env python


import rpyc
import time
import pxssh
import os.path
import getpass
import threading
import subprocess
import numpy as np

portNum = 18812

SSH_COMMAND = "ssh"
session_address = "ssh.ufr-info-p6.jussieu.fr"
machine_address = "ppti-14-%d%02d-%02d.ufr-info-p6.jussieu.fr"
dbPath      = "/users/Etu9/3404759/Workspace/Sources/BioDB/"

class Master(object):
    def __init__(self,slaveList,uname,port,queryFile,outputFileName):
        self.uname    = uname
        self.poolSize = 0
        self.slaveList= slaveList
        self.connPool = []
        self.port     = port
        self.dbPath   = dbPath
        self.queryFile= queryFile
        self.outputFileName = outputFileName
        self.slaveCMDPath = '/users/nfs/Etu9/3404759/Workspace/Semestre02/BimProjet/Part2/ScriptVthread/slave.py'
                
    def isAlive(self,slaveAddr):
        try:
            rst = subprocess.check_output(["ping", "-c", "1", "-W", "1", slaveAddr])
        except:
            return False
        return True

    def slaveConnect(self,slaveAddress):
        conn = pxssh.pxssh()
        if conn.login(slaveAddress,self.uname,self.pwd):
            print ">> Connection %s is login. << "%slaveAddress
            return conn
        else:
            raise ValueError("Password Refuse by %s"%slaveAddress)

    def slaveConnectAll(self):
        self.slaveConnList = []
        slaveList = []
        self.pwd = 'Lu@81823334' # getpass.getpass("Password for Slave Connections : ")
        for slaveAddress in self.slaveList:
            if self.isAlive(slaveAddress):
                slaveList.append(slaveAddress)
                try:
                    self.slaveConnList.append(self.slaveConnect(slaveAddress))
                except:
                    print "Failed to connect ", slaveAddress
                    continue
        self.slaveAliveList = slaveList
                    
    def slaveDisconnect(self,conn):
        if conn.isalive():
            conn.logout()
            print ">> Connection %s is logout.<< "%conn.args[-1]

    def slaveDisconnectAll(self):
        for conn in self.slaveConnList:
            self.slaveDisconnect(conn)
        
    def slaveServiceUp(self,conn):
        conn.sendline("python " + self.slaveCMDPath + ' &')
        
    def slaveServiceAllUp(self):
        for conn in self.slaveConnList:
            self.slaveServiceUp(conn)
        
    def slaveServiceDown(self,slaveConn):
        try: 
            slaveConn.sendline("killall python &")
        except:
            print "kill service failed, probably service is not up. "
        
    def slaveServiceAllDown(self):
        for conn in self.slaveConnList:
            self.slaveServiceDown(conn)
            
    def connect_slave(self,slaveID):
        try:
            conn = rpyc.connect(self.slaveAliveList[slaveID],self.port)
        except:
            print "process [%d] first time fail"%slaveID
            time.sleep(5)
            try:
                conn = rpyc.connect(self.slaveAliveList[slaveID],self.port)
            except:
                print "process [%d] second time fail"%slaveID
                raise
        print "slave service on ", conn.root.slaveHostname(), ' is connected. '
        return conn

    def connectAllServices(self):
        self.services = []
        for slaveID in range(len(self.slaveAliveList)):
            serv = self.connect_slave(slaveID)
            self.services.append(serv)
            
    def disconnect_slave(self,conn):
        if not conn.closed:
            print "Slave Service on ", conn.root.slaveHostname(), ' is disconnecting. '
            conn.close()
        else:
            print "Slave Service on ", conn.root.slaveHostname(), ' is disconnected. '
            
    def disconnectAllServices(self):
        for conn in self.services:
            self.disconnect_slave(conn)
            
    def reset(self):
        print 'Connect to slaves...'
        self.slaveConnectAll()
        print 'Shutdown all slaves services.'
        self.slaveServiceAllDown()
        print 'Disconnect from slaves.'
        self.slaveDisconnectAll()

    def job(self,sid):
        serv = self.services[sid]
        while True:
            self.tacheIDlock.acquire()
            tid = self.tacheID
            self.tacheID += 1
            self.tacheIDlock.release()
            if tid >= self.tacheNum:
                break
            serv.root.hitsToNivExpr(self.rstName(self.queryFile,tid))
        serv.root.sortNivExprDict()
        self.servUpdateNGList(sid)

    def rstName(self,queryFile,tacheID):
        return queryFile+"-rst%03d"%(tacheID)
        
    def getTacheNum(self):
        iCount = -1
        while os.path.exists(self.rstName(self.queryFile,iCount+1)):
            iCount += 1
        if iCount < 0:
            raise ValueError("No Tache Exists.")
        return iCount + 1

    def servUpdateNGList(self,sid):
        serv = self.services[sid]
        rst = serv.root.getNivExpr()
        serv.root.indiceAddOne()
        self.gNameList[sid]   = rst[0]
        self.nivExprList[sid] = rst[1]
            
    def startup(self):
        self.jobOnServList = []
        self.nivExprList   = [0  for i in xrange(len(self.services))]
        self.gNameList     = ['' for i in xrange(len(self.services))]
        self.tacheID       = 0
        self.tacheIDlock   = threading.Lock()
        self.tacheNum      = self.getTacheNum()
        print "getTacheNum(): ", self.getTacheNum()
        
        for sid in range(len(self.services)):
            self.jobOnServList.append(threading.Thread(target=self.job,args=[sid]))

        for job in self.jobOnServList:
            job.start()
            
        for job in self.jobOnServList:
            job.join()

        handle = open(self.outputFileName,'w')
        count = 0
        while sum(self.nivExprList) > 0:
            sid = np.argmax(self.nivExprList)
            niv = self.nivExprList[sid]
            gN  = self.gNameList[sid]
            count += niv
            handle.write("%s\t%d\n"%(gN,niv))
            self.servUpdateNGList(sid)
            '''
            neMax = max(self.nivExprList)
            for sid,value in enumerate(self.nivExprList):
                if value == neMax:
                    niv = self.nivExprList[sid]
                    gN  = self.gNameList[sid]
                    count += niv
                    handle.write("%s\t%d\n"%(gN,niv))
                    self.servUpdateNGList(sid)
            '''
                    
        handle.close()

        print "%d records are writed..."%count
        
    def run(self):
        try:
            print 'Connect to slaves...'
            self.slaveConnectAll()
            print 'start up slave services...'
            self.slaveServiceAllUp()
            print 'connect to all services...' 
            self.connectAllServices()
            
            if len(self.services) > 0:
                print 'starting...'
                self.startup()
            else:
                print 'no salve ready.'
                
            print 'disconnect from all services...' 
            self.disconnectAllServices()
            print 'Shutdown all slaves services.'
            self.slaveServiceAllDown()
            print 'Disconnect from slaves...'
            self.slaveDisconnectAll()
        except:
            print 'Shutdown all slaves services.'
            self.slaveServiceAllDown()
            print 'Disconnect from slaves.'
            self.slaveDisconnectAll()
            raise
        

slaveList = [machine_address%(5,j,i) for i in range(2,17) for j in [3] ]
uname     = '3404759'

queryFile1 = '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS828_DNA_N3/ERR420353_54_MERGED_FASTQ_readsWithpCDS.fasta'
queryFile2 = '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS857_DNA_Anemone/ERR420357_58_MERGED_FASTQ_readsWithpCDS.fasta'
# queryFile1 = '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS821_DNA_Marshmallow/ERR420344_48_MERGED_FASTQ_readsWithpCDS.fasta'
# queryFile2 = '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS861_DNA_Marker113/ERR420363_68_MERGED_FASTQ_readsWithpCDS.fasta'

outputFileName1 = '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS828_DNA_N3/ERR420353_54_MERGED_FASTQ_readsWithpCDS.fasta-nivExpr'
outputFileName2 = '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS857_DNA_Anemone/ERR420357_58_MERGED_FASTQ_readsWithpCDS.fasta-nivExpr'
# outputFileName1 = '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS821_DNA_Marshmallow/ERR420344_48_MERGED_FASTQ_readsWithpCDS.fasta-nivExpr'
# outputFileName2 = '/users/nfs/Etu9/3404759/Workspace/Sources/BioDB/Axial_FS861_DNA_Marker113/ERR420363_68_MERGED_FASTQ_readsWithpCDS.fasta-nivExpr'

m = Master(slaveList,uname,portNum,queryFile1,outputFileName1)
m.run()

m = Master(slaveList,uname,portNum,queryFile2,outputFileName2)
m.run()
