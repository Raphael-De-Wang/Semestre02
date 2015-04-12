#!env python

import rpyc
import time
import pxssh
import os.path
import getpass
import subprocess
import numpy as np
from ctypes import Structure, c_char_p, c_char
from multiprocessing import Process, Value, Array, Lock

portNum = 18812

SSH_COMMAND = "ssh"
session_address = "ssh.ufr-info-p6.jussieu.fr"
machine_address = "ppti-14-%d%02d-%02d.ufr-info-p6.jussieu.fr"

class Master(object):
    def __init__(self,slaveList,uname):
        self.uname    = uname
        self.poolSize = 0
        self.slaveList= slaveList
        self.connPool= []
        
    def initConnPool(self,RESET=True):
        if RESET:
            self.shutdown()
        self.connPool = [ Process(target=self.job) for i in xrange(len(self.slaveAliveList))]
        
    def startup(self):
        for slave in self.connPool:
            slave.start()
            
    def shutdown(self):
        for conn in self.connPool:
            if conn.is_alive():
                conn.terminate()
                
    def join(self):
        for slave in self.connPool:
            if slave.is_alive():
                slave.join()
                
    def isAlive(self,slaveAddr):
        return 0 == subprocess.call(["ping", "-c", "1", "-W", "1", slaveAddr]) 

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
        self.pwd = getpass.getpass("Password for Slave Connections : ")
        for slaveAddress in self.slaveList:
            if self.isAlive(slaveAddress):
                slaveList.append(slaveAddress)
                try:
                    self.slaveConnList.append(self.slaveConnect(slaveAddress))
                except:
                    print "Failed to connect ", slaveAddress
                    continue
        self.slaveAliveList = Array(c_char_p, slaveList)
                    
    def slaveDisconnect(self,conn):
        if conn.isalive():
            conn.logout()
            print ">> Connection %s is logout.<< "%conn.args[-1]

    def slaveDisconnectAll(self):
        for conn in self.slaveConnList:
            self.slaveDisconnect(conn)
        
    def slaveServiceUp(self,conn):
        conn.sendline("python "+self.slaveCMDPath)
        
    def slaveServiceAllUp(self):
        for conn in self.slaveConnList:
            self.slaveServiceUp(conn)
        
    def slaveServiceDown(self,slaveConn):
        try: 
            slaveConn.sendline("killall python")
        except:
            print "kill service failed, probably service is not up. "
            pass
        
    def slaveServiceAllDown(self):
        for conn in self.slaveConnList:
            self.slaveServiceDown(conn)
            
    def connect_slave(self,slaveID):
        try:
            conn = rpyc.connect(self.slaveAliveList[slaveID],self.port.value)
        except:
            print "process [%d] first time fail"%slaveID
            time.sleep(5)
            try:
                conn = rpyc.connect(self.slaveAliveList[slaveID],self.port.value)
            except:
                print "process [%d] second time fail"%slaveID
                raise
        print conn.root.slaveHostname(), ' is connected. '
        return conn
    
    def disconnect_slave(self,conn):
        if not conn.closed:
            print "Slave Service Connection with ", conn.root.slaveHostname(), ' is disconnecting. '
            conn.close()
        else:
            print "Slave Service Connection with ", conn.root.slaveHostname(), ' is disconnected. '
        
    def preJob(self):
        raise NotImplementedError("preJob non implemente.")
    
    def job(self,jobArgs):
        raise NotImplementedError("job non implemente.")
    
    def persistResultat(self):
        raise NotImplementedError("persistResultat non implemente.")

    def reset(self):
        print 'Connect to slaves...'
        self.slaveConnectAll()
        print 'Shutdown all slaves services.'
        self.slaveServiceAllDown()
        print 'Shutdown all process in connection pool. '
        self.shutdown()
        print 'Disconnect from slaves.'
        self.slaveDisconnectAll()
        
    def run(self):
        try:
            print 'Connect to slaves...'
            self.slaveConnectAll()
            print 'start up slave services...'
            self.slaveServiceAllUp()
            print 'init Connection pool...'
            self.initConnPool()
            print 'preparing Jobs...'
            self.preJob()
            print 'starting...'
            self.startup()
            print 'joining...'
            self.join()
            print 'persisting resultat.'
            self.persistResultat()
            print 'Shutdown all slaves services.'
            self.slaveServiceAllDown()
            print 'Disconnect from slaves...'
            self.slaveDisconnectAll()
            print 'Shutdown all process in connection pool. '
            self.shutdown()
        except:
            print 'Shutdown all slaves services.'
            self.slaveServiceAllDown()
            print 'Shutdown all process in connection pool. '
            self.shutdown()
            print 'Disconnect from slaves.'
            self.slaveDisconnectAll()
            raise
        
