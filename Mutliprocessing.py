# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 13:01:34 2017

@author: HP
"""

###############################################################
# Run all Models in Parallel using MultiThreading
###############################################################

#-------------------------------------------------------------
# Step 1: Library Inclusion                             
#-------------------------------------------------------------
import os, sys, time, threading, multiprocessing


#--------------------------------------------------------------
# Step 2: Model program file name
#--------------------------------------------------------------
# Add new folder in the list
srno=[0,1,2,3]


#--------------------------------------------------------------
# Step 3: Checking for correct number of parameters
#--------------------------------------------------------------
# Getting command line arguments
args=sys.argv
if len(args)!=1:
  print "\nError !!! Wrong number of parameters"
  #print "\nUsages: $python runAllModelsInParallel.py <dataFileName.csv> <trainingPercentage> <acceptableError>"
  #print "\nExample: $python runAllModelsInParallel.py dataFileName.csv 70 1\n"
  exit()  
  

#--------------------------------------------------------------
# Step 4: Variable Declaration
#--------------------------------------------------------------
#dataFileName=args[1]       # Data FileName
#training = int(args[2])         # Training Percentage; Testing = 100 - Training
#acceptableError = int(args[3])  # acceptableError


#-------------------------------------------------------------
# Step 5: Getting Number of Cores/CPUs
#-------------------------------------------------------------
numberOfCores=multiprocessing.cpu_count()


#--------------------------------------------------------------
# Step 6: Getting the starting time
#--------------------------------------------------------------
print "\n.........All Running Model....."
startTime = time.time()



#-------------------------------------------------------------
# Step 7: Function Defination
#-------------------------------------------------------------

def executeModel(cmd):
    os.system(cmd)
    return


#-------------------------------------------------------------
# Step 8: Start Program
#-------------------------------------------------------------

c=0
for i in srno:
    cmd='py %s %d >tmp%d'%("final.py",i,c)
    t = threading.Thread(target=executeModel , args=(cmd,))
    t.start()
    c=c+1
    print "Running ...... ", c,"/",len(srno) 
    time.sleep(1)
    while True:
        if threading.activeCount() <= numberOfCores:
            break
        time.sleep(2)
    

while True:
  if threading.activeCount() == 1:
    break
  time.sleep(5)
  print "Model Left ... ",threading.activeCount() - 1
#dele
os.system('del tmp*')
#--------------------------------------------------------------
# Step 9: Merge all result file into one file
#--------------------------------------------------------------
cmd="type *protein_dataset.csv >>final_dataset.csv"
os.system(cmd)




#--------------------------------------------------------------
# Step 10: Grand Total Running Time
#--------------------------------------------------------------
totalTime = time.time() - startTime
print "\n\nTotal Running Time:", totalTime, " sec"
print "\nFinished\n"

