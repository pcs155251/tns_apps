import os
#for creating directory  

#requesting resources
nselect = 1
ncpus =  8
memString = ''
#memString = ':mem=180gb'

#qName = 'cf160'
qName = 'dc20190028'
#qName = 'dc20190021'
#other queues
#qName = 'ct160'
#qName = 'ctest'
#cn0439~0443, cn0633~cn0637

#hostString = ':host=cn0443'
#hostString = ':host=cn0633'
hostString = ''

#fields   = [1, 1.2, 1.3, 1.6, 1.9, 2.0, 2.2, 3.4, 3.7, 4.0]
fields   = []
fieldStart = 3.01
fieldEnd   = 3.091
fieldDelta = 0.01
fieldLoad  = 3.08
prefix     = 'fub'
wvLoadFoldStr = '-wvLoadFolder fub/'
saveFoldStr = '-saveFolder ffa/'
#Dstart     = 14
#Dend       = 15 
#Ddelta     = 1
#ifload     = 0
#loadFolder     = 'evoStop_-10/'
#enLoadFolder   = 'cEn/'
#enSaveFolder   = 'cEn/'
#dataSaveFolder = 'cEn/'
scriptFolder = './scriptsTaiwania/'

def generate():
  createFolder( scriptFolder )
  appendFields( fields, fieldStart, fieldEnd, fieldDelta )
  for field in fields:
    fieldLoad = field
    runFuTaiwania( field, fieldLoad )
  os.system('qstat -u pcs155251')

def createFolder( dir ):
  try:
    if not os.path.exists(dir):
      os.makedirs(dir)
  except OSError:
    print ('Error: Creating directory. ' + dir)
 
def getSign( value ):
  if value > 0:
    return '+'
  elif value == 0:
    return '0'
  else:
    return '-'

def appendFields( fields, fieldStart, fieldEnd, fieldDelta ):
  if len(fields)==0:
    field = fieldStart
    while field<fieldEnd:
      fields.append(field)
      field = field+fieldDelta

#def runFuTaiwania( g_JoverQ, q_sign, prefix, Dstart, Dend, Ddelta, ifload ):
def runFuTaiwania( field, fieldLoad ):
  jobName  = prefix + str(abs(int(round(1000*field)))).zfill(4)
  scriptFile = jobName + '_Taiwania.sh'
  scriptObj = open( scriptFolder + scriptFile, 'w' )
  scriptObj.write('#!/bin/bash\n' )##
  scriptObj.write('#PBS -N %s\n' %(jobName) )##
  scriptObj.write('#PBS -P MST107027\n')
  scriptObj.write('#PBS -q %s\n' %(qName) )
  scriptObj.write('#PBS -l select=%i:ncpus=%i%s%s\n' %(nselect,ncpus,memString,hostString))##
  scriptObj.write('#PBS -M pcs155251@gmail.com\n')
  scriptObj.write('#PBS -m bea\n')
  scriptObj.write('#allow checking error and output at real time (home directory)\n')
  scriptObj.write('#PBS -W sandbox=PRIVATE\n')
  scriptObj.write('\n')
  scriptObj.write('cd $PBS_O_WORKDIR\n')
  scriptObj.write('cd ..\n')
  scriptObj.write('\n')
  scriptObj.write('module load intel/2017_u4\n')
  scriptObj.write('\n')
  scriptObj.write('./main.e ./inputFu.rc  -field %f -fieldLoad %f %s %s \n' %( field, fieldLoad, wvLoadFoldStr, saveFoldStr ))
  scriptObj.write('\n')
  print( scriptFile + ' has been generated' )
  scriptObj.close()
  os.chdir(scriptFolder)
  os.system('qsub ' + scriptFile )
  os.chdir('..' )

generate()

