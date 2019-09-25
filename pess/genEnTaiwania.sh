import os
#for creating directory  

#requesting resources
nselect = 1
ncpus = 40
#qName = 'dc20190028'
qName = 'dc20190021'
#qName = 'cf160'
#qName = 'ct160'
#qName = 'ctest'

mems = 180

#cn0439~0443, cn0633~cn0637
#hostString = ':host=cn0443'
#hostString = ':host=cn0633'
hostString = ''

'''
fields   = [0.19]
fieldStart = 0.19
fieldEnd   = 0.201
fieldDelta = 0.01
q_sign     = 1.0
'''

prefix     = 'test'
saveFolder = 'test'
ifsvd      = 1
enviFolder = saveFolder
Dstart     =  2
Dend       =  5 
scriptFolder = './scriptsTaiwania/'

def generate():
  createFolder( scriptFolder )
  runEnTaiwania( prefix, Dstart, Dend )
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

def runEnTaiwania( prefix, Dstart, Dend ):
  jobName  = prefix + 's' + str(Dstart) + 'e' + str(Dend)
  scriptFile =  jobName + '_Taiwania.sh'
  scriptObj = open( scriptFolder + scriptFile, 'w' )
  scriptObj.write('#!/bin/bash\n' )##
  scriptObj.write('#PBS -N %s\n' %(jobName) )##
  scriptObj.write('#PBS -P MST107027\n')
  scriptObj.write('#PBS -q %s\n' %(qName) )
  scriptObj.write('#PBS -l select=%i:ncpus=%i:mem=%igb%s\n' %(nselect,ncpus,mems,hostString))##
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
  scriptObj.write('./en.e input.rc  -dimDstart %i -dimDend %i -saveFolder %s -enviFolder %s -ifsvd %i\n'\
                                   %(   Dstart,       Dend,    saveFolder,    enviFolder,    ifsvd )
                 )
  scriptObj.write('\n')
  print( scriptFile + ' has been generated' )
  scriptObj.close()
  os.chdir(scriptFolder )
  os.system('qsub ' + scriptFile )
  os.chdir('..' )

generate()

