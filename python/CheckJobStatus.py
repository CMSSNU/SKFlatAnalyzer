import os

#commands.sh  job_0_GetEffLumi.e6409440  job_0_GetEffLumi.o6409440  run.C  submitlog.log

def GetEventDone(l):
  # [SKFlatNtuple::Loop RUNNING] 1185000/1207307 (98.1523 %) @ 2018-5-16 11:0:9
  w = l.split()[2]
  nums = w.split('/')

  if len(nums)<2:
    print nums
    return "0:1"

  return str(nums[0])+':'+str(nums[1])

def GetJobID(logfiledir, cycle, jobnumber, hostname):

  ## SNU : Your job 7223628 ("job_0_GetEffLumi") has been submitted
  ## KNU : 3145702.cluster118.knu.ac.kr

  jobid = ""

  return jobid

def GetLogLastLine(lines):
  index = -1
  for i in range(0,len(lines)):
    l = lines[ len(lines)-1-i ]
    if "LHAPDF" in l:
      continue
    elif "lhapdf" in l:
      continue
    elif "Eur.Phys.J." in l:
      continue
    else:
      return l



def CheckJobStatus(logfiledir, cycle, jobnumber, hostname):
  FinishString = "JOB FINISHED"

  path_log_e = ""
  path_log_o = ""

  if hostname=="KISTI" or hostname=="TAMSA1" or hostname=="TAMSA2" or hostname=="KNU":
    path_log_e = logfiledir+"/job_"+str(jobnumber)+".err"
    path_log_o = logfiledir+"/job_"+str(jobnumber)+".log"

  if (not os.path.exists(path_log_e)) or (not os.path.exists(path_log_o)):
    return "BATCH JOB NOT STARTED"

  log_e = open(path_log_e).readlines()
  length_log_e = 0
  is_not_mounting_err = False
  for e_l in log_e:
    if "**** Following environment variables are going to be unset." in e_l:
        continue
    elif "PROJECT_MULTIARCH_TARGET" in e_l:
        length_log_e -= 1
    elif "WARNING: Not mounting" in e_l:
      length_log_e -= 1
      is_not_mounting_err = True
    else:
      length_log_e += 1
    
  if length_log_e > 0:
    out = 'ERROR\n'
    out += '--------------------------------------\n'
    out += 'logfile : '+path_log_o+'\n'
    out += '--------------------------------------\n'
    for l in log_e:
      out = out+l
    return out

  log_o = open(path_log_o).readlines()

  ## XX.oXX not created
  if len(log_o) == 0:
    return "BATCH LOG NOT CREATED"

  IsCycleRan = False
  for l in log_o:
    if "Processing " in l:
      IsCycleRan = True
      break

  ## XX.oXX exists

  ## "Processing run.C" not yet done
  if not IsCycleRan:
    return "ANALYZER NOT STARTED"

  LASTLINE = GetLogLastLine( log_o )
  if "Processing run.C" in LASTLINE:
    return "EVENT NOT STARTED"

  if "Event Loop Started" in LASTLINE:
    return "EVENT NOT STARTED"

  line_JobStart = ""
  for l in log_o:
    # [SKFlatNtuple::Loop] Event Loop Started 2018-05-17 19:51:10
    if "Event Loop Started" in l:
      line_JobStart = l.replace("[SKFlatNtuple::Loop] Event Loop Started ","")
      break
  ForTimeEst = LASTLINE

  ## 2) Job Finished
  if FinishString in LASTLINE:

    for i in range(0,len(log_o)):
      l = log_o[len(log_o)-1-i]
      if "[SKFlatNtuple::Loop RUNNING]" in l:
        ForTimeEst = l
        break

    # [SKFlatNtuple::Loop] JOB FINISHED 2018-12-06 04:10:37
    line_JobFinished = LASTLINE.replace("[SKFlatNtuple::~SKFlatNtuple] JOB FINISHED ","")
    EventDone = GetEventDone(ForTimeEst)
    return "FINISHED"+"\tEVDONE:"+EventDone+"\t"+line_JobStart+"\t"+line_JobFinished

  ## 4) [SKFlatNtuple::Loop] Event Loop Started 2018-06-04 18:37:57
  elif "Event Loop Started" in LASTLINE:
    return "RUNNING\t"+str(0)+"\tEVDONE:"+str(0)+"\t"+line_JobStart

  ## 3) Running
  elif "[SKFlatNtuple::Loop RUNNING]" in LASTLINE:
    # [SKFlatNtuple::Loop RUNNING] 2011000/38777460 (5.186 %)
    perct =  LASTLINE.split()[3].strip('(')
    EventDone = GetEventDone(ForTimeEst)
    return "RUNNING\t"+perct+"\tEVDONE:"+EventDone+"\t"+line_JobStart
  else:

    for it_l in range(0,len(log_o)):
      l = log_o[len(log_o)-1-it_l]
      if ("[SKFlatNtuple::Loop RUNNING]" in l) and ("@" in l):
        perct =  l.split()[3].strip('(')
        EventDone = GetEventDone(l)
        return "RUNNING\t"+perct+"\tEVDONE:"+EventDone+"\t"+line_JobStart

      return LASTLINE

