##Generates config files, and batch submission batch.src files with appropriate directories
from  optparse  import OptionParser
import sys, os

parser = OptionParser()

parser.add_option("-r", "--raw", dest="rawlist",
		                    help="list of raw files to analyze",
		                    action="store",type="string")

parser.add_option("-c","--config",dest="config",
		                    help="hlt config modified to have replacements done for raw list",
		                    action="store",type="string")

parser.add_option("-o", "--output", dest="output",
		                    help="output destination",
		                    action="store",type="string")

parser.add_option("-j", "--jobs", dest="jobs",
		                    help="number of jobs to run at once",
		                    action="store",type="string")

parser.add_option("-p", "--rmprod", dest="rmprod",
		                    help="remove the producer files generated",
		                    action="store_true",default=False)


(options, args) = parser.parse_args()

parser.print_help()

#get some names and paths
pwd = os.getenv("PWD")
output_dir = options.output
cfg_name = options.config

#OPEN THE FILES
raw = open(options.rawlist,"r")
cfg = open(cfg_name,"r")
#READ THE LINES IN
#reco_lines = reco.readlines()
raw_lines = raw.readlines()

#Count how many config files we will need
#there better not be extra lines in the lists!
#n_reco = len(reco_lines)
n_raw = len(raw_lines)

#make the raw list we will be inserting
raw_list = ""
for ii in raw_lines:
        raw_list += "'"+ii.rstrip("\n")+"',\n"

raw_list=raw_list.rstrip(",\n")+"\n"

#remove old instances of the output directory
os.system("rm -r " + output_dir)
#build the directories
os.system("mkdir " + output_dir)
os.system("mkdir " + output_dir+"/cfgs")
os.system("mkdir " + output_dir+"/logs")
os.system("mkdir " + output_dir+"/src")
os.system("mkdir " + output_dir+"/res")
os.system("cp ~/josh_scripts/build_hlt_tree.py %s" % output_dir)

#strings to be replaced in the config file
#reco_replace_string = "#$REPLACE_RECO$#"
raw_replace_string = "#$REPLACE_RAW$#"
output_replace_string = "#$REPLACE_OUTPUT$#"

commands = []


for ii in range(n_raw):
	#make the new cfg file	
        new_file_name = output_dir+"/cfgs/"+cfg_name.rstrip(".py")+"_"+str(ii)+".py"
        cfg_file = open(cfg_name, "r")
        cfg_file_lines = cfg_file.readlines()
        new_file = open(new_file_name,"a")
	output_file_path = output_dir+"/res/ohlt_output_"+str(ii)+".root" 
	
        #read all the lines of the old cfg file and replace the replace text
        for jj in range(len(cfg_file_lines)):
		#dont replace the reco for now (version is too old 53X)
                #if reco_replace_string in cfg_file_lines[jj]:
                #        cfg_file_lines[jj] = "'"+reco_lines[ii].rstrip("\n")+"'"
                if raw_replace_string in cfg_file_lines[jj]:
                        cfg_file_lines[jj] = cfg_file_lines[jj].replace(raw_replace_string,raw_lines[ii].rstrip("\n"))
                if output_replace_string in cfg_file_lines[jj]:
                        replace_string_out= output_file_path
                        cfg_file_lines[jj] = cfg_file_lines[jj].replace(output_replace_string,replace_string_out)
        #write the cfg file
        new_file.writelines(cfg_file_lines)

	#generate the commands to run on the raw and convert to trees
	output_tree_path = output_dir+"/res/hlt_tree_%i.root" % ii
	cmd_cmsrun = "cmsRun %s\n" % new_file_name

	cmd_tree = "python %s/build_hlt_tree.py -f %s -o %s \n" % (output_dir,output_file_path,output_tree_path)
	cmd_rm = "rm %s \n" % output_file_path

	commands.append((cmd_cmsrun,cmd_tree,cmd_rm))

njobs = options.jobs
job = 0
file_counter = 0
cmds_per_job = int(len(commands)) / int(njobs) #integer division
bsub_cmds = []
bsub_file = None

#loop over all commands
for ii in range(len(commands)):
	#make sure the extra jobs are added to the last file
	if (file_counter == 0) and (job != njobs):		
		#write a new job
		bsub_file_name = output_dir+"/src/bsub_%i.src" % job
		bsub_file = open(bsub_file_name,"a")
		bsub_file.write('#!/bin/bash\n')
		bsub_file.write("cd /afs/cern.ch/user/h/hardenbr/2013/HIGGS_DIPHOTON_HLT/NEW_OPEN_HLT/CMSSW_5_2_8_patch1/src\n")
		bsub_file.write("export SCRAM_ARCH=slc5_amd64_gcc462 \n")
		bsub_file.write("eval `scramv1 ru -sh`\n")

		cmd="bsub -q 1nd "
		cmd+= "-o %s/logs/job_%i.log " % (output_dir ,job)
		cmd+="source "+ bsub_file_name

		bsub_cmds.append(cmd)

		#increment the counter
		job+=1

	bsub_file.write(commands[ii][0])
	bsub_file.write(commands[ii][1])
	if options.rmprod:
		bsub_file.write(commands[ii][2])
	
	#if we are at the number of commands per job reset the counter
	file_counter+=1
	if file_counter == cmds_per_job:file_counter = 0		
	if job == njobs:
		continue

#print out the submission commands
for ii in bsub_cmds:print ii

