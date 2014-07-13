from multiprocessing import Pool
import os

def run_command(command):
    os.system(command)

def run_all(command_queue):
    workers = Pool()
    workers.map(run_command, command_queue)
    workers.join()

def split_by_channel(filelist, targetdir):
    command_queue = []
    for filename in filelist:
        for ch in range(4):
            basename = os.path.basename(filename)[:-6]
            newname = targetdir + "/" + basename + str(ch) + ".fits"
            call_str = "./mcolor.py -c " + str(ch) + " -o " + newname + \
                " " + filename
            command_queue.append(call_str)
    run_all(command_queue)
