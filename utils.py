from multiprocessing import Pool
import os

def run_command(command):
    os.system(command)

def run_all(command_queue):
    workers = Pool()
    workers.map(run_command, command_queue)
    workers.close()
    workers.join()

def split_by_channel(filelist, targetdir, parallel=True):
    command_queue = []
    for filename in filelist:
        for ch in range(4):
            basename = os.path.basename(filename)[:-6]
            newname = targetdir + "/" + basename + str(ch) + ".fits"
            call_str = "./mcolor.py -c " + str(ch) + " -o " + newname + \
                " " + filename
            command_queue.append(call_str)
    if parallel:
        run_all(command_queue)
    else:
        for command in command_queue:
            os.system(command)
