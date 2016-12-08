#!/usr/bin/env python
import sys, os, logging
import subprocess
logging.basicConfig(level = logging.INFO, format = '%(asctime)s - %(levelname)s - %(message)s')

def rsync(target, destination, node1, node2):
    #rsync target on node1 to destination on node2
    #command = "rsync -av {0}:{1} {2}:{3} ".format(node1, target, node2, destination)
    command = "scp {0}:{1} {2}:{3} ".format(node1, target, node2, destination)
    logging.info('executing {0}'.format(command))
    return subprocess.Popen(command, shell = True)

def main():
    if len(sys.argv) <= 5:
	print("Usage: {0} <target> <destination> <sourceNode> <node1,node2, ...>".format(sys.argv[0]))
	raise ValueError()
    target = sys.argv[1]
    destination = sys.argv[2]
    nodes = sys.argv[3:]
    print(nodes)
    logging.info("Broadcasting {0} on {1} to {2} on {3} using rsync.".format(target, nodes[0], destination, ','.join(nodes[1:])))
    logging.warning("[WARNING] Assume {0} is a folder".format(destination))

    #targetCopy is the target copied to destination folder
    targetCopy = os.path.join(destination, os.path.basename(target))
    end = 0 #end index of senders, inclusive
    step = 1 #end + step is receiver index, inclusive
    while step < len(nodes):
	#0...end nodes are senders
	#step...end+step nodes are receivers
	#there is a 1-to-1 relationship
	childProcess = []
	for i in range(0,end + 1):
	    if i + step >= len(nodes):
		continue
	    if i == 0:
	        childProcess.append(rsync(target, destination, nodes[i], nodes[i + step]))
	    else:
	        childProcess.append(rsync(targetCopy, destination, nodes[i], nodes[i + step]))
	#wait for all children
	for i in childProcess:
	    if i.wait() != 0:
		raise SystemError()
	end += step
	step *= 2
    logging.info("All done")

if __name__ == '__main__':
    main()
