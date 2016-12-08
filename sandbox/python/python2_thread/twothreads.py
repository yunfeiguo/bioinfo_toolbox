import cProfile
import threading
def twothreads(count):
    while count > 0:
	count -= 1
def run():
    t1 = threading.Thread(target=twothreads,args=(5e7,))
    t2 = threading.Thread(target=twothreads,args=(5e7,))
    t1.start()
    t2.start()
    t1.join()
    t2.join()
cProfile.run('run()')    
