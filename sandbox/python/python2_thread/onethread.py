import cProfile
def onethread(count):
    while count > 0:
	count -= 1
cProfile.run('onethread(1e8)')
