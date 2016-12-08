import cProfile
import pstats
cProfile.run(open('dictOpt_noClass.py','rb'),'dictOpt.stats')
dictP = pstats.Stats('dictOpt.stats')
print('dictionary optimized')
print(dictP.sort_stats('cumulative').print_stats(10))
cProfile.run(open('strOpt_noClass.py','rb'),'strOpt.stats')
strP = pstats.Stats('strOpt.stats')
print('string optimized')
print(strP.sort_stats('cumulative').print_stats(10))
