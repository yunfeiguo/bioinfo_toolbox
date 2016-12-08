#!/usr/bin/env python
import glob
import logging
import logging.handlers
LOG_FILENAME = 'logging_rotatingfile_example.out'
logger = logging.getLogger('MyLogger')
logger.setLevel(logging.DEBUG)
handler = logging.handlers.RotatingFileHandler(LOG_FILENAME,maxBytes=20,backupCount=5)
logger.addHandler(handler)
for i in range(20):
    logger.debug('i = %d' % i)    
