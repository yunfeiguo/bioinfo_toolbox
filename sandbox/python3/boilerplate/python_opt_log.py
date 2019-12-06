#!/usr/bin/env python

import argparse
import logging
import os
import sys
import re
logger = None

def my_function(blah):
    return
if __name__ == "__main__":
    INFO_FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    DEBUG_FORMAT = '%(levelname)s %(asctime)-15s %(name)-15s %(funcName)-20s %(message)s'

    parser = argparse.ArgumentParser(description="program name", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("input1", type = file)
    parser.add_argument("input2", type = file)
    parser.add_argument("--selection", type = str, default = 'a', choices = ['a', 'b', 'c'], help = 'choose from a,b,c')
    parser.add_argument("--cutoff", type = int, default = 1, help = 'cutoff score')
    parser.add_argument("--variable_args", type = float, action = 'append', nargs = 3, 
	                default = [1.0,2.0,1.2], help = '3 scores')
    parser.add_argument("--variable_args_but_no_need_for_multiple_option_name", type = float, nargs = '+', help = '3 scores')
    parser.add_argument("--log_to_console", action = 'store_true', help = 'log to console instead of a file')
    parser.add_argument("--verbose","-v", action = 'count',  help='increase verbosity')
    args = parser.parse_args()

    if args.verbose >= 1:
        if args.log_to_console:
            logging.basicConfig(level=logging.DEBUG, format = DEBUG_FORMAT)
        else:
            logging.basicConfig(level=logging.DEBUG, format = DEBUG_FORMAT, filename=args.prefix + '.log')
    else:
        if args.log_to_console:
            logging.basicConfig(level=logging.INFO, format= INFO_FORMAT)
        else:
            logging.basicConfig(level=logging.INFO, format= INFO_FORMAT, filename=args.prefix + '.log')
    logger = logging.getLogger(__name__)

    logger.info("working hard ...")
    my_function(args.input1, args.input2)
    logger.info("Done.")
