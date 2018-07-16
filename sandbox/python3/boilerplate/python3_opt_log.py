#!/usr/bin/env python

import argparse
import logging
import os
import sys
import re
logger = None

def process_args(args):
    print(args.password_file)
    return
if __name__ == "__main__":
    INFO_FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    DEBUG_FORMAT = '%(levelname)s %(asctime)-15s %(name)-15s %(funcName)-20s %(message)s'

    parser = argparse.ArgumentParser(description="get data for processing", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("password_file", type = argparse.FileType('r'))
    parser.add_argument("--prefix", type = str, default = 'output', help = 'output prefix')
    parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')
    parser.add_argument("--verbose","-v", action = 'count',  help='increase verbosity')
    args = parser.parse_args()

    if (args.verbose is not None) and (args.verbose >= 1):
        logging.basicConfig(level=logging.DEBUG, format = DEBUG_FORMAT)
    else:
        logging.basicConfig(level=logging.INFO, format= INFO_FORMAT)
    logger = logging.getLogger(__name__)

    logger.info("working hard ...")
    process_args(args)
    logger.info("Done.")
