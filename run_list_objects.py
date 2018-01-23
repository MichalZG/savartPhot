#!/usr/bin/env python

from src.pipeline.list_objects import List_objects

from glob import glob
import argparse
from astropy.io import fits
import os
import sys

class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


if __name__ == '__main__':
	parser = MyParser()
	parser.add_argument('config', help='Configuration file path')
	parser.add_argument('work_path', help='Files to process path')
	args = parser.parse_args()

	listing = List_objects(args.config, args.work_path)
	listing.show_objects()
	listing.create_stars_list()