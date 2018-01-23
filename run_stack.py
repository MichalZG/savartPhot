#!/usr/bin/env python

from src.pipeline.stack import Stack

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
	parser.add_argument('work_path', help='Files to process directory')
	parser.add_argument('stack_size', help='Number of files to stack')
	parser.add_argument('output_directory', help='Output path')

	args = parser.parse_args()

	stacking = Stack(args.config, args.work_path, int(args.stack_size), args.output_directory)
	stacking.process()

