#!/usr/bin/python
# used by Display.sh to get the list of directory sorted
import sys

dir_list_pre = sys.argv[1].split()
print ' '.join(sorted(dir_list_pre, key=lambda x: float(x))),
