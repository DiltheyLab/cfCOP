#!/usr/bin/env python3

import os, sys
from datetime import datetime

# exchanging irregular characters in species names
def replace_characters(input_string):
    name = input_string.replace(" ", "_")
    for c in ",.-#@:":
        name = name.replace(c, "_")
    for c in "[](){}<>%&/$§´`'!^#":
        name = name.replace(c, "")
    return name

def log_time(message):
    time_stamp = datetime.now()
    time_string = time_stamp.strftime("%Y/%m/%d, %H:%M:%S")
    sys.stdout.write('%s: %s\n'%(message,time_string))
    sys.stdout.flush()
    return

