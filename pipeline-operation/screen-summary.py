#!/usr/bin/env python3

"""
Prints the last ten lines of all active screen sessions
"""

import os
import re
import time
import tempfile
import subprocess

COLOR_CYAN = '\x1b[0;36m'
COLOR_MAGENTA = '\x1b[0;35m'
COLOR_END = '\x1b[0m'

# Regular expression that matches the user's prompt.  If this matches the last
# line then the command that was running probably completed, and we should flag
# it.
PROMPT_PATTERN = ".* [$]$"

def list_screens():
    try:
        subprocess.check_output(["screen", "-ls"])
    except subprocess.CalledProcessError as e:
        if e.returncode == 1:
            output = e.output # expected behavior
        else:
            raise

    return [line.split()[0]
            for line in output.decode('utf-8').split("\n")
            if "(Detached)" in line]

def start():
    for screen in list_screens():
        with tempfile.TemporaryDirectory() as workdir:
            tmpfname = os.path.join(workdir, "tmp.txt")

            subprocess.check_call([
                "screen",
                "-S", screen,
                "-X", "hardcopy",
                tmpfname])
            
            # wait for screen to dump like we asked
            while not os.path.exists(tmpfname):
                time.sleep(0.01)

            with open(tmpfname) as inf:
                s = inf.read()

            # delete all trailing empty lines
            lines = re.sub("\n*$", "", s).split("\n")

            screen_label_color = COLOR_CYAN
            if lines and re.match(PROMPT_PATTERN, lines[-1]):
                screen_label_color = COLOR_MAGENTA

            print(screen_label_color + screen + ":" + COLOR_END)

            for line in lines[-10:]:
                print(line)

if __name__ == "__main__":
    start()
