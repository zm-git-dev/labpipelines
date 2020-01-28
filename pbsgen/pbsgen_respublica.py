#!/usr/bin/env python

from pbsgen_uge import pbsgen_main, Default

def set_queue(job):
    job.queue = 'all.q'

default = Default()
default.pbsdir = '/mnt/isilon/zhoulab/tmp/pbs'
default.scriptdir = default.pbsdir
default.stdoutdir = default.pbsdir
default.stderrdir = default.pbsdir
default.hour = 12
default.ppn = 1
default.memG = 5

pbsgen_main(default, set_queue)
