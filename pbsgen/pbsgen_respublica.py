#!/usr/bin/env python

from pbsgen_uge import pbsgen_main, Default

def set_queue(job):
    job.queue = 'all.q'

default = Default()
default.pbsdir = '/mnt/isilon/zhoulab/tmp/pbs'
default.scriptdir = '%s/pbs/' % default.pbsdir
default.stdoutdir = '%s/stdout/' % default.pbsdir
default.stderrdir = '%s/stderr/' % default.pbsdir
default.hour = 12
default.ppn = 1
default.memG = 2

pbsgen_main(default, set_queue)
