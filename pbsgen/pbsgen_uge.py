#!/usr/bin/env python

import subprocess
import sys
import os
import argparse
import types

template="""
#!/bin/bash
###
#$ -S /bin/bash
#$ -N {self.jobname}
#$ -e {self.stderr}
#$ -o {self.stdout}
#$ -q {self.queue}
#$ -l m_mem_free={self.memG}G
#$ -l h_vmem={self.memG}G
#$ -V
#$ -wd {self.workd}
#$ -pe smp {self.ppn}

### time is not used in UGE
###$ -l time={self.time}:00:00

{self.depend}

set -xe

{self.commands}
        
"""

class Job():
    
    
    def __init__(self, args):

        self.pbsdir     = os.path.abspath(args.pbsdir)
        self.nameroot   = args.nameroot
        self.ppn        = args.ppn
        self.hour       = args.hour
        self.memG       = args.memG
        #self.stderrdir  = args.stderr
        #self.stdoutdir  = args.stdout
        self.queue      = args.queue
        self.workd      = args.workd
        if args.depend:
            self.depend = '#$ -hold_jid %s' % args.depend
        else:
            self.depend = ''

    def next_index_auto(self):
        batch_index = 0
        while (True):
            if not os.path.exists('%s/j%d_%s.pbs' %
                                  (self.pbsdir, batch_index, self.nameroot)):
                return batch_index
            batch_index += 1

    def gen(self, batch_index, jobname=None):

        self.time = "%d:00:00" % self.hour

        if batch_index < 0:
            batch_index = self.next_index_auto()

        if jobname is None: ## auto generate at PBSDIR
            self.jobname = ('j%d_'+self.nameroot) % batch_index
            pbs_filename = '%s/%s.pbs' % (self.pbsdir, self.jobname)
        elif jobname[0] == '/': ## jobname is an absolute path
            self.jobname = os.path.basename(jobname)
            pbs_filename = jobname
        else: ## jobname is just a name
            self.jobname = jobname
            pbs_filename = '%s/%s.pbs' % (self.pbsdir, self.jobname)

        self.stdout = '%s.stdout' % (pbs_filename,)
        self.stderr = '%s.stderr' % (pbs_filename,)

        with open(pbs_filename, 'w') as fout:
            fout.write(template.format(self=self))
            
        return pbs_filename

class Indices:

    def __init__(self):
        self.spans = []

    def extend(self, start, end):
        self.spans.append((start, end))

    def extract(self, lst):
        result = []
        for start, end in self.spans:
            if not end:
                end = len(lst)
            result.extend([lst[_] for _ in xrange(start, end)])

        return result

def parse_indices(indstr):
    rgs = indstr.split(',')
    indices = Indices()
    for rg in rgs:
        if rg.find('-') >= 0:
            pair = rg.split('-')
            if not pair[0]:
                pair[0] = 0
            if not pair[1]:
                pair[1] = None
            indices.extend(int(pair[0])-1 if pair[0] else 0,
                           int(pair[1]) if pair[1] else None)
        else:
            indices.extend(int(rg)-1, int(rg))

    return indices

def main(args):

    job = Job(args)
    job.commands = "echo `date`\n\n"

    ## data from STDIN
    if args.command == '-':
        args.command = sys.stdin.read() 

    job.commands += args.command
    job.commands += "\n\necho `date` Done.\n"
    pbs_filename=job.gen(int(args.index), jobname=args.name)
    if not args.silent:
        sys.stderr.write('==================\n\n')
        sys.stderr.write('pbs path: %s\n' % pbs_filename)
        sys.stderr.write('pbs stdout: %s\n' % job.stdout)
        sys.stderr.write('pbs stderr: %s\n' % job.stderr)
        # sys.stderr.write('hour: %d\n' % job.hour)
        sys.stderr.write('memG: %d\n' % job.memG)
        sys.stderr.write('ppn: %d\n' % job.ppn)
        sys.stderr.write('==================\n\n')
    if args.submit:
        subprocess.check_call(['qsub', '-terse', pbs_filename])

def add_default_settings(parser, d):

    """ default setting """

    default_pbsdir = os.environ['PBSDIR'] if 'PBSDIR' in os.environ else os.getcwd()
    parser.add_argument("-pbsdir", default=default_pbsdir,
                        help="pbs directory [%s]" % default_pbsdir)
    parser.add_argument("-hour", default=d.hour, type=int,
                        help="wall time in hour [%d]" % d.hour)
    # parser.add_argument("-stdout", default=d.stdoutdir,
    #                     help="stdout [%s]" % d.stdoutdir)
    # parser.add_argument('-stderr', default=d.stderrdir,
    #                     help='stderr [%s]' % d.stderrdir)
    parser.add_argument("-ppn", default=d.ppn, type=int,
                        help="ppn [%d]" % d.ppn)
    parser.add_argument("-workd", default=os.getcwd(),
                        help="workd [%s]" % d.workd)
    default_nameroot = os.environ['NAMEROOT'] if 'NAMEROOT' in os.environ else 'LabJob'
    parser.add_argument("-nameroot", default=default_nameroot,
                        help="job name root [%s]" % default_nameroot)
    parser.add_argument("-memG", default=d.memG, type=int,
                        help="memory in G [%d]" % d.memG)
    parser.add_argument('-queue', default='all.q', help='queue to force to')
    parser.add_argument('-depend', default='', help='dependency')
    parser.add_argument('-silent', action='store_true', help='suppress info print')

def pbsgen_main(parser, setting):

    parser = argparse.ArgumentParser(description="generate pbs files")

    parser.add_argument('command', nargs='?', type=str, default='-', help='command to run')
    parser.add_argument("-index", type=int, default=-1, help="index of next pbs script [inferred from existing file names]")
    # I merged -dest into -name since they are largely redundant
    # parser.add_argument('-dest', default=None, help='destination pbs file')
    parser.add_argument('-name', default=None, help='job name (if specified, I will ignore index), if given an absolute path, use the basename as job name, use the dirname as pbsdir')
    add_default_settings(parser, Default())
    parser.add_argument('-submit', action='store_true', help='submit pbs')

    args = parser.parse_args()
    main(args)
    return

class Default:

    def __init__(self):

        self.pbsdir = os.getcwd()
        # self.stdoutdir = '%s/stdout/' % self.pbsdir
        # self.stderrdir = '%s/stderr/' % self.pbsdir
        self.hour = 96
        self.ppn = 1
        self.memG = 5
        self.workd = os.getcwd()

if __name__ == "__main__":

    def set_queue_respublica(job):
        job.queue = 'all.q'

    pbsgen_main(Default(), set_queue_respublica)
