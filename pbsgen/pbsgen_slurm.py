#!/usr/bin/env python

import subprocess
import sys
import os
import argparse

template="""#!/bin/bash
###
#SBATCH -J {self.jobname}
#SBATCH -e {self.stderr}
#SBATCH -o {self.stdout}
#SBATCH --mem={self.memG}G
{self.gpu}
#SBATCH -c {self.ppn}
#SBATCH --export=ALL
#SBATCH -D {self.workd}
#SBATCH -n {self.ntasks}
#SBATCH -t {self.days}-{self.hours}:{self.minutes}:00
{self.depend}
{self.nodes}

set -xe

{self.commands}
"""

class Job():
    def __init__(self, args):
        self.pbsdir     = os.path.abspath(args.pbsdir)
        self.nameroot   = args.nameroot
        self.ppn        = args.ppn
        self.ntasks      = args.ntasks #number of tasks
        self.days      = int(args.days/1)
        self.hours    = int((args.days % 1) * 24)
        self.minutes    = int(((args.days % 1) * 24 % 1) * 60)
        self.memG       = args.memG
        #self.stderrdir  = args.stderr
        #self.stdoutdir  = args.stdout
        self.workd      = args.workd
        if args.depend:
            self.depend = '#SBATCH -d %s' % args.depend
        else:
            self.depend = ''

        if args.nodes:
            self.nodes = '#SBATCH -w %s' % args.nodes
        else:
            self.nodes = ''
			
        if args.gpu > 0:#--cpus-per-gpu is not compatible with -c
            self.gpu = '#SBATCH -p gpuq \n#SBATCH --gres=gpu:1'
        else:
            self.gpu = ''

    def next_index_auto(self):
        batch_index = 0
        while (True):
            if not os.path.exists('%s/j%d_%s.pbs' %
                                  (self.pbsdir, batch_index, self.nameroot)):
                return batch_index
            batch_index += 1

    def gen(self, batch_index, jobname=None):

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
            result.extend([lst[_] for _ in range(start, end)])

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
    pbs_filename=job.gen(int(args.index),jobname=args.name)
    if not args.silent:
        sys.stderr.write('==================\n\n')
        sys.stderr.write('pbs path: %s\n' % pbs_filename)
        sys.stderr.write('pbs stdout: %s\n' % job.stdout)
        sys.stderr.write('pbs stderr: %s\n' % job.stderr)
        sys.stderr.write('memG: %d\n' % job.memG)
        sys.stderr.write('ntasks: %d\n' % job.ntasks)
        sys.stderr.write('time limit: %d-%d:%d\n' %(job.days,job.hours,job.minutes))
        sys.stderr.write('ppn: %d\n' % job.ppn)
        sys.stderr.write('==================\n\n')
    if args.submit:
        subprocess.check_call(['sbatch',pbs_filename])

def add_default_settings(parser, d):

    """ default setting """

    # default_pbsdir = os.environ['PBSDIR'] if 'PBSDIR' in os.environ else os.getcwd()
    default_pbsdir = os.path.join(os.getcwd(),'pbs')
    parser.add_argument("-pbsdir", default=default_pbsdir,
                        help="pbs directory [%s]" % default_pbsdir)
    parser.add_argument("-ppn", default=d.ppn, type=int,
                        help="ppn [%d]" % d.ppn)
    parser.add_argument("-ntasks", default=d.ntasks, type=int,
                        help="ntasks [%d]" % d.ntasks)

    parser.add_argument("-days", default=d.days, type=float,
                        help="time limit (days), could be a int or float [%f]" % d.days)
    parser.add_argument("-workd", default=os.getcwd(),
                        help="workd [%s]" % d.workd)
    default_nameroot = os.environ['NAMEROOT'] if 'NAMEROOT' in os.environ else 'LabJob'
    parser.add_argument("-nameroot", default=default_nameroot,
                        help="job name root [%s]" % default_nameroot)
    parser.add_argument("-memG", default=d.memG, type=int,
                        help="memory in G [%d]" % d.memG)
    parser.add_argument('-depend', default='', help='dependency')
    parser.add_argument('-nodes', default='', help='nodes, for example: -nodes c-0-08, separated by comma if more than one nodes requested.')
    parser.add_argument('-gpu', default=0, help='gpu',type=int)
    parser.add_argument('-silent', action='store_true', help='suppress info print')

def pbsgen_main(parser):
    suffix_maping_dict={'r':'Rscript','py':'python','sh':'sh'}
    parser = argparse.ArgumentParser(description="generate pbs files")
    parser.add_argument('command', nargs='*', type=str, default='-', help='command to run')
    parser.add_argument("-index", type=int, default=-1, help="index of next pbs script [inferred from existing file names]")
    parser.add_argument('-name', default=None, help='job name (if specified, I will ignore index), if given an absolute path, use the basename as job name, use the dirname as pbsdir')
    add_default_settings(parser, Default())
    parser.add_argument('-submit', action='store_true', help='submit pbs')
    args = parser.parse_args()
    # print(args)
    if args.command != '-':
        if len(args.command)==1 and os.path.isfile(args.command[0]): #if the commond is a file, add prefix automatically.
            suffix=args.command[0].split('.')[-1].lower()
            if suffix in suffix_maping_dict:
                if args.name is None:
                    args.name=args.command[0]
                args.command=[suffix_maping_dict[suffix]]+args.command
        args.command=' '.join(args.command)
    if not os.path.exists(args.pbsdir):
        os.mkdir(args.pbsdir)
    main(args)
    return

class Default:
    def __init__(self):
        self.pbsdir = os.getcwd()
        # self.stdoutdir = '%s/stdout/' % self.pbsdir
        # self.stderrdir = '%s/stderr/' % self.pbsdir
        self.ppn = 1
        self.ntasks = 1
        self.days = 2.25
        self.memG = 10
        self.workd = os.getcwd()

if __name__ == "__main__":
    pbsgen_main(Default())

