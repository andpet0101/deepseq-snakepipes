#!/usr/bin/env python
'''
The MIT License (MIT)

Copyright (c) <2018> <Mathias Lesche>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

Check out https://www.python.org/dev/peps/pep-0008/#naming-conventions for naming
conventions

contact: mathias.lesche(at)tu-dresden.de
'''

''' python modules '''
import logging

from argparse import ArgumentParser as ArgumentParser
from argparse import RawDescriptionHelpFormatter


class MainLogger(object):
    def __init__(self, logtitle, streamh = True, fileh = False, logfilename = ''):
        self.__logfilename = logfilename
        self.__ch = ''
        self.__streamh = streamh
        self.__fh = ''
        self.__fileh = fileh
        self.set_logger(logtitle, streamh, fileh)

    def set_logger(self, logtitle, streamh, fileh):
        self.logger = logging.getLogger(logtitle)
        self.logger.setLevel(logging.DEBUG)
        self.formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s', datefmt= '%m/%d/%Y %I:%M:%S %p')

        if fileh:
            fh = logging.FileHandler(self.__logfilename)
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(self.formatter)
            self.logger.addHandler(fh)
            self.__fh = fh

        if streamh:
            ch = logging.StreamHandler()
            ch.setLevel(logging.INFO)
            ch.setFormatter(self.formatter)
            self.logger.addHandler(ch)
            self.__ch = ch

    def add_filelogger(self, filename):
        self.__logfilename = filename
        fh = logging.FileHandler(self.__logfilename)
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(self.formatter)
        self.logger.addHandler(fh)
        self.__fh = fh
        self.__fileh = True  
        
    def add_streamlogger(self):
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(self.formatter)
        self.logger.addHandler(ch)
        self.__ch = ch
        self.__streamh = True

    def get_logger(self):
        return self.logger
    
    def close(self):
        if self.__fileh: self.logger.removeHandler(self.__fh)
        if self.__streamh: self.logger.removeHandler(self.__ch)

    def get_fileh(self):
        return self.__fileh

    fileh = property(get_fileh)

class Parser(object):
    def __init__(self):
        self.__parser = ArgumentParser(description="""
        Wrapper function for the snakemake pipeline
        """, formatter_class=RawDescriptionHelpFormatter, 
        add_help = False, prog = 'snakemake_wrapper.py')
        
        self.initialiseParser()
        self.__drmaa = '--drmaa'
        self.__submission = ''

        self.__logger = logging.getLogger('snakemake.wrapper')
        self.parse()

    def initialiseParser(self):
        self.__parser.add_argument('-h', '--help', dest = 'h', help = 'Display help message', action= 'store_true')
        self.__wrapper = self.__parser.add_argument_group("snakemake specific arguments")
        self.__wrapper.add_argument('-d', '--drmaa', dest = 'use_drmaa', action = 'store_true', help = 'use DRMAA for submitting/controlling cluster jobs')
        #self.__parser.add_argument('-s', '--submission', dest='submission', metavar='STRING', choices = ('sge', 'slurm'), required = True, help = 'job submission system')
    #     self.__parser.add_argument('-f', '--files', type=str, metavar='FILE', dest='file', nargs= '+', help="list of fq-files(' ' separated)")
        #self.__parser.add_argument('-o', '--output', type=str, metavar='STRING', dest='output', required=True, help='id of the run which will be used as output folder')
        #self.__parser.add_argument('-u', '--undet', dest='undet', action='store_true', help='transfer undetermined')
        #self.__parser.add_argument("-e", "--email", dest = 'email', action= 'store_true', help= 'activates email notification')
        #self.__parser.add_argument('-q', '--qstat', type=str, metavar='STRING', dest='qstatname', default='tsp_fq', help='qstat name of the qsub job (default: tsp_fq)')

    def parse(self, inputstring = None):
        if inputstring == None:
            self.__options = self.__parser.parse_args()
        else:
            self.__options = self.__parser.parse_args(inputstring)

    def show_log(self, level, message):
        if level == 'debug':
            self.__logger.debug(message)
        elif level == 'info':
            self.__logger.info(message)
        elif level == 'warning':
            self.__logger.warning(message)
        elif level == 'error':
            self.__logger.error(message)
        elif level == 'critical':
            self.__logger.critical(message)

    def main(self):
        print(self.__options)
        print(self.__wrapper)

        if self.__options.h: self.__parser.print_help()

        
        
        #if not self.__options.use_drmaa: self.__drmaa = ''
        #print(self.__drmaa)
        #self.__submission = self.__options.submission
        #print(self.__submission)
        #self.__qstatname = self.__options.qstatname

        #for dirname in self.check_directory(self.__options.directory):
        #if dirname  != '': self.__files.extend(add_FilestoList_recursive(dirname, [], ('fastq.gz', 'fastq', 'fq.gz', 'fq', 'bam', 'sam')))
    ##     remove empty lines
        #self.__files = sorted(set(self.__files))
        #self.__files = [i for i in self.__files if i != '']

        #if len(self.__files) == 0:
        #self.show_log('error', "valid directory (-d) wasn't provided")
        #exit(2)
        #else:
        #self.show_log('info', '{0} files for transfer'.format(len(self.__files)))

        #self.__output = self.__options.output
        #self.show_log('info', 'storage folder on {0}: {1}'.format(Information.FQSTORAGE, self.__output))
        #self.__undet = self.__options.undet
        #self.show_log('info', 'transfer undetermined fq-files: {0}'.format(self.__undet))
        #sleep(3)


if __name__ == '__main__':
    mainlog = MainLogger('snakemake')
    parser = Parser()
    parser.main()

    mainlog.close()
    logging.shutdown()