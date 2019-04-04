#!/usr/bin/env python

from __future__ import print_function

import os
import glob
from argparse import ArgumentParser

from . Base import Base


class Parser(Base):

    def __init__(self, verbose=2):
        super(Parser, self).__init__(verbose=verbose, log_name='Parser')

    def parse_args(self):
        self.Generate_Parser()
        args = self.parser.parse_args()

        # Determine if RE input options are valid
        if (
                not os.path.exists("%s/%s.npz" % (args.NAME, args.NAME))
                and args.RE is None and args.INDEX is None):
            self.logger.error("A saved Segmenter file, RE bed file or a "
                              "bowtie index must be passed")
            return None
        if args.RE is not None and not os.path.exists(args.RE):
            self.logger.error("Could not locate %s" % args.RE)

        # Determine if only one input type is passed
        if args.FASTQ is not None and args.BAM is not None:
            self.logger.warn("Both fastq and bam/sam file arguments were "
                             "passed. Ignoring the fastq arguments")
            args.FASTQ = None
        if args.BAM is not None and args.COUNTS is not None:
            self.logger.warn("Both bam/sam and count bed file arguments were "
                             "passed. Ignoring the bam/sam arguments")
            args.BAM = None

        # If an bowtie index is provided, make sure the path is correct
        if args.INDEX is not None:
            if not os.path.exists('%s.1.ebwt' % args.INDEX):
                self.logger.error("Could not locate the bowtie index %s" %
                                  args.INDEX)
                return None

        # If passing fastq files, determine if necessary arguments are present
        if args.FASTQ is not None:
            for i in range(2):
                args.FASTQ[i] = args.FASTQ[i].split(',')
                # Check if it could be a regular expression
                if len(args.FASTQ[i]) == 1:
                    args.FASTQ[i] = glob.glob(args.FASTQ[i][0])
                for fname in args.FASTQ[i]:
                    if not os.path.exists(fname):
                        self.logger.error("Could not locate %s" % fname)
                        return None
            if args.INDEX is None:
                self.logger.error("A bowtie index must be provided to map "
                                  "fastq reads")
                return None

        # If passing bam files, determine if necessary arguments are present
        if args.BAM is not None:
            for i in range(2):
                args.BAM[i] = args.BAM[i].split(',')
                # Check if it could be a regular expression
                if len(args.BAM[i]) == 1:
                    args.BAM[i] = glob.glob(args.BAM[i][0])
                for fname in args.BAM[i]:
                    if not os.path.exists(fname):
                        self.logger.error("Could not locate %s" % fname)
                        return None

        # If passing count files, determine if necessary arguments are present
        if args.COUNTS is not None:
            for i in range(2):
                args.COUNTS[i] = args.COUNTS[i].split(',')
                # Check if it could be a regular expression
                if len(args.COUNTS[i]) == 1:
                    args.COUNTS[i] = glob.glob(args.COUNTS[i][0])
                for fname in args.COUNTS[i]:
                    if not os.path.exists(fname):
                        self.logger.error("Could not locate %s" % fname)
                        return None

        # Format remaining arguments correctly
        args.ADAPTERS = args.ADAPTERS.encode('UTF-8').split(b',')
        args.SPLIT = args.SPLIT.encode('UTF-8').split(b',')
        if os.path.dirname(args.NAME) == '':
          args.OUTDIR = '.'
        else:
          args.OUTDIR = os.path.dirname(args.NAME)
        return args

    def Generate_Parser(self):
        desc = ("A start to finish analysis program for taking sequencing "
                "data, mapped data, or scores produced by DamID and calling "
                "differentially marked regions.")
        self.parser = ArgumentParser(prog='pLADetector', description=desc,
                                     add_help=True)
        self.parser.add_argument('-v', '--verbose', dest='VERBOSE',
                                 type=int, action='store', default=1,
                                 help='Level of log reporting.')
        self.add_RE_options()
        self.add_fastq_options()
        self.add_count_options()
        self.add_segmentation_options()
        self.add_output_options()

    def add_RE_options(self):
        group = self.parser.add_argument_group("Restriction Site Options")
        group.add_argument("--re", dest="RE", type=str, action='store',
                           default=None, help="Restriction enzyme sites or "
                           "fragments file in BED or BEDGRAPH format.")
        group.add_argument('--unalignable', dest='UNALIGN', type=str,
                           action='store', default=None, help="BED file "
                           "containing blacklisted regions to exclude from "
                           "analysis.")
        group.add_argument('--focus', dest='FOCUS', type=str, action='store',
                           default=None, choices=['fragments', 'sites', None],
                           help="Focus of analysis, RE sites, RE fragments, "
                           "or default to the format of the RE file.")
        group.add_argument('--index', dest='INDEX', type=str, action='store',
                           default=None, help='Bowtie index.')

    def add_fastq_options(self):
        group = self.parser.add_argument_group('Fastq Options')
        group.add_argument('--fastq', dest='FASTQ', nargs=2, type=str,
                           action='store', help="A pair of comma-separated "
                           "lists or regular expressions of treatment and "
                           "control fastq files, respectively.")
        group.add_argument('--threads', dest='THREADS', action='store',
                           default=1, type=int, help="Number of threads to "
                           "use for bowtie mapping.")
        group.add_argument('--split', dest='SPLIT', type=str, action='store',
                           help="A comma-separated list of sequences to split "
                           "reads by.",
                           default='CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA'
                           )
        group.add_argument('--adapters', dest='ADAPTERS', type=str,
                           action='store',
                           default='GATCCTCGGCCGCGACCGGTCGCGGCCGAGGATC',
                           help="A comma-separated list of adapter sequences "
                           "to strip from read ends.")
        group.add_argument('--quality-trim', dest='QTRIM', action='store_true',
                           help="Trim ends of reads based on a windowed three "
                           "base average quality.")
        group.add_argument('--quality-cutoff', dest='QCUTOFF', type=float,
                           action='store', default=30., help="The minimum "
                           "quality average used by quality trimming.")
        group.add_argument('--min-read', dest='MINREAD', type=int,
                           action='store', default=20, help="The minimum "
                           "read size to keep after preprocessing.")
        group.add_argument('--multimapping', dest='MULTIMAPPING',
                           action='store_true', help='Add two additional '
                           "rounds of bowtie mapping with 3' and 5' "
                           'trimming to recover more reads.')
        group.add_argument('--trim-size', dest='TRIMSIZE', type=int,
                           action='store', default=10, help="The number "
                           "of bases to trim from 3' and 5' read ends in "
                           "multimapping approach.")
        group.add_argument('--keep-preprocessed', dest='KEEPPRE',
                           action='store_true', help='Keep the preprocessed '
                           'fastq read file.')
        group.add_argument('--seed', dest='SEED', type=int,
                           action='store', default=None, help="A seed for the "
                           "random number generator.")

    def add_count_options(self):
        group = self.parser.add_argument_group('Count Options')
        group.add_argument('--bam', dest='BAM', nargs=2, type=str,
                           action='store', help="A pair of comma-separated "
                           "lists or regular expressions of treatment and "
                           "control bam/sam files, respectively.")
        group.add_argument('--counts', dest='COUNTS', nargs=2, type=str,
                           action='store', help="A pair of comma-separated "
                           "lists or regular expressions of treatment and "
                           "control bed/bedgraph count files, respectively.")
        group.add_argument('--max-re-dist', dest='MAXDIST', type=int,
                           action='store', default=-1, help="The maximum "
                           "distance from a restriction site a read can map "
                           "and still be considered valid. A negative value "
                           "indicates no maximum.")
        group.add_argument('--count-overlaps', dest='OVERLAP',
                           action='store_true', help="Count reads overlapping "
                           "restriction sites as counting towards both RE "
                           "fragments.")

    def add_segmentation_options(self):
        group = self.parser.add_argument_group('Segmentation Options')
        group.add_argument('--binsize', dest='BINSIZE', action='store',
                           type=int, default=0, help="Binsize to bin data "
                           "into prior to segmentation.")
        group.add_argument('--mindip', dest='MINDIP', action='store',
                           type=int, default=2000, help="Minimum size to "
                           "qualify as a DIP.")
        group.add_argument('--maxdip', dest='MAXDIP', action='store',
                           type=int, default=10000, help="Maximum size to"
                           "qualify as a DIP/gap to stitch together adjacent "
                           "LADs.")

    def add_output_options(self):
        group = self.parser.add_argument_group('Output Options')
        group.add_argument('--name', dest='NAME', type=str, action='store',
                           default='out', help="Prefix for results files and "
                           "name of output directory.")
        group.add_argument('--bedgraph', dest='BEDGRAPH', action='store_true',
                           help="Save RE treament and control counts to "
                           "bedgraph files.")
        group.add_argument('--plot', dest='PLOT', action='store_true',
                           help='Plot data to PDF file.')
        group.add_argument('--plot-binsize', dest='PBINSIZE', action='store',
                           type=int, default=10000, help="Bin data to this "
                           "resolution before plotting.")
        group.add_argument('--lads', dest='LADS', action='store_true',
                           help='Add LADs to plot.')
        group.add_argument('--dips', dest='DIPS', action='store_true',
                           help='Add DIPs to plot.')
        group.add_argument('--means', dest='MEANS', action='store_true',
                           help='Add segmentation means to plot.')
        group.add_argument('--segments', dest='SEGMENTS', action='store_true',
                           help='Add segmentation boundaries to plot.')
