#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import subprocess
import logging
import gzip
import time
import multiprocessing
if sys.version_info[0] < 3:
    import Queue
    import string
else:
    import queue as Queue

import numpy

from . Base import Base
from . import string_trimming


class RE(Base):
    """Class for creating and handling restriction enzyme fragment maps
    and counts.
    """
    arguments = ['bowtie', '-r', '--norc', '-v', '0', '-a', '-S',
                 '--sam-nohead', '--sam-nosq']
    data = None
    chroms = None

    def __init__(self, fname, index=None, focus=None, verbose=0, logger=None):
        super(RE, self).__init__(verbose=verbose, log_name='RE')
        if focus is not None and focus not in ['fragments', 'sites']:
            self.logging.error('The focus must be either fragments or sites')
            return None
        self.focus = focus
        if fname is not None:
            self.load_data(fname)
        else:
            self.map_RE(index)

    def load_data(self, fname):
        """Load RE fragment data and possibly scores from a bed or
        bedgraph file
        """
        self.logger.info("Loading RE data from file")
        if not os.path.exists(fname):
            self.logger.error("The RE fragment bed/bedgraph file %s could "
                              "not be found" % fname)
            return None
        data = []
        with open(fname) as f:
            for line in f:
                if line[0] == '#':
                    continue
                try:
                    line = line.decode('UTF-8').rstrip('\n').split('\t')
                    chrom, start, stop = line[0], int(line[1]), int(line[2])
                except (IndexError, ValueError):
                    continue
                try:
                    # First assume file is a bedgraph
                    score = float(line[3])
                except (IndexError, ValueError):
                    try:
                        # Next see if it is a bed file with valid score data
                        score = float(line[4])
                    except (IndexError, ValueError):
                        score = 0.
                data.append((chrom, (start, stop), score))
        data = numpy.array(data, dtype=numpy.dtype([
            ('chr', 'S20'), ('coords', numpy.int32, (2,)),
            ('score', numpy.float64)]))
        if data.shape[0] == 0:
            self.logger.warn("No bed data could be loaded from the "
                             "specified file")
            return None
        # Determine if data are RE sites or fragments
        size = numpy.mean(data['coords'][:, 1] - data['coords'][:, 0])
        if size <= 4:
            data_focus = 'sites'
            # if the coordinates are the center of the cut sites, shift to
            # start and end of recognition site
            if size == 0:
                data['coords'][:, 0] -= 2
                data['coords'][:, 1] += 2
        else:
            data_focus = 'fragments'
        data = data[numpy.argsort(data['coords'][:, 0])]
        # Order chromosomes by name/number
        chroms = numpy.unique(data['chr'])
        chrints = []
        for i in range(chroms.shape[0]):
            try:
                chrints.append((
                    str(int(chroms[i].lstrip('chr'))).rjust(2, '0'),
                    chroms[i]))
            except ValueError:
                chrints.append((chroms[i], chroms[i]))
        chrints.sort()
        chroms = []
        for i in range(len(chrints)):
            chroms.append(chrints[i][1])
        self.chroms = numpy.array(chroms)
        self.chr_indices = numpy.zeros(self.chroms.shape[0] + 1,
                                       dtype=numpy.int32)
        # Determine if requested focus matches current data focus
        if self.focus is None:
            self.focus = data_focus
            self.logger.info("Defaulting to %s-focused analysis" % data_focus)
        if data_focus == self.focus:
            N = data.shape[0]
        elif data_focus == 'fragments':
            N = data.shape[0] + self.chroms.shape[0]
        else:
            N = data.shape[0] - self.chroms.shape[0]
        self.data = numpy.zeros(N, dtype=numpy.dtype([
            ('chr', numpy.int32), ('coords', numpy.int32, (2,)),
            ('treatment', numpy.int32), ('control', numpy.int32),
            ('score', numpy.float64), ('alignable', numpy.bool)]))
        self.data['alignable'].fill(True)
        for i in range(self.chroms.shape[0]):
            where = numpy.where(data['chr'] == self.chroms[i])[0]
            start = self.chr_indices[i]
            if data_focus == self.focus:
                self.chr_indices[i + 1] = start + where.shape[0]
                stop = self.chr_indices[i + 1]
                self.data['coords'][start:stop, :] = data['coords'][where, :]
                # The only case where we keep scores is if the requested and
                # data focuses match
                self.data['score'][start:stop] = data['score'][where]
            elif data_focus == 'fragments':
                span = (data['coords'][where[1], 0]
                        - data['coords'][where[0], 1])
                self.chr_indices[i + 1] = start + where.shape[0] + 1
                stop = self.chr_indices[i + 1]
                self.data['coords'][start, 0] = (data['coords'][where[0], 0]
                                                 - span)
                self.data['coords'][start:(stop - 1), 1] = (
                    data['coords'][where, 0])
                self.data['coords'][(start + 1):stop, 0] = (
                    data['coords'][where, 1])
                self.data['coords'][stop - 1, 1] = (
                    data['coords'][where[-1], 1] + span)
            else:
                self.chr_indices[i + 1] = start + where.shape[0] - 1
                self.data['coords'][start:stop, 0] = (
                    data['coords'][where[:-1], 1])
                self.data['coords'][start:stop, 1] = (
                    data['coords'][where[1:], 0])
            self.data['chr'][start:stop] = i

    def map_RE(self, index):
        """Find RE sites by mapping recognition sequence to genome"""
        if index is None:
            self.logger.error("The bowtie genome index must be specified to "
                              "map restriction enzyme sites")
            return None
        self.logger.info("Mapping restriction enyzme recognition sites")
        # Start bowtie as a subprocess
        self.logger.debug("Running command '%s'" % ' '.join(self.arguments
                                                            + [index, '-']))
        mapping = subprocess.Popen(
            self.arguments + [index, '-'], stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Send the raw sequence of the DpnII recognition site
        mapping.stdin.write(b'GATC')
        mapping.stdin.close()
        bed = {}
        total = 0
        # Retrieve the alignments from bowtie
        with mapping.stdout as f:
            for line in f:
                line = line.decode('UTF-8').split('\t')
                chrom, start = line[2], int(line[3])
                stop = start + 4
                if chrom not in bed:
                    bed[chrom] = []
                bed[chrom].append((start, stop))
                total += 1
        # Log mapping results
        with mapping.stderr as f:
            for line in f:
                if line[0] == '#':
                    continue
                self.logger.debug(line.decode('UTF-8').rstrip('\n'))
        # Sort chromosome list by name/number
        chroms = numpy.array(list(bed))
        chrints = []
        for i in range(chroms.shape[0]):
            try:
                chrints.append((
                    str(int(chroms[i].lstrip('chr'))).rjust(2, '0'),
                    chroms[i]))
            except ValueError:
                chrints.append((chroms[i], chroms[i]))
        chrints.sort()
        chroms = []
        for i in range(len(chrints)):
            chroms.append(chrints[i][1])
        self.chroms = numpy.array(chroms)
        self.chr_indices = numpy.zeros(self.chroms.shape[0] + 1,
                                       dtype=numpy.int32)
        if self.focus is None:
            self.logger.info("Defaulting to a fragment-focused analysis")
            self.focus = 'fragments'
        if self.focus == 'fragments':
            N = total - self.chroms.shape[0]
        else:
            N = total
        # Arrange data into single array with indexed chromosomes
        self.data = numpy.zeros(N, dtype=numpy.dtype([
            ('chr', numpy.int32), ('coords', numpy.int32, (2,)),
            ('treatment', numpy.int32), ('control', numpy.int32),
            ('score', numpy.float64), ('alignable', numpy.bool)]))
        self.data['alignable'].fill(True)
        for i in range(self.chroms.shape[0]):
            chrom = self.chroms[i]
            bed[chrom] = numpy.array(bed[chrom])
            bed[chrom] = bed[chrom][numpy.argsort(bed[chrom][:, 0]), :]
            start = self.chr_indices[i]
            if self.focus == 'fragments':
                self.chr_indices[i + 1] = start + bed[chrom].shape[0] - 1
                stop = self.chr_indices[i + 1]
                self.data['coords'][start:stop, 0] = bed[chrom][:-1, 1]
                self.data['coords'][start:stop, 1] = bed[chrom][1:, 0]
            else:
                self.chr_indices[i + 1] = start + bed[chrom].shape[0]
                stop = self.chr_indices[i + 1]
                self.data['coords'][start:stop, :] = bed[chrom]
            self.data['chr'][start:stop] = i


class Fastq(Base):
    """Class for creating a BAM alignment file from a FASTQ file,
    processing reads as requested.
    """
    arguments = ['bowtie', '-q', '-n', '2', '-y', '--best', '--strata',
                 '-m', '1', '-S']

    def __init__(self, fastq, index, threads=1, outdir='./', seed=None,
                 verbose=0, multimapping=False, trim_size=10):
        super(Fastq, self).__init__(verbose=verbose,
                                    log_name='Fastq-%s' % fastq)
        if not os.path.exists(fastq):
            self.logger.error("The fastq file %s could not be found" % fastq)
            return None
        if index is None:
            self.logger.error("The bowtie genome index must be specified")
            return None
        name = fastq.split('/')[-1]
        self.logger.info("Mapping %s" % name)
        # Set output names
        self.outdir = outdir
        bam_fname = "%s/%s.bam" % (
            self.outdir, name.split('.fastq')[0].replace('.preprocessed', ""))
        stats_fname = "%s/%s.stats" % (
            self.outdir, name.split('.fastq')[0].replace('.preprocessed', ""))
        # Check whether a seed was provided
        if seed is None:
            seed = numpy.random.randint(10000)
        args = [['--seed', str(seed), '--un',
                 fastq.replace('.fastq', '.un.fastq'), '--max',
                 fastq.replace('.fastq', '.mult.fastq')],
                ['--seed', str(seed + 1), '--un',
                 fastq.replace('.fastq', '.un1.fastq'), '-5', str(trim_size),
                 '--max', fastq.replace('.fastq', '.mult1.fastq')],
                ['--seed', str(seed + 2), '--un',
                 fastq.replace('.fastq', '.un2.fastq'), '-3', str(trim_size)]]
        fastq_fnames = [fastq,
                        fastq.replace('.fastq', '.un.fastq'),
                        fastq.replace('.fastq', '.un1.fastq')]
        bam_fnames = [bam_fname,
                      bam_fname.replace('.bam', '.1.bam'),
                      bam_fname.replace('.bam', '.2.bam')]
        for i in range(1 + 2 * bool(multimapping)):
            bamfile = open(bam_fnames[i], 'wb')
            statsfile = open(bam_fnames[i].replace('.bam', '.stats'), 'wb')
            # Start bowtie subprocess
            self.logger.debug("Running command '%s | samtools view -b - > %s'" % (' '.join(
                self.arguments + args[i] + ['--threads', '%i' % threads,
                index, fastq_fnames[i]]), bam_fnames[i]))
            mapping_task = subprocess.Popen(
                self.arguments + args[i] + ['--threads', '%i' % threads,
                                            index, fastq_fnames[i]],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            # Start samtools subprocess
            samtools_task = subprocess.Popen(['samtools', 'view', '-b',
                                              '-'],
                                             bufsize=1,
                                             stdin=mapping_task.stdout,
                                             stdout=bamfile)
            # Start grep subprocess to remove stats output from bowtie that
            # we don't care about
            grep_task = subprocess.Popen(['grep', '-v', 'Warning'],
                                         bufsize=1,
                                         stdin=mapping_task.stderr,
                                         stdout=statsfile)
            samtools_task.wait()
            grep_task.wait()
            bamfile.close()
            statsfile.close()


class Preprocessor(Base):

    def __init__(self, fastq, outdir='./', threads=1, quality_trim=True,
                 quality_cutoff=30, minsize=20, verbose=0,
                 split_by=[b'CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA'],
                 adapters=[b'GATCCTCGGCCGCGACCGGTCGCGGCCGAGGATC',
                           b'CTAATACGACTCACTATAGGGCAGCGTGGTCGCGGCCGAGGA']):
        super(Preprocessor, self).__init__(verbose=verbose,
                                           log_name='Preprocessor-%s' % fastq)
        if not os.path.exists(fastq):
            self.logger.error("The fastq file %s could not be found" % fastq)
        self.logger.info("Preprocessing %s" % fastq.split('/')[-1])
        # Set output names
        self.outdir = outdir
        preprocessed_fname = "%s/%s.preprocessed.fastq" % (
            self.outdir, fastq.split('/')[-1].split('.fastq')[0])
        # Set python version-specific functions
        if sys.version_info[0] < 3:
            self.maketrans = string.maketrans
            self.find = string.find
        else:
            self.maketrans = bytes.maketrans
            self.find = bytes.find
        # Preprocess all adapter sequences and get reverse complements
        # if not already in the list
        adapter_tables = []
        for adapter in adapters:
            adapter_tables.append([
                [adapter[::-1], numpy.zeros(len(adapter), dtype=numpy.int32)],
                [adapter, numpy.zeros(len(adapter), dtype=numpy.int32)]])
            string_trimming.make_kmp_table(adapter_tables[-1][0][1],
                                           adapter_tables[-1][0][0])
            string_trimming.make_kmp_table(adapter_tables[-1][1][1],
                                           adapter_tables[-1][1][0])
            revcomp = self.revcomp(adapter)
            if revcomp != adapter:
                adapter_tables.append([
                    [revcomp[::-1], numpy.zeros(len(adapter),
                                                dtype=numpy.int32)],
                    [revcomp, numpy.zeros(len(adapter), dtype=numpy.int32)]])
                string_trimming.make_kmp_table(adapter_tables[-1][0][1],
                                               adapter_tables[-1][0][0])
                string_trimming.make_kmp_table(adapter_tables[-1][1][1],
                                               adapter_tables[-1][1][0])
        # Add reverse complements to split_by list if not already present
        N = len(split_by)
        for i in range(N):
            revcomp = self.revcomp(split_by[i])
            if revcomp not in split_by:
                split_by.append(revcomp)
        if sys.version_info[0] >= 3:
            for i in range(len(split_by)):
                split_by[i] = split_by[i].decode('UTF-8')
        self.quality_trim = quality_trim
        self.quality_cutoff = quality_cutoff
        self.split_by = split_by
        if len(split_by) > 0:
            self.split = True
        else:
            self.split = False
        self.adapter_tables = adapter_tables
        if len(adapter_tables) > 0:
            self.strip = True
        else:
            self.strip = False
        # Set python version-specific functions
        if sys.version_info[0] < 3:
            self.maketrans = string.maketrans
            self.find = string.find
        else:
            self.maketrans = bytes.maketrans
            self.find = bytes.find
        # Create output queue
        out_queue = multiprocessing.JoinableQueue(maxsize=(1000 * threads))
        count_queue = multiprocessing.JoinableQueue()
        # Define Process function

        def preprocess(fname, fs_func, pattern, out_queue, start, stop,
                       minsize, process_read, total_reads):
            infile = fs_func(fname)
            infile.seek(start)
            new_pos = start
            pos = start
            data = []
            count = 0
            while pos < stop:
                # Get next chunk of data
                pos = new_pos
                new_pos = min(pos + 10000, stop)
                chunk = new_pos - pos
                new_data = infile.read(chunk).split('\n')
                # If this is the first chunk, get rid of partial reads
                if pos == start:
                    while (new_data
                           and new_data[0][:min(len(new_data[0]),
                                                len(pattern))] != pattern):
                        new_data.pop(0)
                if len(data) > 0:
                    data[-1] += new_data.pop(0)
                    if data[0] == '':
                        data.pop(0)
                data += new_data
                # Break lines into reads (4 lines per)
                dpos = 0
                for i in range(0, len(data) - 3, 4):
                    read = [data[i], data[i + 1], data[i + 2], data[i + 3]]
                    # make sure last line isn't partial
                    if len(read[3]) == len(read[1]):
                        count += 1
                        reads = process_read(read)
                        # Put processed read (which may be broken into
                        # multiple reads) into outgoing queue
                        for read in reads:
                            if len(read[1]) >= minsize:
                                out_queue.put(read, True)
                        dpos += 4
                data = data[dpos:]
            # Check if last read wasn't completely read
            if len(data) > 0:
                new_data = infile.read(10000).split('\n')
                data[-1] += new_data.pop(0)
                data += new_data
                if len(data) >= 4:
                    read = [data[0], data[1], data[2], data[3]]
                    count += 1
                    # make sure last line isn't partial
                    if len(read[3]) == len(read[1]):
                        reads = process_read(read)
                        # Put processed read (which may be broken into
                        # multiple reads) into outgoing queue
                        for read in reads:
                            if len(read[1]) >= minsize:
                                out_queue.put(read, True)
            infile.close()
            out_queue.put(None, True)
            total_reads.put(count)

        # Determine if fastq file is gzipped
        try:
            infile = gzip.GzipFile(fastq, 'rb')
            line = infile.readline()
            pattern = line[:3]
            fs_func = gzip.GzipFile
            # Determine fastq.gz file size
            while True:
                infile.seek(100000000, 1)
                if infile.read(1) == '':
                    break
            fsize = infile.tell()
            infile.close()
        except IOError:
            infile = open(fastq)
            line = infile.readline()
            pattern = line[:3]
            infile.close()
            fs_func = open
            # Determine fastq file size
            s = subprocess.Popen(['ls', '-l', fastq], stdout=subprocess.PIPE)
            fsize = int(s.stdout.readline().split()[4])
        self.logger.debug("File size %i" % fsize)
        outfile = open(preprocessed_fname, 'while')
        # Determine process byte ranges
        proc_range = numpy.round(numpy.linspace(
            0, fsize, max(2, threads))).astype(numpy.int64)
        # Start Processes
        processes = []
        for i in range(max(1, threads - 1)):
            processes.append(multiprocessing.Process(
                target=preprocess, args=(fastq, fs_func, pattern, out_queue,
                                         proc_range[i], proc_range[i + 1],
                                         minsize, self.preprocess_read,
                                         count_queue)))
            processes[-1].daemon = True
            processes[-1].start()
        # Start writing finished reads to file
        count = 0
        finished = 0
        while True:
            read = out_queue.get(True)
            count += 1
            if read is None:
                finished += 1
                if finished == max(1, threads - 1):
                    break
            else:
                outfile.write(read[0] + '\n')
                outfile.write(read[1] + '\n')
                outfile.write(read[2] + '\n')
                outfile.write(read[3] + '\n')
        outfile.close()
        count = 0
        finished = 0
        for i in range(max(1, threads - 1)):
            c = count_queue.get(True)
            count += c
        self.logger.debug("Processed %i reads" % count)

    def preprocess_read(self, read):
        reads = [read]
        if self.quality_trim:
            self.trim_quality(reads)
        if self.strip:
            self.strip_adapters(reads)
        if self.split:
            self.split_string(reads)
        return reads

    def revcomp(self, seq):
        """Find the reverse complement of a sequence."""
        tab = self.maketrans(b'ACNGT', b'TGNCA')
        return seq.translate(tab)[::-1]

    def trim_quality(self, reads):
        """Trim off ends of read using a sliding 3-base window for
        finding average quality.
        """
        cut = self.quality_cutoff * 3
        start = 0
        qscores = reads[0][3]
        qual = ord(qscores[0]) + ord(qscores[1]) + ord(qscores[2]) - 99
        while qual < cut:
            start += 1
            try:
                qual += ord(qscores[start + 2]) - ord(qscores[start - 1])
            except IndexError:
                break
        stop = len(qscores)
        qual = ord(qscores[-1]) + ord(qscores[-2]) + ord(qscores[-3]) - 99
        while qual < cut:
            stop -= 1
            try:
                qual += ord(qscores[stop - 3]) - ord(qscores[stop])
            except IndexError:
                break
        reads[0][1] = reads[0][1][start:stop]
        reads[0][3] = reads[0][3][start:stop]

    def strip_adapters(self, reads):
        """Using a Knuth-Morris-Pratt type algorith, identify and strip
        overlapping adapter sequences from the ends of read.
        """
        for a_set in self.adapter_tables:
            M = len(reads[0][1])
            N = min(M, len(a_set[0][0]))
            start = string_trimming.overlap(
                a_set[0][0], reads[0][1][:N][::-1], a_set[0][1])
            stop = M - string_trimming.overlap(
                a_set[1][0], reads[0][1][-N:], a_set[1][1])
            if stop - start < M:
                reads[0][1] = reads[0][1][start:stop]
                reads[0][3] = reads[0][3][start:stop]

    def split_string(self, reads):
        """Split reads based on target sequences, removing target sequence."""
        r = 0
        while r < len(reads):
            for pattern in self.split_by:
                index = reads[r][1].find(pattern)
                if index >= 0:
                    pos = index + len(pattern)
                    reads.append([])
                    reads[-1].append(reads[r][0])
                    reads[-1].append(reads[r][1][pos:])
                    reads[-1].append(reads[r][2])
                    reads[-1].append(reads[r][3][pos:])
                    reads[r][1] = reads[r][1][:index]
                    reads[r][3] = reads[r][3][:index]
            r += 1
        if len(reads) > 1:
            for r in range(len(reads))[::-1]:
                if len(reads[r][1]) < 25:
                    reads.pop(r)
        if len(reads) > 1:
            for r in range(len(reads)):
                reads[r][0] += ('.%i' % r)
