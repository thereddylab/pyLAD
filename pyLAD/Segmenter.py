#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import logging
import time
import multiprocessing
if sys.version_info[0] < 3:
    import Queue
else:
    import queue as Queue

import numpy
from scipy.stats import norm

from . Mapping import RE
from . Base import Base
from . import dnacopy


class Segmenter(Base):
    """Class for scoring and segmenting mapped DamID data"""
    LADs = None
    DIPs = None
    sbdry = None
    unalignable = None
    data = None
    binned = None
    chroms = None
    partitions = None

    def __init__(self, verbose=0):
        super(Segmenter, self).__init__(verbose=verbose, log_name='Segmenter')

    def save(self, fname):
        """Save Segmenter object to a numpy npz file"""
        self.logger.info("Saving segmentation file")
        data = {}
        for key in self.__dict__:
            if key == 'verbose':
                continue
            # ints and floats need to be converted to numpy arrays to be saved
            if isinstance(self[key], int) or isinstance(self[key], float):
                data[key] = numpy.array([self[key]])
            # only ints, floats, and arrays should be saved
            elif isinstance(self[key], numpy.ndarray):
                data[key] = self[key]
        numpy.savez(fname, **data)

    def load(self, fname):
        """Load a Segmenter object from a numpy npz file"""
        self.logger.info("Loading segmentation file")
        if not os.path.exists(fname):
            self.logger.error("The Segmentation file %s could not be found"
                              % fname)
            sys.exit(0)
        try:
            data = numpy.load(fname)
        except IOError:
            self.logger.error("The specified file does not appear to be a "
                              "numpy npz file")
            sys.exit(0)
        for key in data:
            # The only length-1 arrays should be from ints and floats
            # converted during saving
            if len(data[key].shape) == 1 and data[key].shape[0] == 1:
                self[key] = data[key][0]
            else:
                self[key] = data[key]
        # Determine whether the data are binned, site-focused, or
        # fragment-focused
        if self.binned is not None:
            self.focus = 'binned'
        else:
            sizes = self.data['coords'][:, 1] - self.data['coords'][:, 0]
            if numpy.mean(sizes) <= 4:
                self.focus = 'sites'
            else:
                self.focus = 'fragments'
        # If RE data has been loaded, make our lookup table to chrom ints
        if 'chroms' in self.__dict__:
            self.chr2int = {}
            for i in range(len(self.chroms)):
                self.chr2int[self.chroms[i]] = i

    def load_RE(self, fname=None, index=None, focus='fragments'):
        """Load RE fragment data from a bed or bedgraph file, or map DpnII
        sites against the reference genome.
        """
        data = RE(fname=fname, index=index, focus=focus, verbose=self.verbose)
        if data.data is None:
            # Getting the data did not work
            return None
        self.data = data.data
        self.focus = data.focus
        self.chr_indices = numpy.r_[
            0, numpy.cumsum(numpy.bincount(self.data['chr']))]
        # Determine if blacklisted data (and with it chrom info) has already
        # been loaded
        if self.unalignable is not None:
            # Renumber unaligned chrom data
            uchroms = self.chroms
            self.chroms = data.chroms
            uindices = numpy.r_[
                0, numpy.cumsum(numpy.bincount(self.unalignable['chr']))]
            self.unalignable['chr'] = -1
            for i in range(self.chroms.shape[0]):
                where = numpy.where(uchroms == self.chroms[i])[0]
                if where.shape[0] == 0:
                    continue
                self.unalignable['chr'][uindices[where[0]]:
                                        uindices[where[0] + 1]] = i
            # Remove unalignable data from chroms not in RE data
            self.unalignable = (
                self.unalignable[numpy.where(self.unalignable['chr'] >= 0)[0]])
            # Reorder by new chromosome numbering
            self.unalignable = (
                self.unalignable[numpy.argsort(self.unalignable['chr'])])
            self.remove_unalignable()
        else:
            self.chroms = data.chroms
        self.chr2int = {}
        for i in range(self.chroms.shape[0]):
            self.chr2int[self.chroms[i]] = i

    def load_unalignable(self, fname):
        """Load blacklisted genomic regions from a bed file"""
        self.logger.info("Loading unalignable regions")
        if not os.path.exists(fname):
            self.logger.error("The unalignable sequence bed file %s could "
                              "not be found" % fname)
            return None
        # If no chrom data is present, this will define the current chrom list
        bed = []
        chroms = []
        chr2int = {}
        with open(fname) as f:
            for line in f:
                if line[0] == '#':
                    continue
                temp = line.rstrip('\n').split('\t')
                # If no chrom data is already loaded, use the bed file's
                # chrom data
                if self.data is None:
                    if temp[0] not in chroms:
                        chr2int[temp[0]] == len(chroms)
                        chroms.append(temp[0])
                    chrint = chr2int[temp[0]]
                elif temp[0] not in self.chr2int:
                    continue
                else:
                    chrint = self.chr2int[temp[0]]
                bed.append((chrint, (int(temp[1]), int(temp[2]))))
        bed.sort()
        bed = numpy.array(bed, dtype=numpy.dtype([
            ('chr', numpy.int32), ('coords', numpy.int32, (2,))]))
        # If getting chrom data from bed, find correct chrom ordering
        if self.data is None:
            indices = numpy.r_[0, numpy.cumsum(numpy.bincount(bed['chr']))]
            chrints = []
            for chrom in chroms:
                try:
                    chrints.append((
                        str(int(chroms[i].lstrip('chr'))).rjust(2, '0'),
                        chroms[i]))
                except ValueError:
                    chrints.append((chrom, chrom))
            chrints.sort()
            chroms = []
            for i in range(len(chrints)):
                chroms.append(chrints[i][1])
            self.chroms = numpy.array(chroms)
            for i in range(len(self.chroms)):
                chrint = chr2int[self.chroms[i]]
                bed['chr'][indices[chrint]:indices[chrint + 1]] = i
            bed = bed[numpy.argsort(bed['chr'])]
        self.unalignable = bed
        # If RE data is already loaded, remove blacklisted regions
        if self.data is not None:
            self.remove_unalignable()

    def load_reads(self, fname, condition='treatment', count_overlaps=True,
                   maxdist=-1):
        """Load reads from a SAM or BAM file for the specified condition."""
        if condition not in ['treatment', 'control']:
            self.logger.warn("Condition must be either 'treatment' or "
                             "'control'")
            return None
        if self.data is None:
            self.logger.warn("Restriction sites must be loaded prior to "
                             "loading reads")
            return None
        try:
            import pysam
        except ImportError:
            self.logger.warn("The package pysam must be installed to read "
                             "sam or bam files")
            return None
        self.logger.info('Counting reads from %s' % fname)
        infile = pysam.AlignmentFile(fname)
        seqs = {}
        # Get conversion from idx to chromosome name
        idx2chrom = {}
        for i in range(len(infile.header['SQ'])):
            idx2chrom[i] = infile.header['SQ'][i]['SN']
        for idx in idx2chrom:
            seqs[idx] = {}
        # Load mapped reads, chrom and coords
        for seq in infile:
            if seq.flag != 16 and seq.flag != 0:
                continue
            if seq.reference_id not in seqs:
                continue
            seqs[seq.reference_id][(seq.reference_start, seq.reference_start
                                    + len(seq.query_sequence))] = None
        # Parse reads by RE location data
        for idx in seqs:
            if idx2chrom[idx] not in self.chr2int:
                continue
            seqs[idx] = numpy.array(list(seqs[idx]))
            if seqs[idx].shape[0] == 0:
                continue
            chrint = self.chr2int[idx2chrom[idx]]
            start = self.chr_indices[chrint]
            stop = self.chr_indices[chrint + 1]
            if self.focus == 'fragments':
                # Find overlap between reads and fragments on a per
                # chromosome basis
                if count_overlaps:
                    # Count all fragments touching reads
                    starts = numpy.searchsorted(
                        self.data['coords'][start:stop, 1], seqs[idx][:, 0],
                        side='right') + start
                    stops = numpy.searchsorted(
                        self.data['coords'][start:stop, 0], seqs[idx][:, 1],
                        side='right') + start
                    for i in range(starts.shape[0]):
                        self.data[condition][starts[i]:stops[i]] += 1
                else:
                    # Count only largest overlapping fragment
                    indices = numpy.searchsorted(
                        numpy.r_[self.data['coords'][start:stop, 0],
                                 self.data['coords'][stop - 1, 1]],
                        (seqs[idx][:, 0] + seqs[idx][:, 1]) / 2,
                        side='right') - 1
                    valid = numpy.where((indices >= 0)
                                        & (indices < stop - start))[0]
                    indices += self.chr_indices[chrint]
                    for i in valid:
                        self.data[condition][indices[i]] += 1
            else:
                # Find nearest RE site
                mids = (seqs[idx][:, 0] + seqs[idx][:, 1]) / 2
                lindices = numpy.searchsorted(self.data['coords'][1:, 0] - 4,
                                              seqs[idx][:, 0])
                rindices = numpy.searchsorted(self.data['coords'][:-1, 1] + 4,
                                              seqs[idx][:, 1])
                left_dists = (seqs[idx][:, 0] + 4
                              - self.data['coords'][lindices, 0])
                right_dists = (self.data['coords'][rindices, 1]
                               - seqs[idx][:, 1] + 4)
                left = numpy.logical_or(
                    numpy.logical_and(
                        numpy.less(left_dists, right_dists),
                        left_dists >= 0),
                    right_dists < 0)
                if maxdist >= 0:
                    # Determine if read end is close enough to RE site
                    # to be counted
                    for i in numpy.where(left)[0]:
                        if left_dists[i] <= maxdist:
                            self.data[condition][lindices[i]] += 1
                    for i in numpy.where(numpy.logical_not(left))[0]:
                        if right_dists[i] <= maxdist:
                            self.data[condition][rindices[i]] += 1
                else:
                    # Count regardless of distance to RE site
                    for i in numpy.where(left)[0]:
                        self.data[condition][lindices[i]] += 1
                    for i in numpy.where(numpy.logical_not(left))[0]:
                        self.data[condition][rindices[i]] += 1
            seqs[idx] = None
        infile.close()
        # If counts are present for both treatment and control,
        # (re)calculate scores
        if (numpy.amax(self.data['treatment']) > 0
                and numpy.amax(self.data['control']) > 0):
            self.calculate_scores()

    def load_counts(self, fname, condition='treatment'):
        """Load reads from a SAM or BAM file for the specified condition."""
        if condition not in ['treatment', 'control']:
            self.logger.warn("Condition must be either 'treatment' or "
                             "'control'")
            return None
        if self.data is None:
            self.logger.warn("Restriction sites must be loaded prior to "
                             "loading reads")
            return None
        self.logger.info('Reading counts from %s' % fname)
        data = []
        with open(fname) as f:
            for line in f:
                if line[0] == '#':
                    continue
                try:
                    line = line.rstrip('\n').split('\t')
                    chrom, start, stop = line[0], int(line[1]), int(line[2])
                except (IndexError, ValueError):
                    continue
                try:
                    # First assume file is a bedgraph
                    count = int(line[3])
                except (IndexError, ValueError):
                    try:
                        # Next see if it is a bed file with valid score data
                        count = int(line[4])
                    except (IndexError, ValueError):
                        continue
                if chrom in self.chr2int:
                    data.append((self.chr2int[chrom], (start, stop), count))
        data = numpy.array(data, dtype=numpy.dtype([
            ('chr', numpy.int32), ('coords', numpy.int32, (2,)),
            ('count', numpy.int32)]))
        if data.shape[0] == 0:
            self.logger.warn("No count data could be loaded from the "
                             "specified file")
            return None
        for i in range(self.chroms.shape[0]):
            where = numpy.where(data['chr'] == i)[0]
            if where.shape[0] == 0:
                continue
            start = self.chr_indices[i]
            stop = self.chr_indices[i + 1]
            starts = numpy.searchsorted(self.data['coords'][start:stop, 0],
                                        data['coords'][where, 0]) + start
            stops = numpy.searchsorted(self.data['coords'][start:stop, 1],
                                       data['coords'][where, 1]) + start
            if (numpy.sum(self.data['coords'][starts, 0]
                          != data['coords'][where, 0])
                + numpy.sum(self.data['coords'][stops, 1]
                            != data['coords'][where, 1])) > 0:
                self.logger.warn("The counts file does not appear to use the "
                                 "same RE site mapping as was loaded into "
                                 "the Segmenter. Ignoring counts")
                return None
            self.data[condition][starts] += data['count'][where]
        # If counts are present for both treatment and control,
        # (re)calculate scores
        if (numpy.amax(self.data['treatment']) > 0
                and numpy.amax(self.data['control']) > 0):
            self.calculate_scores()

    def load_bed(self, fname, datatype='LADs'):
        if datatype not in ['LADs', 'DIPs', 'partitions']:
            self.logger.warn("Datatype must be LADs, DIPs, or partitions")
            return None
        if not os.path.exists(fname):
            self.logger.error("The bed file %s could not be found" % fname)
            return None
        data = []
        for line in open(fname):
            try:
                line = line.split('\t')
                if datatype == 'partitions':
                    chrom, start, stop, score = (line[0], int(line[1]),
                                                 int(line[2]), float(line[4]))
                    if chrom in self.chr2int:
                        chrint = self.chr2int[chrom]
                        data.append((chrint, score, (start, stop)))
                    data.append((c))
                else:
                    chrom, start, stop = (line[0], int(line[1]), int(line[2]))
                    if chrom in self.chr2int:
                        chrint = self.chr2int[chrom]
                        data.append((chrint, (start, stop)))
            except VaueError:
                pass
        if datatype == 'partitions':
            data = numpy.array(data, dtype=numpy.dtype([
                ('chr', numpy.int32), ('score', numpy.float32),
                ('coords', numpy.int32, (2,))]))
        else:
            data = numpy.array(data, dtype=numpy.dtype([
                ('chr', numpy.int32), ('coords', numpy.int32, (2,))]))
        self[datatype] = data

    def calculate_scores(self):
        if self.focus == 'binned':
            data = self.binned
        else:
            data = self.data
        valid = numpy.where(numpy.logical_not(numpy.isnan(data['score']))
                            & ((data['treatment'] > 0)
                               | (data['control'] > 0)))[0]
        treatment = numpy.maximum(
            0.5, data['treatment'][valid].astype(numpy.float32))
        control = numpy.maximum(
            0.5, data['control'][valid].astype(numpy.float32))
        data['score'].fill(numpy.nan)
        data['score'][valid] = numpy.log2((treatment / numpy.sum(treatment))
                                          / (control / numpy.sum(control)))

    def remove_unalignable(self):
        """Mark blacklisted genomic regions as unalignable in the RE data"""
        if self.unalignable is None or self.data is None:
            self.logger.warn("Both unalignable regions and RE data must be "
                             "loaded prior to running this function")
            return None
        self.logger.info("Removing unalignable regions")
        data = [self.data]
        indices = [self.chr_indices]
        if self.binned is not None:
            data.append(self.binned)
            indices.append(self.bin_indices)
        for h in range(len(data)):
            for i in range(self.chroms.shape[0]):
                if indices[h][i] == indices[h][i + 1]:
                    continue
                start = indices[h][i]
                stop = indices[h][i + 1]
                where = numpy.where(self.unalignable['chr'] == i)[0]
                if where.shape[0] == 0:
                    continue
                mids = (data[h]['coords'][start:stop, 0]
                        + data[h]['coords'][start:stop, 1]) // 2
                starts = numpy.searchsorted(
                    mids, self.unalignable['coords'][where, 0]) + start
                stops = numpy.searchsorted(
                    mids, self.unalignable['coords'][where, 1]) + start
                for j in range(starts.shape[0]):
                    data[h]['alignable'][starts[j]:stops[j]] = False
                    data[h]['score'][starts[j]:stops[j]] = numpy.nan
        # If counts are present for both treatment and control,
        # (re)calculate scores
        if (numpy.amax(self.data['treatment']) > 0
                and numpy.amax(self.data['control']) > 0):
            self.calculate_scores()

    def bin_data(self, binsize=10000):
        """Bin counts/scores into 'binsize'-sized bins."""
        if self.data is None:
            self.logger.warn("RE data must be loaded prior to binning")
            return None
        try:
            self.binsize = int(binsize)
        except ValueError:
            self.logger.warn("The binsize must be an integer")
            return None
        self.logger.info("Binning data at %i bp resolution" % self.binsize)
        self.focus = 'binned'
        bounds = numpy.zeros((self.chroms.shape[0], 2), dtype=numpy.int32)
        self.bin_indices = numpy.zeros(self.chroms.shape[0] + 1,
                                       dtype=numpy.int32)
        for i in range(self.chroms.shape[0]):
            bounds[i, 0] = (self.data['coords'][self.chr_indices[i], 0]
                            // self.binsize) * self.binsize
            bounds[i, 1] = (
                (self.data['coords'][self.chr_indices[i + 1] - 1, 1] - 1)
                // self.binsize + 1) * self.binsize
            self.bin_indices[i + 1] = (
                self.bin_indices[i] + (bounds[i, 1] - bounds[i, 0])
                // self.binsize)
        self.binned = numpy.zeros(self.bin_indices[-1], dtype=numpy.dtype([
            ('chr', numpy.int32), ('coords', numpy.int32, (2,)),
            ('treatment', numpy.int32), ('control', numpy.int32),
            ('score', numpy.float64), ('alignable', numpy.bool)]))
        self.binned['score'] = numpy.nan
        for i in range(self.chroms.shape[0]):
            start = self.chr_indices[i]
            stop = self.chr_indices[i + 1]
            bstart = self.bin_indices[i]
            bstop = self.bin_indices[i + 1]
            valid = (numpy.where(numpy.logical_not(numpy.isnan(
                self.data['score'][start:stop])))[0] + start)
            indices = (
                (self.data['coords'][valid, 0]
                 + self.data['coords'][valid, 1]) // 2
                - bounds[i, 0]) // self.binsize
            self.binned['treatment'][bstart:bstop] = numpy.bincount(
                indices, weights=self.data['treatment'][valid],
                minlength=(bstop - bstart))
            self.binned['control'][bstart:bstop] = numpy.bincount(
                indices, weights=self.data['control'][valid],
                minlength=(bstop - bstart))
            self.binned['score'][start:stop] = numpy.bincount(
                indices, weights=self.data['score'][valid],
                minlength=(stop - start))
            self.binned['chr'][start:stop] = i
            self.binned['coords'][start:stop, 0] = (
                numpy.arange(stop - start) * self.binsize + bounds[i, 0])
            self.binned['coords'][start:stop, 1] = (
                numpy.arange(stop - start) * self.binsize + bounds[i, 0]
                + self.binsize)
            self.binned['alignable'] = numpy.bincount(
                indices, minlength=(self.bin_indices[i + 1]
                                    - self.bin_indices[i])) > 0
        # If counts are present for both treatment and control,
        # (re)calculate scores
        if (numpy.amax(self.binned['treatment']) > 0
                and numpy.amax(self.binned['control']) > 0):
            self.calculate_scores()

    def segment_data(self, alpha=0.001, nperm=10000, eta=0.05, trim=0.025,
                     min_width=2, kmax=25, nmin=200, ngrid=100, tol=1e-6,
                     region=10, outlier_scale=4, smooth_scale=2, threads=1):
        """Using circular binary segmentation, find partitioning of data."""
        if self.focus != 'binned':
            if (self.data is None or numpy.nanmin(self.data['score']) == 0):
                self.logger.warn("Cannot partition data that is not loaded")
                return None
        elif (self.binned is None or numpy.nanmin(self.binned['score']) == 0):
            self.logger.warn("Cannot partition binned data that is not loaded")
            return None
        self.alpha = alpha
        self.nperm = nperm
        self.eta = eta
        self.trim = trim
        self.min_width = min_width
        self.kmax = kmax
        self.nmin = nmin
        self.ngrid = ngrid
        self.tol = tol
        self.region = region
        self.outlier_scale = outlier_scale
        self.smooth_scale = smooth_scale
        if self.sbdry is None:
            self.get_sbdry()
        smoothed = self.smooth_data(region, outlier_scale, smooth_scale, trim)
        valid = numpy.where(numpy.logical_not(numpy.isnan(smoothed)))[0]
        trimmed_SD = (self.trimmed_variance(smoothed[valid], trim)) ** 0.5
        chr_queue = multiprocessing.JoinableQueue()
        partition_queue = multiprocessing.JoinableQueue()
        partition_threads = []

        def find_changepoints(partition_queue, data, trimmed_SD, q):
            while True:
                try:
                    index = q.get(False)
                    if index is None:
                        partitions = []
                    else:
                        partitions = self.find_changepoints(data, index,
                                                            trimmed_SD)
                    partition_queue.put(partitions)
                    q.task_done()
                except Queue.Empty:
                    break

        for i in range(self.chr_indices.shape[0] - 1):
            chr_queue.put(i)
        for i in range(min(self.chroms.shape[0], threads)):
            partition_threads.append(multiprocessing.Process(
                target=find_changepoints,
                args=(partition_queue, smoothed, trimmed_SD, chr_queue)))
            partition_threads[-1].start()
        self.logger.info("Segmenting score data")
        chr_queue.join()
        partitions = []
        for i in range(self.chroms.shape[0]):
            partitions += partition_queue.get(True)
        self.partitions = numpy.array(partitions, dtype=numpy.dtype([
                ('chr', numpy.int32), ('score', numpy.float32),
                ('coords', numpy.int32, (2,))]))
        self.partitions = self.partitions[numpy.lexsort((
            self.partitions['coords'][:, 0], self.partitions['chr']))]

    def find_changepoints(self, data, index, trimmed_SD):
        """Find segmentation boundaries in a single chromosome, iteratively
        splitting each segment.
        """
        if self.focus == 'binned':
            start = self.bin_indices[index]
            stop = self.bin_indices[index + 1]
            coords = self.binned['coords']
        else:
            start = self.chr_indices[index]
            stop = self.chr_indices[index + 1]
            if self.focus == 'fragments':
                coords = self.data['coords']
            else:
                coords = numpy.zeros((self.data.shape[0], 2),
                                     dtype=numpy.int32)
                for i in range(self.chroms.shape[0]):
                    start2 = self.chr_indices[i]
                    stop2 = self.chr_indices[i + 1]
                    mids = (self.data['coords'][start2:(stop2 - 1), 1]
                            + self.data['coords'][(start2 + 1):stop2, 0]) // 2
                    coords[(start2 + 1):stop2, 0] = mids
                    coords[start2:(stop2 - 1), 1] = mids
                    coords[start2, 0] = self.data['coords'][start2, 0]
                    coords[stop2 - 1, 1] = self.data['coords'][stop2 - 1, 1]
        valid = numpy.where(
            numpy.logical_not(numpy.isnan(data[start:stop])))[0] + start
        if valid.shape[0] == 0:
            return []
        seg_ends = [0, valid.shape[0]]
        k = len(seg_ends)
        change_loc = []
        while k > 1:
            ncpt = numpy.zeros(1, numpy.int32)
            icpt = numpy.zeros(2, numpy.int32)
            current_n = seg_ends[k - 1] - seg_ends[k - 2]
            if current_n >= 2 * self.min_width:
                curr_data = data[valid[seg_ends[k - 2]:seg_ends[k - 1]]]
                curr_data -= numpy.mean(curr_data)
                if self.nmin < current_n:
                    hybrid = True
                    delta = (self.kmax + 1.) / current_n
                else:
                    hybrid = False
                    delta = 0
                if numpy.amin(curr_data) != numpy.amax(curr_data):
                    fcurrent_n = float(current_n)
                    current_tss = numpy.sum(curr_data ** 2)
                    dnacopy.fndcpt(
                        x=curr_data,
                        tss=current_tss,
                        nperm=self.nperm,
                        cpval=self.alpha,
                        ncpt=ncpt,
                        icpt=icpt,
                        ibin=False,
                        hybrid=hybrid,
                        al0=self.min_width,
                        hk=self.kmax,
                        delta=delta,
                        ngrid=self.ngrid,
                        sbdry=self.bdry,
                        tol=self.tol)
            if ncpt[0] == 0:
                change_loc.append(seg_ends[-1])
                seg_ends = seg_ends[:-1]
            elif ncpt[0] == 1:
                seg_ends = (
                    seg_ends[:-1] + [seg_ends[-2] + icpt[0], seg_ends[-1]])
            else:
                seg_ends = (
                    seg_ends[:-1] + [seg_ends[-2] + icpt[0], seg_ends[-2]
                                     + icpt[1], seg_ends[-1]])
            k = len(seg_ends)
        seg_ends = numpy.r_[0, change_loc[::-1]]
        partitions = []
        for i in range(len(seg_ends) - 1):
            start = valid[seg_ends[i]]
            stop = valid[seg_ends[i + 1] - 1] + 1
            curr_data = data[start:stop]
            mu = numpy.mean(curr_data[numpy.where(
                numpy.logical_not(numpy.isnan(curr_data)))[0]])
            partitions.append((index, mu,
                              (coords[start, 0], coords[stop - 1, 1])))
        self.logger.info("Finished segmenting chromosome %s"
                         % self.chroms[index])
        return partitions

    def smooth_data(self, region=10, outlier_scale=4, smooth_scale=2,
                    trim=0.025):
        """Perform DNAcopy-based smoothing of data."""
        if self.data is None or self.focus == 'binned' and self.binned is None:
            self.logger.warn("Cannot smooth data that is not loaded")
            return None
        self.logger.info("Smoothing score data")
        if self.focus == 'binned':
            data = self.binned
        else:
            data = self.data
        valid = numpy.where(numpy.logical_not(numpy.isnan(data['score'])))[0]
        chr_sizes = numpy.bincount(data['chr'][valid])
        chr_sizes = chr_sizes[numpy.where(chr_sizes > 0)[0]]
        smoothed = numpy.zeros(valid.shape[0], dtype=numpy.float64)
        trimmed_SD = (self.trimmed_variance(data['score'][valid], trim)) ** 0.5
        outlier_SD = outlier_scale * trimmed_SD
        smooth_SD = smooth_scale * trimmed_SD
        dnacopy.smoothlr(
            data['score'][valid], chr_sizes, smoothed, region,
            outlier_SD, smooth_SD)
        data = numpy.zeros(data.shape[0], dtype=numpy.float64)
        data.fill(numpy.nan)
        data[valid] = smoothed
        return data

    def trimmed_variance(self, data, trim):
        """Find trimmed variance of data"""
        n = int(round((1 - 2 * trim) * (data.shape[0] - 1)))
        a = norm.ppf(1 - trim)
        x = numpy.linspace(-a, a, 10001)
        x = (x[1:] + x[:-1]) / 2
        inflfact = (1 / (numpy.sum(x ** 2 * norm.pdf(x) / (1 - 2 * trim))
                    * (2 * a / 10000.)))
        diffs = numpy.abs(data[1:] - data[:-1])
        diffs = diffs[numpy.argsort(diffs)]
        return inflfact * numpy.sum((diffs[:n])**2 / (2 * n))

    def get_sbdry(self):
        """Find permuted data cutoffs."""
        self.logger.info("Finding permuted boundary cutoffs")
        n = int(numpy.floor(self.nperm * self.alpha) + 1)
        m = n * (n + 1) // 2
        self.bdry = numpy.zeros(m, dtype=numpy.int32)
        etastr = numpy.zeros(n, dtype=numpy.float64)
        dnacopy.getbdry(self.eta, self.nperm, self.bdry, etastr, self.tol)

    def find_LADs_from_segmentation(self, mindip=5000, maxdip=10000):
        """Stitch together segments into LADs and NonLADs, defining DIP
        regions in LADs.
        """
        if self.partitions is None:
            self.logger.warn("Cannot determine LADs before segmentation has "
                             "been performed")
            return None
        self.logger.info("Finding LADs and DIPs from segmentation")
        lads = []
        dips = []
        indices = numpy.r_[0, numpy.cumsum(numpy.bincount(
            self.partitions['chr'], minlength=self.chroms.shape[0]))]
        for i in range(self.chroms.shape[0]):
            start = indices[i]
            stop = indices[i + 1]
            if start == stop:
                continue
            breaks = numpy.r_[start, numpy.where(
                numpy.sign(self.partitions['score'][(start + 1):stop])
                != numpy.sign(self.partitions['score'][start:(stop - 1)])
            )[0] + 1 + start, stop]
            signs = numpy.sign(self.partitions['score'][breaks[:-1]])
            # Identify DIPs
            for j in range(1, breaks.shape[0] - 2):
                if signs[j] < 0:
                    size = (self.partitions['coords'][breaks[j], 1]
                            - self.partitions['coords'][breaks[j + 1] - 1, 0])
                    if size >= mindip and size <= maxdip:
                        dips.append(
                            (i, (self.partitions['coords'][breaks[j], 0],
                             self.partitions['coords'][breaks[j + 1] - 1, 1])))
            # Stitch together LADs
            raw = []
            for j in range(breaks.shape[0] - 1):
                if signs[j] > 0:
                    raw.append([
                        i, self.partitions['coords'][breaks[j], 0],
                        self.partitions['coords'][breaks[j + 1] - 1, 1]])

            for j in list(range(len(raw) - 1))[::-1]:
                if raw[j + 1][1] - raw[j][2] <= maxdip:
                    raw[j][2] = raw[j + 1][2]
                    del raw[j + 1]
            for lad in raw:
                lads.append((lad[0], (lad[1], lad[2])))
        lads.sort()
        dips.sort()
        self.LADs = numpy.array(lads, dtype=numpy.dtype([
            ('chr', numpy.int32), ('coords', numpy.int32, (2,))]))
        self.DIPs = numpy.array(dips, dtype=numpy.dtype([
            ('chr', numpy.int32), ('coords', numpy.int32, (2,))]))

    def write_bed(self, target, fname):
        """Write requested data to file in BED format."""
        if target not in ['LADs', 'DIPs', 'partitions', 'score',
                          'treatment', 'control']:
            self.logger.warn("Target must be LADs, DIPs, partitions, score, "
                             "treatment, or control")
            return None
        if target in ['score', 'treatment', 'control']:
            if self.focus == 'binned':
                key = 'binned'
            else:
                key = 'data'
            if self[key] is None:
                self.logger.warn('No data for %s' % target)
                return None
            scores = self[key][target]
            if target == 'score':
                valid = numpy.where(numpy.logical_not(numpy.isnan(scores)))[0]
            else:
                valid = numpy.where(scores > 0)[0]
        else:
            key = target
            if self[key] is None:
                self.logger.warn('No data for %s' % target)
                return None
            if 'score' in self[key].dtype.names:
                scores = self[key]['score']
                valid = numpy.where(numpy.logical_not(numpy.isnan(scores)))[0]
            else:
                scores = ['.'] * self[key].shape[0]
                valid = numpy.arange(self[key].shape[0])
        with open(fname, 'wb') as outfile:
            for i in valid:
                outfile.write(("%s\t%i\t%i\t.\t%s\t.\n" % (
                               self.chroms[self[key]['chr'][i]],
                               self[key]['coords'][i, 0],
                               self[key]['coords'][i, 1],
                               str(scores[i]))).encode('UTF-8'))

    def write_bedgraph(self, target, fname):
        """Write requested data to file in BEDGRAPH format."""
        if target not in ['partitions', 'score', 'treatment', 'control']:
            self.logger.warn("Target must be partitions, score, treatment, "
                             "or control")
            return None
        if target in ['score', 'treatment', 'control']:
            if self.focus == 'binned':
                key = 'binned'
            else:
                key = 'data'
            if self[key] is None:
                self.logger.warn('No data for %s' % target)
                return None
            scores = self[key][target]
            if target == 'score':
                valid = numpy.where(numpy.logical_not(numpy.isnan(scores)))[0]
            else:
                valid = numpy.where(scores > 0)[0]
        else:
            key = target
            if self[key] is None:
                self.logger.warn('No data for %s' % target)
                return None
            scores = self[key]['score']
            valid = numpy.where(numpy.logical_not(numpy.isnan(scores)))[0]
        with open(fname, 'wb') as outfile:
            for i in valid:
                    outfile.write(("%s\t%i\t%i\t%s\n" % (
                                   self.chroms[self[key]['chr'][i]],
                                   self.data['coords'][i, 0],
                                   self.data['coords'][i, 1],
                                   str(scores[i]))).encode('UTF-8'))

    def plot_scores(self, fname, binsize=10000, show_lads=False,
                    show_dips=False, show_means=False, show_partitions=False):
        """Plot a PDF of score data with optional indicators of partitions, "
        partition means, LADs, and DIPs."""
        # Make sure that the PYX module is available
        try:
            from pyx import canvas, path, document, color, text, style
        except ImportError:
            self.logger.warn("The package pyx must be installed to plot data")
            return None
        # Make sure desired data is present
        if self.focus != 'binned':
            if (self.data is None or numpy.nanmin(self.data['score']) == 0):
                self.logger.warn("Requested data is not available for "
                                 "plotting")
                return None
        elif (self.binned is None or numpy.nanmin(self.binned['score']) == 0):
            self.logger.warn("Requested binned data is not available for "
                             "plotting")
            return None
        self.logger.info("Plotting data")
        # Determine which data to use
        if self.focus == 'binned':
            data = self.binned
            chr_indices = self.bin_indices
        else:
            data = self.data
            chr_indices = self.chr_indices
        # Plot each chromosome on its own page
        pages = []
        for i in range(self.chroms.shape[0]):
            valid = numpy.where(numpy.logical_not(numpy.isnan(
                data['score'][chr_indices[i]:chr_indices[i + 1]])))[0]
            valid += chr_indices[i]
            # Skip chromosomes without valid data
            if valid.shape[0] == 0:
                continue
            mids = (data['coords'][valid, 0] + data['coords'][valid, 1]) // 2
            if binsize > 0:
                # Bin data to the resolution requested
                start = (mids[0] // binsize) * binsize
                stop = (mids[-1] // binsize + 1) * binsize
                indices = (mids - start) // binsize
                counts = numpy.bincount(indices)
                Ys = numpy.bincount(
                    indices,
                    weights=data['score'][valid]) / numpy.maximum(1, counts)
                mids = (start + numpy.arange((stop - start) // binsize)
                        * binsize + binsize / 2)
                valid = numpy.where(counts > 0)[0]
                Ys = Ys[valid]
                mids = mids[valid]
                coords = numpy.zeros((Ys.shape[0], 2), dtype=numpy.float64)
                coords[:, 0] = start + valid * binsize
                coords[:, 1] = start + valid * binsize + binsize
            else:
                Ys = data['score'][valid]
                start = data['coords'][valid, 0]
                stop = data['coords'][valid[-1], 1]
                coords = data['coords'][valid, :].astype(numpy.float64)
            c = canvas.canvas()
            width = (stop - start) / 5000000.
            height = 5.
            maxscore = numpy.amax(numpy.abs(Ys))
            if show_means and self.partitions is not None:
                where = numpy.where(self.partitions['chr'] == i)[0]
                maxscore = max(
                    maxscore,
                    numpy.amax(numpy.abs(self.partitions['score'][where])))
            Ys /= maxscore
            Ys *= height * 0.5
            Xs = (mids - start) / float(stop - start) * width
            coords -= start
            coords /= (stop - start)
            coords *= width
            lpath = path.path(path.moveto(0, 0))
            for j in range(valid.shape[0]):
                if j == 0 or valid[j] - valid[j - 1] > 1:
                    lpath.append(path.lineto(coords[j, 0], 0))
                lpath.append(path.lineto(Xs[j], Ys[j]))
                if j == Xs.shape[0] - 1 or valid[j + 1] - valid[j] > 1:
                    lpath.append(path.lineto(coords[j, 1], 0))
            lpath.append(path.lineto(width, 0))
            lpath.append(path.closepath())

            # add lads if requests and already determined
            if show_lads and self.LADs is not None:
                where = numpy.where(self.LADs['chr'] == i)[0]
                if where.shape[0] == 0:
                    continue
                for j in where:
                    X0 = ((self.LADs['coords'][j, 0] - start)
                          / float(stop - start) * width)
                    X1 = ((self.LADs['coords'][j, 1] - start)
                          / float(stop - start) * width)
                    c.fill(path.rect(X0, -height / 2, X1 - X0, height),
                           [color.gray(0.85)])

            # add dips if requests and already determined
            if show_dips and self.DIPs is not None:
                print(self.DIPs.shape)
                where = numpy.where(self.DIPs['chr'] == i)[0]
                if where.shape[0] == 0:
                    continue
                for j in where:
                    X0 = ((self.DIPs['coords'][j, 0] - start)
                          / float(stop - start) * width)
                    X1 = ((self.DIPs['coords'][j, 1] - start)
                          / float(stop - start) * width)
                    c.fill(path.rect(X0, -height / 2, X1 - X0, height),
                           [color.rgb.red])

            # add signal track
            c.fill(lpath)
            c.stroke(path.line(0, -height / 2, width, -height / 2))

            # add partition mean line if requested and already determined
            if show_means and self.partitions is not None:
                where = numpy.where(self.partitions['chr'] == i)[0]
                coords = ((self.partitions['coords'][where, :] - start)
                          / float(stop - start) * width)
                Ys = ((self.partitions['score'][where] / maxscore)
                      * height * 0.5)
                lpath = path.path(path.moveto(0, 0))
                for j in range(Ys.shape[0]):
                    if j == 0 or coords[j, 0] != coords[j - 1, 1]:
                        lpath.append(path.lineto(coords[j, 0], 0))
                    lpath.append(path.lineto(coords[j, 0], Ys[j]))
                    lpath.append(path.lineto(coords[j, 1], Ys[j]))
                    if (j == Ys.shape[0] - 1
                            or coords[j, 1] != coords[j + 1, 0]):
                        lpath.append(path.lineto(coords[j, 1], 0))
                lpath.append(path.lineto(width, 0))
                c.stroke(lpath, [color.rgb.blue, style.linewidth.THIN,
                                 color.transparency(0.5)])

            # add partition lines if requested and already determined
            if show_partitions and self.partitions is not None:
                where = numpy.where(self.partitions['chr'] == i)[0]
                coords = ((self.partitions['coords'][where, :] - start)
                          / float(stop - start) * width)
                for j in range(coords.shape[0] - 1):
                    c.stroke(
                        path.line(coords[j, 1], -height * 0.5, coords[j, 1],
                                  height * 0.5), [style.linewidth.THin])
                    if coords[j, 1] != coords[j + 1, 0]:
                        c.stroke(
                            path.line(coords[j + 1, 0], -height * 0.5,
                                      coords[j + 1, 0], height * 0.5),
                            [style.linewidth.THin])

            # add coordinates
            for j in range(int(numpy.ceil(start / 10000000.)),
                           int(numpy.floor(stop / 10000000.) + 1)):
                X = (j * 10000000. - start) / (stop - start) * width
                c.stroke(path.line(X, -height / 2, X, -height / 2 - 0.2))
                c.text(X, -height / 2 - 0.25, "%i Mb" % (j * 10),
                       [text.halign.center, text.valign.top])
            c.text(width * 0.5, -height / 2 - 0.6,
                   "%s" % self.chroms[i].replace('_', ' '),
                   [text.halign.center, text.valign.top])
            pages.append(document.page(c))
        doc = document.document(pages)
        doc.writePDFfile(fname)
