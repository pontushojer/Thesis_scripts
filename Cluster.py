from pysam import AlignedSegment
from itertools import count
from statistics import stdev, mean


class Cluster:
    counts = dict()
    duplication_rates = []
    localization_rates = []
    total_interactions_dist = []
    total_phasing_dist = []
    total_insert_size_filtered = []
    total_insert_size_unique = []
    total_interactions = []

    interaction_threshold = 1000
    pair_threshold = 2000
    positions_merge_threshold = 5  # Possible to have 9nt overlap from tagmentation
    header_printed = False

    counts_keys = set()

    def __init__(self, cluster_id):
        Cluster.count_cluster('Clusters', 1)
        self.cluster_id = cluster_id
        self.cluster_nr = Cluster.counts['Clusters']

        self.occupied_chromosomes = set()
        self.occupied_positions = dict()
        self.occupied_positions_filtered = dict()
        self.occupied_positions_phased = dict()
        self.phased_regions = []
        self.alignments_single = dict()
        self.alignments_paired = dict()
        self.duplicates_removed = []
        self.other_removed = []
        self.interactions_dist = []
        self.interactions = []
        self.insert_sizes_filtered = []
        self.insert_sizes_unique = []
        self.phasing_dist = []
        self.read_names_for_pairing = dict()

        self.counts = {}
        self.includes = set()
        # self.is_status = set()

        self.is_overlapping = False
        self.is_filtered = False
        self.is_dedup = False
        self.is_phased = False
        self.is_localized = False
        self.includes_singles = False
        self.includes_pairs = False
        self.includes_alignments = False
        self.includes_unmapped = False
        self.includes_phasing = False
        self.includes_interactions = False
        self.includes_interactions_trans = False
        self.no_true_alignments = False

    @classmethod
    def count_cluster(cls, key, value=1):
        try:
            cls.counts[key] += value
        except KeyError:
            cls.counts[key] = value

    @classmethod
    def print_summary(cls, width=30, pretty=False):
        # Update values
        if 'Interactions trans' in cls.counts.keys():
            if cls.counts['Interactions trans'] != 0:
                Cluster.count_cluster('Cis-Trans ratio', cls.counts['Interactions']/cls.counts['Interactions trans'])

        if Cluster.total_insert_size_filtered:
            Cluster.count_cluster('Mean inserts size filtered', mean(Cluster.total_insert_size_filtered))
            Cluster.count_cluster('Std dev inserts size filtered', stdev(Cluster.total_insert_size_filtered))
        if Cluster.total_insert_size_unique:
            Cluster.count_cluster('Mean inserts size unique', mean(Cluster.total_insert_size_unique))
            Cluster.count_cluster('Std dev inserts size unique', stdev(Cluster.total_insert_size_unique))

        value_max = 10
        temp = []
        for key, value in cls.counts.items():
            if pretty:
                if len(key) > width:
                    key = key[:width-3] + '...'

                if len(str(value)) > value_max:
                    if type(value) == float:
                        value = '{0:.7}'.format(value)
                    else:
                        value = '{0:.5e}'.format(value)

                temp.append(key.ljust(width) + ':' + str(value).rjust(value_max))
            else:
                temp.append(str(key) + '\t' + str(value))
        print('\n'.join(sorted(temp)))

    def run_standard_analysis(self, debug=False, debug_all=False, include_singles=True):
        if debug:
            print('=' * 30, '\nCLUSTER = ', self.cluster_id, '\n')
            print('###\n### START ###\n###')
            self.print_status()

        if include_singles:
            self.include_single()

        self.filter_empty()

        if debug and debug_all:
            print('###\n### AFTER INCLUDING SINGLES AND FILTERING CHROMOSOMES ###\n###')
            self.print_status()

        self.dedup_alignments()

        if debug and debug_all:
            print('###\n### AFTER DEDUPLICATION ###\n###')
            self.print_status()

        self.merge_overlapping()

        if debug:
            print('###\n### END ###\n###')
            self.print_status()
            print('=' * 30)

    def count_instance(self, key, value=1):
        try:
            self.counts[key] += value
        except KeyError:
            self.counts[key] = value
        Cluster.counts_keys.add(key)

    def print_csv(self, include_header=False, separator=','):
        if include_header and not Cluster.header_printed:
            print('Cluster nr,Cluster id,' + separator.join(sorted(Cluster.counts_keys)))

        list_to_print = [str(self.counts[key]) for key in sorted(Cluster.counts_keys)]
        list_to_print = [str(self.cluster_nr), str(self.cluster_id)] + list_to_print

        print(separator.join(list_to_print))

        Cluster.header_printed = True

    def retrun_csv(self, include_header=False, separator=','):
        line = str()
        if include_header and not Cluster.header_printed:
            line = 'Cluster nr,Cluster id,' + separator.join(sorted(Cluster.counts_keys)) + '\n'

        list_to_print = [str(self.counts[key]) for key in sorted(Cluster.counts_keys)]
        list_to_print = [str(self.cluster_nr), str(self.cluster_id)] + list_to_print

        line += separator.join(list_to_print) + '\n'
        Cluster.header_printed = True
        return line

    def print_interactions(self, separator=','):
        if self.includes_interactions:
            for interaction in self.interactions:
                list_to_print = [str(item) for item in interaction]
                print(separator.join(list_to_print) + ',' + str(self.cluster_id))

    def print_status(self):
        print('INCLUDE VALUES: \n', sorted(self.includes), '\n')
        print('COUNTS: \n', sorted([(key, value) for key, value in self.counts.items()]), '\n')
        print('read_names_for_pairing:', len(self.read_names_for_pairing), '\n', self.read_names_for_pairing, '\n')
        print('duplicates_removed:', len(self.duplicates_removed), '\n', self.duplicates_removed, '\n')
        print('other_removed:', len(self.other_removed), '\n', self.other_removed, '\n')
        print('alignments_paired:', len(self.alignments_paired), '\n', self.alignments_paired, '\n')
        print('alignments_single:', len(self.alignments_single), '\n', self.alignments_single, '\n')
        print('occupied_positions:', len(self.occupied_positions), '\n', self.occupied_positions, '\n')
        print('occupied_positions_filtered:', len(self.occupied_positions_filtered), '\n',
              self.occupied_positions_filtered, '\n')

    def add_cluster_stats(self):

        if self.counts['Fragments filtered'] > 1:
            Cluster.count_cluster('Clusters not singlett', 1)
        elif self.counts['Fragments filtered'] == 1:
            Cluster.count_cluster('Clusters singlett', 1)
        else:
            Cluster.count_cluster('Clusters empty', 1)

        if self.counts['Chromosomes'] > 1:
            Cluster.count_cluster('Cluster multi chromsomes', 1)
        elif self.counts['Chromosomes'] == 1:
            Cluster.count_cluster('Cluster single chromosomes', 1)

    def add_pysam_alignment(self, aln, pair_threshold=2000, include_all_chrom=False):
        """
        Add alingment to cluster
        :param aln: Pysam fetch object.
        :param pair_threshold: Maximum distance between paired read.
        :param include_all_chrom: Remove all chromosomes whos names include 'chrUn', 'random', 'chrM', 'chrEBV'.
        """
        query_name = aln.query_name

        Cluster.count_cluster('Alignments')

        self.count_instance('Alignments')
        self.count_instance('Alignments removed', 0)
        self.count_instance('Alignments unmapped', 0)
        self.count_instance('Alignments paired', 0)
        self.count_instance('Fragments', 0)

        if aln.is_unmapped:
            # self.includes_unmapped = True
            self.includes.add('Unmapped reads')
            Cluster.count_cluster('Alignments unmapped')
            Cluster.count_cluster('Alignments', -1)

            self.count_instance('Alignments unmapped')
            self.count_instance('Alignments removed')
            self.count_instance('Alignments', -1)

            self.other_removed.append([query_name, '-', '-', '-', 'Unmapped'])
            return

        self.includes.add('Mapped alignments')

        chromosome = aln.reference_name
        start = aln.reference_start
        end = aln.reference_end

        if chromosome not in self.occupied_chromosomes:
            if not include_all_chrom and any(s in chromosome for s in ['chrUn', 'random', 'chrM', 'chrEBV']):
                Cluster.count_cluster('Alignments removed', 1)
                Cluster.count_cluster('Alignments', -1)
                self.count_instance('Alignments removed', 1)
                self.count_instance('Alignments', -1)
                self.other_removed.append([query_name, chromosome, start, end, 'Chromosome filtered'])
                return

            self.occupied_chromosomes.add(chromosome)
            self.occupied_positions[chromosome] = []
            self.occupied_positions_filtered[chromosome] = []
            self.occupied_positions_phased[chromosome] = []
            self.alignments_paired[chromosome] = dict()

            self.count_instance('Chromosomes')

        if query_name not in self.read_names_for_pairing.keys():
            self.read_names_for_pairing[query_name] = (chromosome, start, end)
        else:
            (next_chromosome, next_start, next_end) = self.read_names_for_pairing[query_name]

            if chromosome == next_chromosome:   # Required to be same chrom
                fragment_start = min([start, next_start])
                fragment_end = max([end, next_end])
                if abs(fragment_start - fragment_end) <= pair_threshold:     # Require pairs to be within threshold
                    try:
                        self.alignments_paired[chromosome][(fragment_start, fragment_end)].append(query_name)
                    except KeyError:
                        self.alignments_paired[chromosome][(fragment_start, fragment_end)] = [query_name]
                    Cluster.count_cluster('Alignments paired', 2)
                    Cluster.count_cluster('Fragments')

                    self.count_instance('Alignments paired', 2)
                    self.count_instance('Fragments')

                    self.includes.add('Paired alignments')
                    # self.includes_pairs = True
                    reason = 'Paired_with_mate'
                else:
                    Cluster.count_cluster('Alignments removed', 2)

                    self.count_instance('Alignments removed', 2)
                    self.count_instance('Alignments', -2)
                    self.includes.add('Pair separated by distance')
                    reason = 'Pair_separated'
            else:
                Cluster.count_cluster('Alignments removed', 2)

                self.count_instance('Alignments removed', 2)
                self.count_instance('Alignments', -2)
                self.includes.add('Pair separated by chromosome')
                reason = 'Pair_dif_chrom'

            del self.read_names_for_pairing[query_name]

            if self.read_names_for_pairing:
                self.includes.add('Singles')
            elif 'Singles' in self.includes:
                self.includes.remove('Singles')

            self.other_removed.append([query_name, chromosome, start, end, reason])
            self.other_removed.append([query_name, next_chromosome, next_start, next_end, reason])

    def filter_empty(self):
        chrom_to_remove = []
        for chromosome in self.occupied_chromosomes:
            in_paired = False
            in_single = False

            if chromosome in self.alignments_paired.keys() and not self.alignments_paired[chromosome]:
                del self.alignments_paired[chromosome]
                in_paired = True

            if 'Singles' in self.includes:
                if chromosome in self.alignments_single.keys() and not self.alignments_single[chromosome]:
                    del self.alignments_single[chromosome]
                    in_single = True
                elif chromosome in self.alignments_single.keys():
                    continue

            if in_paired or in_single:
                chrom_to_remove.append(chromosome)

        self.count_instance('Chromosomes', -len(chrom_to_remove))

        for chromosome in chrom_to_remove:
            self.occupied_chromosomes.remove(chromosome)
            del self.occupied_positions[chromosome]
            del self.occupied_positions_filtered[chromosome]
            del self.occupied_positions_phased[chromosome]

    def include_single(self):
        self.count_instance('Alignments single', 0)

        if not self.read_names_for_pairing:
            return

        else:
            # self.includes_singles = True
            self.includes.add('Singles')

        for query_name, (chromosome, start, end) in self.read_names_for_pairing.items():
            Cluster.count_cluster('Alignments single', 1)
            Cluster.count_cluster('Fragments', 1)

            self.count_instance('Alignments single')
            self.count_instance('Fragments')
            if chromosome not in self.alignments_single:
                self.alignments_single[chromosome] = {'start': {}, 'end': {}}

            try:
                self.alignments_single[chromosome]['start'][start].append((query_name, end))
                self.alignments_single[chromosome]['end'][end].append((query_name, start))
            except KeyError:
                self.alignments_single[chromosome]['start'][start] = [(query_name, end)]
                self.alignments_single[chromosome]['end'][end] = [(query_name, start)]

    def dedup_alignments(self, include_single=True, insert_size_analysis=False):
        self.count_instance('Duplicates', 0)
        self.count_instance('Fragments unique', 0)

        if not any(type_aln in self.includes for type_aln in ['Singles', 'Paired alignments']):
            self.includes.add('No alignments')
            return

        # PAIRS
        if 'Paired alignments' in self.includes:
            for chromosome in self.alignments_paired.keys():
                if not self.alignments_paired[chromosome]:
                    continue

                for (start, end) in self.alignments_paired[chromosome].keys():
                    position_alignments = self.alignments_paired[chromosome][(start, end)]
                    weight = len(position_alignments)

                    Cluster.count_cluster('Fragments unique')

                    self.count_instance('Fragments unique')

                    if weight > 1:
                        self.duplicates_removed += [(name, chromosome, start, end, 'Pair')
                                                    for name in position_alignments[1:]]
                        Cluster.count_cluster('Duplicates', weight - 1)
                        self.count_instance('Duplicates', weight - 1)

                    self.occupied_positions[chromosome].append((start, end, position_alignments[0], weight))

                    # Remove any overlaping singletts
                    # if self.includes_singles and include_single and chromosome in self.alignments_single:
                    #     if start in self.alignments_single[chromosome]['start'].keys():
                    #         single_alignments = [(name, chromosome, start, single_end, 'Single') for (name, single_end)
                    #                              in self.alignments_single[chromosome]['start'][start]]
                    #
                    #         self.duplicates_removed += single_alignments
                    #
                    #         Cluster.count_cluster('Duplicates', len(single_alignments))
                    #
                    #         self.count_instance('Duplicates', len(single_alignments))
                    #         self.total_duplicates += len(single_alignments)
                    #
                    #         for single in single_alignments:
                    #             if single[3] in self.alignments_single[chromosome]['end']:
                    #                 del self.alignments_single[chromosome]['end'][single[3]]
                    #
                    #         del self.alignments_single[chromosome]['start'][start]
                    #
                    #     elif end in self.alignments_single[chromosome]['end'].keys():
                    #         single_alignments = [(name, chromosome, single_start, end, 'Single') for
                    #                              (name, single_start) in self.alignments_single[chromosome]['end'][end]]
                    #
                    #         self.duplicates_removed += single_alignments
                    #
                    #         Cluster.count_cluster('Duplicates', len(single_alignments))
                    #
                    #         self.count_instance('Duplicates', len(single_alignments))
                    #         self.total_duplicates += len(single_alignments)
                    #
                    #         for single in single_alignments:
                    #             if single[2] in self.alignments_single[chromosome]['start']:
                    #                 del self.alignments_single[chromosome]['start'][single[2]]
                    #
                    #         del self.alignments_single[chromosome]['end'][end]

        # SINGLES
        if 'Singles' in self.includes and include_single:
            for chromosome in self.alignments_single.keys():
                for start in self.alignments_single[chromosome]['start'].keys():
                    position_alignments = self.alignments_single[chromosome]['start'][start]
                    weight = len(position_alignments)

                    Cluster.count_cluster('Fragments unique', 1)
                    self.count_instance('Fragments unique')

                    if weight > 1:
                        self.duplicates_removed += [(name, chromosome, start, end, 'Single') for name, end in
                                                    position_alignments[1:]]
                        Cluster.count_cluster('Duplicates', weight - 1)
                        self.count_instance('Duplicates', weight - 1)

                        self.includes.add('Duplicates on position')

                    self.occupied_positions[chromosome].append((start, position_alignments[0][1],
                                                                position_alignments[0][0], weight))

                    for single in position_alignments:
                        if single[1] in self.alignments_single[chromosome]['end']:
                            del self.alignments_single[chromosome]['end'][single[1]]

        self.is_dedup = True

        if insert_size_analysis:
            self.insert_sizes_unique = [int(p[1] - p[0]) for positions in self.occupied_positions.values()
                                          for p in positions]
            Cluster.total_insert_size_unique += self.insert_sizes_unique

            if len(self.insert_sizes_unique) >= 2:
                self.count_instance('Mean inserts size unique', mean(self.insert_sizes_unique))
                self.count_instance('Std dev inserts size unique', stdev(self.insert_sizes_unique))
            elif self.insert_sizes_unique:
                self.count_instance('Mean inserts size unique', self.insert_sizes_unique[0])
                self.count_instance('Std dev inserts size unique', 0)

    def merge_overlapping(self, merge_threshold=5, insert_size_analysis=False):
        self.count_instance('Fragments filtered', 0)

        if not self.is_dedup:
            return

        for chromosome in self.occupied_positions.keys():
            position_list = sorted(self.occupied_positions[chromosome])

            if len(position_list) > 1:
                # First element
                (prev_start, prev_end, prev_name, prev_weight) = position_list[0]

                for (start, end, name, weight) in position_list[1:]:
                    # Check if any of ends of fragment overlap within threshold
                    if abs(prev_start - start) < merge_threshold or abs(prev_end - end) < merge_threshold:
                        # if prev_start != start and prev_end != end:
                        #     print(self.cluster_id, chromosome, prev_start, prev_end, start, end)
                    # Check if fragments overlap
                    # if prev_end > start and (start - prev_end) >= merge_threshold:
                    # if (start - prev_end) >= merge_threshold:
                    # if prev_end > start:                                      # USED UNTIL 180114
                        Cluster.count_cluster('Duplicates')
                        self.count_instance('Duplicates')
                        self.includes.add('Overlapping fragments')
                        self.duplicates_removed += [prev_name, chromosome, prev_start, prev_end, 'Merge']
                        self.duplicates_removed += [name, chromosome, start, end, 'Merge']

                        prev_end = max([prev_end, end])
                        prev_name = 'Merged'
                        prev_weight += weight
                        # If overlap, check if current position had more duplicates i.e. is true.
                        # if weight > prev_weight:
                        #     self.duplicates_removed += [prev_name, chromosome, prev_start, prev_end, 'Merge']
                        #     (prev_start, prev_end, prev_name, prev_weight) = (start, end, name, weight)

                        # Else, keep previous position
                        # else:
                        #     self.duplicates_removed += [name, chromosome, start, end, 'Merge']

                    # If there is no overlap, add the previous to filter position and asign new previous
                    else:
                        Cluster.count_cluster('Fragments filtered')
                        self.count_instance('Fragments filtered')

                        self.occupied_positions_filtered[chromosome].append((prev_start, prev_end, prev_name,
                                                                             prev_weight))
                        (prev_start, prev_end, prev_name, prev_weight) = (start, end, name, weight)

                # Last element shoould always be added
                Cluster.count_cluster('Fragments filtered')
                self.count_instance('Fragments filtered')

                self.occupied_positions_filtered[chromosome].append((prev_start, prev_end, prev_name, prev_weight))

            # If only one fragment on selected chromosome.
            elif position_list:
                Cluster.count_cluster('Fragments filtered')
                self.count_instance('Fragments filtered')

                self.occupied_positions_filtered[chromosome].append(position_list[0])

        self.is_filtered = True

        if insert_size_analysis:
            self.insert_sizes_filtered = [int(p[1]-p[0]) for positions in self.occupied_positions_filtered.values()
                                          for p in positions]
            Cluster.total_insert_size_filtered += self.insert_sizes_filtered

            if len(self.insert_sizes_filtered) >= 2:
                self.count_instance('Mean inserts size filtered', mean(self.insert_sizes_filtered))
                self.count_instance('Std dev inserts size filtered', stdev(self.insert_sizes_filtered))
            elif self.insert_sizes_filtered:
                self.count_instance('Mean inserts size filtered', self.insert_sizes_filtered[0])
                self.count_instance('Std dev inserts size filtered', 0)

    def calculate_duplication(self):
        self.count_instance('Duplication rate', 1)

        if not self.is_dedup and not self.is_filtered:
            return

        elif self.is_filtered:
            self.counts['Duplication rate'] = self.counts['Duplicates'] / self.counts['Fragments']

            Cluster.duplication_rates.append(self.counts['Duplication rate'])

        else:
            self.counts['Duplication rate'] = self.counts['Duplicates'] / self.counts['Fragments']

            Cluster.duplication_rates.append(self.counts['Duplication rate'])

        Cluster.counts['Mean duplication rate'] = sum(Cluster.duplication_rates) / len(Cluster.duplication_rates)

    def analyce_phasing(self, phasing_threshold=100000, only_monochromosomal=False):
        self.count_instance('Phased regions', 0)

        if not self.is_filtered:
            return

        if self.counts['Chromosomes'] > 1 and only_monochromosomal:
            return

        for chromosome in self.occupied_positions_filtered.keys():
            position_list = sorted(self.occupied_positions_filtered[chromosome])

            if len(position_list) == 1:
                self.occupied_positions_phased[chromosome].append(position_list[0] + tuple([1]))
                continue
            elif not position_list:
                continue

            (prev_start, prev_end, prev_name, prev_weight) = position_list[0]
            positions_included = 1
            region_positions = []
            for (start, end, name, weight) in position_list[1:]:
                # Check if fragments overlap within threshold
                if abs(prev_end - start) < phasing_threshold:
                    if region_positions:
                        region_positions.append((start, end, name, weight))
                    else:
                        region_positions = [(prev_start, prev_end, prev_name, prev_weight),
                                            (start, end, name, weight)]
                    prev_end = end
                    prev_weight += weight
                    positions_included += 1

                else:
                    self.occupied_positions_phased[chromosome].append([prev_start, prev_end, prev_name, prev_weight,
                                                                       positions_included])
                    if positions_included > 1:
                        self.phasing_dist.append(abs(prev_end-prev_start))
                        Cluster.total_phasing_dist.append(abs(prev_end-prev_start))

                        self.includes_phasing = True
                        self.includes.add('Phasing')

                        Cluster.count_cluster('Phased regions', 1)
                        self.count_instance('Phased regions')

                        self.phased_regions.append([[chromosome, prev_start, prev_end, prev_name, prev_weight,
                                                     positions_included, region_positions]])

                    (prev_start, prev_end, prev_name, prev_weight) = (start, end, name, weight)
                    positions_included = 1
                    region_positions = []

            self.occupied_positions_phased[chromosome].append([prev_start, prev_end, prev_name, prev_weight,
                                                               positions_included])
            if positions_included > 1:
                self.phasing_dist.append(abs(prev_end - prev_start))
                Cluster.total_phasing_dist.append(abs(prev_end - prev_start))

                self.includes_phasing = True

                Cluster.count_cluster('Phased regions', 1)
                self.count_instance('Phased regions')

                self.phased_regions.append([[chromosome, prev_start, prev_end, prev_name, prev_weight,
                                             positions_included, region_positions]])
        self.is_phased = True
        if self.includes_phasing:
            Cluster.count_cluster('Cluster with phasing', 1)

    def interactions_distances(self, interaction_threshold=interaction_threshold, only_monochromosomal=True,
                               interchromsomal=False):
        self.count_instance('Interactions', 0)
        self.count_instance('Interactions proximate', 0)

        if interchromsomal:
            trans_chroms = []
            self.count_instance('Interactions trans', 0)

        if not self.is_filtered:
            return

        if only_monochromosomal and len(self.occupied_chromosomes) > 1:
            return

        for chromosome in self.occupied_positions_filtered.keys():

            position_list = sorted(self.occupied_positions_filtered[chromosome])

            # TO DO
            #   Add filtering of proximal loci to not call distal interactions to this region more than once.
            #

            if interchromsomal and len(self.occupied_chromosomes) > 1:
                self.includes_interactions_trans = True
                self.includes.add('Interactions trans')
                trans_list = []
                for trans_chrom in self.occupied_positions_filtered.keys():
                    if trans_chrom != chromosome and trans_chrom not in trans_chroms:
                        trans_position_list = sorted(self.occupied_positions_filtered[trans_chrom])

                        for pos in trans_position_list:
                            trans_list.append([trans_chrom, pos[0], pos[1], pos[3]])

            while len(position_list) > 1:
                p1 = position_list.pop()

                for p2 in position_list:
                    dist = abs(p1[1] - p2[0])
                    if dist > interaction_threshold:
                        self.interactions_dist.append(dist)
                        # Interactions: pos1_chr pos1_start pos1_end pos1_weight pos2_chr pos2_start pos2_end
                        #               pos2_weight distance type

                        # Frist direction
                        self.interactions.append([chromosome, p1[0], p1[1], p1[3], chromosome, p2[0], p2[1], p2[3],
                                                  dist, 'cis'])
                        # Second direction
                        self.interactions.append([chromosome, p2[0], p2[1], p2[3], chromosome, p1[0], p1[1], p1[3],
                                                  dist, 'cis'])

                        Cluster.total_interactions_dist.append(dist)
                        self.includes_interactions = True
                        self.includes.add('Interactions')

                        Cluster.count_cluster('Interactions', 1)
                        self.count_instance('Interactions')
                    else:
                        Cluster.count_cluster('Interactions proximate', 1)
                        self.count_instance('Interactions proximate')

                if self.includes_interactions_trans:
                    # print(len(position_list), position_list)
                    for pos in trans_list:
                        self.interactions.append([chromosome, p1[0], p1[1], p1[3]] + pos + ['-', 'trans'])
                        Cluster.count_cluster('Interactions trans', 0.5)
                        self.count_instance('Interactions trans')

            if self.includes_interactions_trans:
                p1 = position_list.pop()
                # print(len(position_list), position_list)
                for pos in trans_list:
                    self.interactions.append([chromosome, p1[0], p1[1], p1[3]] + pos + ['-', 'trans'])
                    Cluster.count_cluster('Interactions trans', 0.5)
                    self.count_instance('Interactions trans')

                trans_chroms.append(chromosome)

        if self.includes_interactions:
            Cluster.count_cluster('Cluster with interactions', 1)
            Cluster.total_interactions += [i + [self.cluster_id] for i in self.interactions]

        if self.includes_interactions_trans:
            Cluster.count_cluster('Cluster with interactions trans', 1)

    def get_localization(self):
        if not self.is_filtered:
            self.count_instance('Localization rate', 0)
            return

        # If all reads located to one chromosome --> totally localized
        if self.counts['Chromosomes'] == 1:
            self.count_instance('Localization rate', 1)

        else:
            chromosome_freq = []
            for chromosome in self.occupied_positions_filtered.keys():
                postion_list = self.occupied_positions_filtered[chromosome]
                chromosome_freq.append(len(postion_list))

            # If on multiple chromosomes, the one with the most alignments is considered local and
            # the rate calcluted by dividing with the total number of positions.
            self.count_instance('Localization rate', max(chromosome_freq) / sum(chromosome_freq))

        Cluster.localization_rates.append(self.counts['Localization rate'])
        Cluster.counts['Mean localization rate'] = sum(Cluster.localization_rates) / len(Cluster.localization_rates)

        if self.counts['Localization rate'] != 0:
            self.is_localized = True


