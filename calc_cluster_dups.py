import sys
import pysam
import argparse
from itertools import combinations
import networkx as nx
from networkx.algorithms import components
from util import verbosity


def get_argument_parser():
    parser = argparse.ArgumentParser(
        description='Reads BAM-file and calc cluster duplicates.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('bam', help='Input bam file', default=False)
    parser.add_argument('-t', '--tag', help='Set barcode tag in bam.', default='bc', type=str)
    parser.add_argument('-s', '--separator', help='Set separator if read name contains pair informations or else to be '
                                                  'removed', type=str, default=None)
    parser.add_argument('-p', '--positions-threshold',
                        help='Number of cluster overlapping positions to required to call duplicate.',
                        type=int, default=2)
    parser.add_argument('--translate', help='Translate cluster duplicates to parent cluster', default=False,
                        action='store_true')
    parser.add_argument('--quiet', help='No verbosity', default=False, action='store_true')
    return parser.parse_args()


def resolve_duplicate_clusters(clusters, pairs):
    if not clusters:
        summary.add('Duplicate components', 0)
        summary.add('Clusters truely duplicate', 0)
        return dict(), dict()

    G = nx.Graph()

    G.add_nodes_from(clusters)
    G.add_edges_from(pairs)
    #
    # print('Nodes:', G.number_of_nodes())
    # print('Edges:', G.number_of_edges())

    components_list = []
    if not components.is_connected(G):
        # verbosity('Graph not connected, components = ' + str(components.number_connected_components(G)), args.quiet)
        for component in components.connected_components(G):
            components_list.append(component)
    else:
        components_list.append(clusters)

    summary.add('Duplicate components', components.number_connected_components(G))

    translation_dict = {}
    duplicate_weights = {}
    for component in components_list:
        main_node = min(component)
        component.remove(main_node)
        try:
            duplicate_weights[len(component)] += 1
        except KeyError:
            duplicate_weights[len(component)] = 1

        for node in component:
            translation_dict[node] = main_node
            summary.add('Clusters truely duplicate')

    return translation_dict, duplicate_weights


def main():
    global summary, args
    args = get_argument_parser()

    verbosity('Starting analysis', args.quiet)

    verbosity('COMMAND: ' + ' '.join(sys.argv), args.quiet, include_time=False)

    verbosity('OPTIONS:\n\t' + '\n\t'.join(sorted([str(option) + ': ' + str(value)
                                                   for option, value in eval(str(vars(args))).items()])),
              args.quiet, include_time=False)

    verbosity('Parsing file', args.quiet)

    file_name = args.bam
    tag = args.tag

    bamfile = pysam.AlignmentFile(file_name, 'rb')

    positions = dict()
    dup_positions = set()

    summary = Summary()

    clusters = set()
    read1_cache = dict()
    read2_cache = dict()
    sep = args.separator
    for aln in bamfile.fetch(until_eof=True):
        summary.add('Alignments total')

        cluster = int(aln.get_tag(tag))
        clusters.add(cluster)

        if not aln.is_proper_pair:
            summary.add('Alignments not pp')
            continue

        summary.add('Alignments pp')

        if sep:
            name = aln.qname.split(sep)[0]
        else:
            name = aln.qname

        aln_mate = None
        if aln.is_read1:
            if name in read2_cache:
                aln_mate = read2_cache[name]
                del read2_cache[name]
            else:
                read1_cache[name] = aln
                continue

        if aln.is_read2:
            if name in read1_cache:
                aln_mate = read1_cache[name]
                del read1_cache[name]
            else:
                read2_cache[name] = aln
                continue

        if aln_mate:
            summary.add('PP')
            chrom = aln.reference_name
            start = min(aln.reference_start, aln_mate.reference_start)
            end = max(aln.reference_start + aln.query_length, aln_mate.reference_start + aln_mate.query_length)
            pos = (chrom, start, end)
            try:
                positions[pos].add(cluster)
            except KeyError:
                positions[pos] = {cluster}

            if aln.is_duplicate and aln_mate.is_duplicate:
                dup_positions.add(pos)
                summary.add('PP duplicate')

    bamfile.close()

    summary.add('Positions', len(positions))
    summary.add('Positions duplicate', len(dup_positions))
    summary.add('Clusters', len(clusters))

    verbosity('Parsing duplicate positions', args.quiet)

    cluster_pair_dups = dict()
    clusters_dup = set()
    for dup_pos in dup_positions:
        clusters_at_pos = positions[dup_pos]
        clusters_at_pos = sorted(list(clusters_at_pos))

        clusters_dup.update(clusters_at_pos)

        cluster_pairs = combinations(clusters_at_pos, 2)

        for c1, c2 in cluster_pairs:
            try:
                cluster_pair_dups[(c1, c2)] += 1
            except KeyError:
                cluster_pair_dups[(c1, c2)] = 1
                summary.add('Cluster pairs on dupl pos')

    summary.add('Clusters on dupl pos', len(clusters_dup))

    pairs_above_threshold = list()
    clusters_above_threshold = set()
    dup_counts = dict()

    for cluster_pair, dup_count in cluster_pair_dups.items():
        if dup_count >= args.positions_threshold:
            summary.add('Cluster pairs above treshold')
            pairs_above_threshold.append(cluster_pair)
            clusters_above_threshold.update(set(cluster_pair))
            # print(cluster_pair, dup_count)

        try:
            dup_counts[dup_count] += 1
        except KeyError:
            dup_counts[dup_count] = 1

    translation_dict, dup_weights = resolve_duplicate_clusters(clusters_above_threshold, pairs_above_threshold)

    if args.translate:
        bamfile = pysam.AlignmentFile(file_name, 'rb')
        outfile = pysam.AlignmentFile('.'.join(file_name.split('.')[:-1]) + '.cluster_rmdup.bam', 'wb',
                                      template=bamfile)

        for aln in bamfile.fetch(until_eof=True):
            cluster = int(aln.get_tag(args.tag))
            chrom = aln.reference_name
            start = min(aln.reference_start, aln_mate.reference_start)
            end = max(aln.reference_start + aln.query_length, aln_mate.reference_start + aln_mate.query_length)
            pos = (chrom, start, end)

            try:
                parent_cluster = translation_dict[cluster]
                aln.set_tag('BC', str(cluster), value_type='Z')
                aln.set_tag(args.tag, str(parent_cluster), value_type='Z')
                summary.add('Alignments translated')
            except KeyError:
                pass
            outfile.write(aln)
        outfile.close()
        bamfile.close()

    summary.print()

    # print('\nDUP COUNTS')
    # print('\n'.join(sorted([str(key) + '\t' + str(value) for key, value in dup_counts.items()],
    #                        key=lambda x: int(x.split('\t')[0]))))
    print('\nDUP COUNTS')
    print('\n'.join(sorted([str(key) + '\t' + str(value) for key, value in dup_weights.items()],
                           key=lambda x: int(x.split('\t')[0]))))


class Summary:
    def __init__(self):
        self.count = dict()

    def add(self, key, value=1):
        try:
            self.count[key] += value
        except KeyError:
            self.count[key] = value

    def print(self):
        print('\n'.join(sorted([key + '\t' + str(value) for key, value in self.count.items()])))


if __name__ == '__main__':
    main()
