import pysam
from Cluster import *
import sys
import argparse
import os
from time import localtime, strftime
from util import verbosity
from multiprocessing import Pool


def get_file_name(input_name):
    cwd = os.getcwd()

    if '/' in input_name:
        output_name = input_name.split('/')[-1]
        output_name = '.'.join(output_name.split('.')[:-1])
    else:
        output_name = '.'.join(input_name.split('.')[:-1])

    return cwd + '/' + output_name


def get_argument_parser():
    parser = argparse.ArgumentParser(
        description='Reads BAM-file and grab alignments with barcode as a SAM-tag to analyces interactions.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('bam',
                        help='Input bam file',
                        default=False)

    parser.add_argument('-t', '--tag',
                        help='Set barcode tag in bam.',
                        default='bc',
                        type=str)

    parser.add_argument('-o', '--output',
                        help='Output file prefix',
                        default=None)

    parser.add_argument('-c', '--cluster_data',
                        help='Output file with cluster data',
                        default=False,
                        action='store_true')

    parser.add_argument('-b', '--bed',
                        help='Output bed file with cluster fragments.',
                        default=False,
                        action='store_true')

    parser.add_argument('--debug',
                        help='Debug',
                        default=False,
                        action='store_true')

    parser.add_argument('-i', '--interactions',
                        help='Include interaction analysis',
                        default=False,
                        action='store_true')

    parser.add_argument('-p', '--phasing',
                        help='Include phasing analysis',
                        default=False,
                        action='store_true')

    parser.add_argument('--stats',
                        help='Include extra cluster stats',
                        default=False,
                        action='store_true')

    parser.add_argument('--dist-only',
                        help='Include interaction analysis but only output distance information',
                        default=False,
                        action='store_true')

    parser.add_argument('--insert',
                        help='Include insert size output file.',
                        default=False,
                        action='store_true')

    parser.add_argument('-s', '--singles',
                        help='Remove single alingments in analysis',
                        default=False,
                        action='store_true')

    parser.add_argument('--separator',
                        help='Output file separator',
                        default=',',
                        type=str)

    parser.add_argument('-m', '--monochromosomal',
                        help='Only include monochromsomal clusters for phasing and interaction analysis.',
                        default=False,
                        action='store_true')

    parser.add_argument('--interchromosomal',
                        help='Include interchromosomal interactions in analysis.',
                        default=False,
                        action='store_true')

    parser.add_argument('-q', '--quiet',
                        help='No verbosity',
                        action='store_true',
                        default=False)

    return parser.parse_args()


def read_bam_file(file_name, tag, quiet, debug=False):
    verbosity('Reading and storing alignment information ' + file_name, quiet)
    file = pysam.AlignmentFile(file_name, 'rb')

    clusters = {}
    aln_count = 0
    for aln in file.fetch(until_eof=True):
        aln_count += 1

        try:
            cluster_id = aln.get_tag(tag=tag)
        except KeyError:
            verbosity('Tag ' + tag + 'not present in bam record', True)
            verbosity('Exiting program', True)
            sys.exit()

        if cluster_id not in clusters:
            clusters[cluster_id] = Cluster(cluster_id)

        clusters[cluster_id].add_pysam_alignment(aln)

        if aln_count % 1e6 == 0:
            verbosity('#' + str(aln_count) + ' alignments analyced', quiet)
            sys.stdout.flush()

        if aln_count % 1e5 == 0 and debug:
            break

    file.close()

    verbosity('Complete', quiet)

    return clusters


def main():
    global cluster_instances

    args = get_argument_parser()

    verbosity('Cluster analysis started', args.quiet)
    verbosity('COMMAND: ' + ' '.join(sys.argv), args.quiet, include_time=False)
    verbosity('OPTIONS:\n\t' + '\n\t'.join(sorted([str(option) + ': ' + str(value)
                                                   for option, value in eval(str(vars(args))).items()])),
              args.quiet, include_time=False)

    if not args.output:
        prefix = get_file_name(args.bam)
    else:
        prefix = args.output

    cluster_instances = read_bam_file(args.bam, args.tag, args.quiet, debug=args.debug)

    # Write output to file
    if not args.debug and args.cluster_data:
        # old_stdout = sys.stdout
        # log_file = open(prefix + ".cluster_info.csv", "w")
        # sys.stdout = log_file
        data_file = open(prefix + ".cluster_info.csv", "a")

    if args.bed and not args.debug:
        bed_file = open(prefix + ".cluster_fragments.bed", "a")

    verbosity('Analycing clusters', args.quiet)

    cluster_count = 0
    for cluster_id, cluster_object in cluster_instances.items():
        cluster_count += 1

        cluster_object.run_standard_analysis(debug=args.debug, debug_all=False, include_singles=not args.singles)

        if args.phasing:
            cluster_object.analyce_phasing(phasing_threshold=10000, only_monochromosomal=args.monochromosomal)

        if args.interactions:
            cluster_object.interactions_distances(only_monochromosomal=args.monochromosomal,
                                                  interchromsomal=args.interchromosomal)

        if args.stats:
            cluster_object.add_cluster_stats()

        if args.cluster_data and not args.debug:
            # cluster_object.print_csv(include_header=True, separator=args.separator)
            data_file.write(cluster_object.retrun_csv(include_header=True, separator=args.separator))

        if args.bed and not args.debug:
            # Filter for only cluster carrying fragments
            if cluster_object.counts['Fragments filtered'] >= 1:
                for chrom, fragments in cluster_object.occupied_positions_filtered.items():
                    fragments = sorted(fragments)
                    for fragment in fragments:
                        line = (chrom, str(fragment[0]), str(fragment[1]), str(cluster_id), str(fragment[3]), '+')
                        bed_file.write('\t'.join(line) + '\n')

    # Close file for saving output
    if not args.debug and args.cluster_data:
        # sys.stdout = old_stdout
        # log_file.close()
        data_file.close()

    if not args.debug and args.bed:
        bed_file.close()

    verbosity('Complete', args.quiet)

    Cluster.print_summary(pretty=False)

    # Write interactions file
    if not args.debug and args.interactions:
        if args.dist_only:
            with open(prefix + ".interaction_dists.csv", "a") as file1:
                file1.write('Index,Interaction distance\n')
                for nr, interaction in enumerate(Cluster.total_interactions_dist):
                    file1.write(str(nr) + args.separator + str(interaction) + '\n')

        else:
            with open(prefix + ".interactions.csv", "a") as file2:
                header = ['Chr1', 'Start1', 'End1', 'Weight1', 'Chr2', 'Start2', 'End2', 'Weight2', 'Dist', 'Type',
                          'Cluster_id']
                file2.write(args.separator.join(header) + '\n')
                for interaction in Cluster.total_interactions:
                    line_to_write = [str(item) for item in interaction]
                    file2.write(args.separator.join(line_to_write) + '\n')

    if not args.debug and args.insert:
        import collections
        counter = collections.Counter(Cluster.total_insert_size_filtered)
        with open(prefix + ".insert_sizes.csv", "a") as file2:
            header = ['Insert size', 'Freq']
            file2.write(args.separator.join(header) + '\n')
            for insert_size, freq in counter.items():
                file2.write(str(insert_size) + args.separator + str(freq) + '\n')

    verbosity('Analysis finished', args.quiet)


if __name__ == '__main__':
    main()
