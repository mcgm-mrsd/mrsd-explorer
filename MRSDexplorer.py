# v6.0 overhauls the analysis to use the tabulated format of input data

import argparse
import gzip
import math
import os
import pickle
import sys


def parse_args():

    is_flagged = False
    splice_dir = os.path.realpath(sys.path[0])
    splice_dir = os.path.join(splice_dir, "..", "files")

    parser = argparse.ArgumentParser(description="Script to generate minimum required sequencing depth (MRSD) " +
                                     "predictions for transcripts, genes and tissues of interest.")
    params = parser.add_argument_group("Model parameters")

    parser.add_argument("gene_list", nargs="+", help="Path to file(s) containing list of GENCODE v19 IDs," +
                                                     "or BED format file with exon boundaries (comma-separated)")
    parser.add_argument("--tissues", nargs="+", default="all",
                        help="Prefix of input files for tissue(s) of interest in splice_dir (see below)")
    parser.add_argument("--output_prefix", action="store", default="MRSD_predictions",
                        help="Prefix for output file names.")
    parser.add_argument("--splice_dir", default=splice_dir,
                        help="Path to directory containing control files, if not in 'files' folder of this tool")
    parser.add_argument("--output_junctions", action="store_true", default=False,
                        help="Include output file with per-1 M read coverage values for individual junctions")

    parser.add_argument("--transcript_model", nargs=1, default=os.path.join(splice_dir,
                                                                            "genome-wide.transcripts.apr20.pkl"),
                        help="Path to dictionary of transcript model.")

    params.add_argument("--sj_prop", action="store", type=float, default=0.75,
                        help="Desired proportion of transcript splice junctions to be covered by number_reads reads")
    params.add_argument("--number_reads", action="store", type=int, default=8,
                        help="Desired number of reads to cover sj_prop of splice junctions in each transcript")
    params.add_argument("--stringency", action="store", type=float, default=0.95,
                        help="Percentile value to be returned from list of all predictions")
    params.add_argument("--read_type", action="store", choices=["total", "unique", "multimap"], default="unique",
                        help="Type of read to be used for MRSD calculation")
    params.add_argument("--dp", action="store", type=int, default=2,
                        help="Number of decimal places used when returning MRSD predictions")

    args = parser.parse_args()

    args.tissues = check_tissues(args.tissues, args.splice_dir)

    if is_flagged:
        sys.exit()

    return args


def check_tissues(tissue_selection, splice_dir):

    """ Ensures that the specified directory (args.splice_dir) contains valid input files for the selected tissues. """

    tissues = {}

    if tissue_selection != "all":
        tissues = {x: [False, False] for x in tissue_selection[0].split(',')}

    for file in os.listdir(splice_dir):
        pfx = file.split('.')[0]

        if pfx not in tissues and tissue_selection == "all":
            tissues[pfx] = [False, False]
        if pfx in tissues:
            if file == pfx + '.junc-counts.txt' or file == pfx + '.junc-counts.txt.gz':
                tissues[pfx][0] = True
            if file == pfx + '.seq-depths.txt':
                tissues[pfx][1] = True

    tissue_success = [[], [], []]
    for tissue in tissues:
        if all(tissues[tissue]):
            tissue_success[0].append(tissue)
        elif any(tissues[tissue]):
            tissue_success[1].append(tissue)
        else:
            tissue_success[2].append(tissue)

    if tissue_success[0]:
        if tissue_success[1] or tissue_success[2]:
            print('The following tissues were successfully found: ' + ', '.join(tissue_success[0]))
        else:
            print('All tissues successfully loaded: (' + ', '.join(tissue_success[0]) + ')\n')
    if tissue_success[1]:
        print('Only one correctly named input file was identified for the following tissues: '
              + ', '.join(tissue_success[1]))
    if tissue_success[2]:
        print('No correctly named input file was identified for the following tissues: ' + ', '.join(tissue_success[2]))
    if tissue_success[1] or tissue_success[2]:
        print('Ensure each desired tissue type has two files in the ' + splice_dir + ' folder with the names ' +
              'PREFIX.junc-counts.txt and PREFIX.seq-depths.txt, where PREFIX is the name of the tissue. Consult' +
              ' the documentation for more information on how to generate these files.')

    if not tissue_success[0]:
        sys.exit('At least one valid pair of tissue files required for analysis.')

    return tissue_success[0]


def pickle_load(file_name):

    """ Loads requested serialised file using Python pickling protocols. """

    with open(file_name, 'rb') as file:
        return pickle.load(file)


def collate_read_depths(tissue, read_depths, read_idx):

    """ Store read depths read from PREFIX.seq-depths.txt """

    with open(os.path.join(args.splice_dir, tissue+".seq-depths.txt")) as ipt:

        depth_dict = {tissue: {}}

        for line in ipt:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            depth_dict[tissue][line[0]] = line[read_idx]

        read_depths.update(depth_dict)
    return read_depths


def parse_input_file_type(input_file):

    """Identifies whether input file is in gene list or BED format"""

    with open(input_file) as ipt:
        for line in ipt:
            line = line.strip().split('\t')
            if len(line) <= 2:
                return "list"
            elif len(line) >= 4:
                return "bed"
            else:
                sys.exit("Error in input file " + input_file + "\nInvalid number of columns. Files should have" +
                         " maximum 2 columns if in gene ID format, or minimum 4 if in BED format.")


def collect_list_sjs(tpt_dict, sj_dict, list_file, read_depths, transcript_model, failed_genes):

    """Creates list of genes and splice junctions of interest based on .txt format entry."""

    match_count = 0

    for line in list_file:

        line = [x.strip() for x in line.strip().split('\t')]
        matched = False

        for item in line:
            if item in transcript_model:

                gene_match = True

                chrm = transcript_model[item][0].replace('chr', '')

                for intron in transcript_model[item][3:]:
                    sj_id = '-'.join([chrm, intron[0], intron[1]])

                    for tissue in read_depths:
                        if line[0] not in tpt_dict:
                            tpt_dict[line[0]] = {tissue: {'ids': line, sj_id: []}}
                        else:
                            if tissue not in tpt_dict[line[0]]:
                                tpt_dict[line[0]][tissue] = {'ids': line, sj_id: []}
                            else:
                                tpt_dict[line[0]][tissue][sj_id] = []

                        samp_no = len(read_depths[tissue])
                        if sj_id not in sj_dict:
                            sj_dict[sj_id] = {}
                        sj_dict[sj_id][tissue] = [0] * samp_no

                matched = True
                break

        if not matched and '/'.join(line) not in failed_genes:
            failed_genes.append('/'.join(line))
        else:
            match_count += 1

    print('\n{}/{} ({}%) gene IDs successfully mapped to transcripts.\n'.format(str(match_count),
                                                                                str(match_count + len(failed_genes)),
                                                                                round(float(match_count)*100 /
                                                                                float(match_count +
                                                                                      len(failed_genes)), 2)))

    print('The following gene IDs cannot be found in our master splice junction set: ' + '; '.join(failed_genes) + '.')
    print('These may be single-exon genes, this gene name may be deprecated, or may be differently' +
          ' named in this genome version (GENCODE v19):\n')

    return tpt_dict, sj_dict, failed_genes


def collect_bed_sjs(tpt_dict, sj_dict, bed, read_depths):
    """Creates list of genes and splice junctions of interest based on BED format entry."""
    for line in bed:

        if line.startswith('#'):
            continue

        line = line.strip().split('\t')
        sj_id = '-'.join((line[0].replace('chr', ''), line[1], line[2]))

        for tissue in read_depths:
            if line[3] not in tpt_dict:
                tpt_dict[line[3]] = {tissue: {'ids': [line[3]], sj_id: []}}
            else:
                tpt_dict[line[3]][tissue][sj_id] = []

            samp_no = len(read_depths[tissue])
            if sj_id not in sj_dict:
                sj_dict[sj_id] = {}
            sj_dict[sj_id][tissue] = [0] * samp_no

    return tpt_dict, sj_dict


def check_id_concordance(read_depths, id_dict, tissue):

    """ Check list of IDs from both read depth and SJ input files and check IDs match up. If not, remove troublesome
    entries and, if necessary, return indices to be skipped when reading master SJ file."""

    depth_ids = set(read_depths[tissue].keys())
    sj_ids = set(id_dict[tissue])

    singleton_depth_ids = depth_ids.difference(sj_ids)
    singleton_sj_ids = sj_ids.difference(depth_ids)
    exc_idxs = []

    if singleton_depth_ids:
        print('The following IDs are present in your seq-depths file but absent from your junc-counts file: ' +
              ', '.join(list(singleton_depth_ids)))
        for sample_id in singleton_depth_ids:
            del read_depths[tissue][sample_id]

    if singleton_sj_ids:
        print('The following IDs are present in your junc-counts file but absent from your seq-depths file: ' +
              ', '.join(list(singleton_sj_ids)))
        exc_idxs = [id_dict[tissue].index(x) for x in singleton_sj_ids]
        id_dict[tissue] = [x for x in id_dict[tissue] if x not in exc_idxs]

    if singleton_sj_ids or singleton_depth_ids:
        print('Removing from analysis. Consider manually cross-referencing IDs between input files.\n')

    print(str(len(read_depths[tissue])) + ' sample IDs successfully loaded for ' + tissue + ' dataset.')

    return read_depths, id_dict, exc_idxs


def collate_master_sjs(read_depths, id_dict, sj_dict, tissue):

    """Reads in master splice junction file and collates junction-supporting read counts for all control samples."""

    try:
        if tissue + ".junc-counts.txt.gz" in os.listdir(args.splice_dir):
            master_sj_file = gzip.open(os.path.join(args.splice_dir, tissue + ".junc-counts.txt.gz"), "rt")
        else:
            master_sj_file = open(os.path.join(args.splice_dir, tissue + ".junc-counts.txt"))

        id_dict[tissue] = master_sj_file.readline().strip().split('\t')[1:]
        read_depths, id_dict, exc_idxs = check_id_concordance(read_depths, id_dict, tissue)

        for line in master_sj_file:
            line = line.strip().split('\t')
            sj = line[0].replace('chr', '')
            if sj in sj_dict:
                sj_dict[sj][tissue] = [int(x) for x in line[1:] if x not in exc_idxs]

    except IOError:
        print('Unable to open reference file ' + master_sj_file)

    master_sj_file.close()

    return read_depths, id_dict, sj_dict


def calculate_gene_threshold(sj_cov, reads, prop):

    """Returns number of reads required for at least given proportion of junctions to be
    covered by given number of reads"""

    if len(sj_cov) == 1:
        position = 0
    else:
        sj_cov = sorted(sj_cov)
        position = int(math.floor(len(sj_cov) * (1-prop)))
    try:
        req_reads = reads/sj_cov[position]
    except ZeroDivisionError:
        # Return arbitrarily high figure to signpost cases where there is no coverage at crucial position in list of
        # splice junction coverage values (which would otherwise yield an infinite MRSD)
        req_reads = 99999999
    return req_reads


def calculate_list_threshold(values, prop):

    """Returns value in a sorted list of values such that a given proportion of values in the list
    lie below the returned value"""

    if len(values) == 1:
        req_read_depth = values
        position = 0
    else:
        req_read_depth = sorted(values)
        position = int(math.ceil(len(values) * prop)) - 1
    return req_read_depth[position]


def calculate_mrsd(tpt_dict, indiv_sj_dict, sj_dict, id_dict, read_depths, tissue,
                   number_reads, sj_prop, stringency):

    """Takes collated read counts for each splice junction and maps their predicted required sequencing depth
    back to respective transcripts."""

    sample_no = len(read_depths[tissue])
    depths = [read_depths[tissue][x] for x in id_dict[tissue]]

    if tissue not in indiv_sj_dict:
        indiv_sj_dict[tissue] = {}

    for gene in tpt_dict:

        if gene not in indiv_sj_dict[tissue]:
            indiv_sj_dict[tissue][gene] = {}

        tpt_dict[gene][tissue]['mrsd'] = []
        for i in range(sample_no):

            new_gene = []
            for sj in tpt_dict[gene][tissue].keys():

                if sj != 'mrsd' and sj != 'ids':
                    # Divide per-junction read count in given sample by its sequencing depth
                    per_read_cov = sj_dict[sj][tissue][i] / (float(depths[i]))
                    # Calculate per M read SJ coverage
                    new_gene.append(per_read_cov * (10 ** 6))

                    if sj not in indiv_sj_dict[tissue][gene]:
                        indiv_sj_dict[tissue][gene][sj] = []
                    indiv_sj_dict[tissue][gene][sj].append(str(per_read_cov * (10 ** 6)))

            gene_thresh = calculate_gene_threshold(new_gene, number_reads, sj_prop)
            tpt_dict[gene][tissue]['mrsd'].append(gene_thresh)

        tpt_dict[gene][tissue]['mrsd'].append(calculate_list_threshold(tpt_dict[gene][tissue]['mrsd'],
                                                                       stringency))

    return tpt_dict, indiv_sj_dict


def print_mrsds(tpt_dict, ordered_tissue_ids, opt_header, failed_genes):

    """Prints results to text file."""

    with open(args.output_prefix + ".results.mrsd.txt", "w") as opt:
        opt.write("INPUT PARAMETERS\n-------------\n\n")
        opt.write("Selected tissue(s): " + "\t".join(ordered_tissue_ids) + "\n")
        opt.write("Selected read type: " + args.read_type + "\n")
        opt.write("Desired coverage of splice junctions: " + str(args.number_reads) + "\n")
        opt.write("Proportion of splice junctions to reach this coverage per gene: " + str(args.sj_prop) + "\n")
        opt.write("Confidence level: " + str(args.stringency) + "\n")
        opt.write("\n\nGENE-LEVEL REQUIRED READ DEPTH\n---------------------\n\nWhere '-' is given as a " +
                  "required read depth, median coverage of splice junctions for the given gene was 0 reads.\n\nNumber" +
                  " of genes successfully cross-referenced: " + str(len(tpt_dict)) + "\n\n")

        opt.write(opt_header)
        for gene in tpt_dict:

            opt_mrsds = [round(tpt_dict[gene][x]['mrsd'][-1], args.dp) for x in ordered_tissue_ids]
            opt_mrsds = [str(x) if x != 99999999 else "-" for x in opt_mrsds]

            # If only one ID was provided per-line in the input file, add '-' to the id2 column
            if len(tpt_dict[gene][ordered_tissue_ids[0]]['ids']) == 1:
                tpt_dict[gene][ordered_tissue_ids[0]]['ids'].append('-')

            opt.write('\t'.join(tpt_dict[gene][ordered_tissue_ids[0]]['ids'] + opt_mrsds) + '\n')

        if failed_genes:
            opt.write('\n\nThe following genes were not able to be identified - these may be single-exon ' +
                      'genes, this gene name may be deprecated, or may be differently named in this genome ' +
                      'version (GENCODE v19):\n' + '; '.join(failed_genes))


def print_juncs(indiv_sj_dict):

    """ Prints MRSD at individual splice junctions in a separate output file. Suppressed by default as these files
    can be very large. """

    with open(args.output_prefix.rstrip('.') + ".junc-coverage.txt", "w") as opt:

        opt.write('id\tgene\ttissue\tchr\tstart\tend\tjunction_no\tper_M_coverage\n')

        row_no = 0

        for tissue in indiv_sj_dict:
            for gene in indiv_sj_dict[tissue]:
                ordered_sjs = sorted(indiv_sj_dict[tissue][gene].keys(),
                                     key=lambda x: (int(x.split('-')[1]), int(x.split('-')[2])))
                count = 0
                for sj in ordered_sjs:
                    count += 1
                    chrm, start, end = sj.split('-')
                    for cov in indiv_sj_dict[tissue][gene][sj]:
                        row_no += 1
                        opt_list = [str(row_no),
                                    gene,
                                    tissue,
                                    chrm,
                                    start,
                                    end,
                                    str(count),
                                    str(cov)]

                        opt.write('\t'.join(opt_list) + '\n')


def main():

    transcript_model = pickle_load(args.transcript_model)

    read_depths = {}

    # Convert read depths to dictionary
    for tissue in args.tissues:
        print('Processing ' + tissue + '.')
        read_depths = collate_read_depths(tissue, read_depths, READ_TYPES[args.read_type])

    tpt_dict = {}
    sj_dict = {}
    indiv_sj_dict = {}
    failed_genes = []

    # Collate all SJs for transcripts of interest, incorporating both BED and text-style entries
    for file in args.gene_list:
        file_type = parse_input_file_type(file)
        with open(file) as ipt:
            if file_type == "list":
                tpt_dict, sj_dict, failed_genes = collect_list_sjs(tpt_dict, sj_dict, ipt,
                                                                   read_depths, transcript_model, failed_genes)

            elif file_type == "bed":
                tpt_dict, sj_dict = collect_bed_sjs(tpt_dict, sj_dict, ipt, read_depths)

            else:
                sys.exit("Unknown error in input file " + file)

    id_dict = {}

    for tissue in read_depths:

        # Collate junction counts for each sample
        read_depths, id_dict, sj_dict = collate_master_sjs(read_depths, id_dict, sj_dict, tissue)

        # Calculate MRSDs for at both transcript and junction level by cross-referencing sequencing depths and
        # junction counts
        tpt_dict, indiv_sj_dict = calculate_mrsd(tpt_dict, indiv_sj_dict, sj_dict, id_dict, read_depths, tissue,
                                                 args.number_reads, args.sj_prop, args.stringency)

    ordered_tissue_ids = sorted(args.tissues)
    opt_header = '\t'.join(['id1\tid2'] + ['MRSD(' + x + ')' for x in ordered_tissue_ids]) + '\n'

    # Print transcript-level results to file
    print_mrsds(tpt_dict, ordered_tissue_ids, opt_header, failed_genes)

    # Print junction-level results to file
    if args.output_junctions:
        print_juncs(indiv_sj_dict)

    print('Analysis complete.')


if __name__ == "__main__":

    READ_TYPES = {'total': 1, 'unique': 2, 'multimap': 3}
    args = parse_args()
    main()
