import argparse
import sys


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument("transcripts", action="store",
                        help="Path to file containing transcripts of interest, one per line")
    parser.add_argument("annotation_gtf", action="store",
                        help="Path to annotation GTF")
    parser.add_argument("output_prefix", action="store",
                        help="Prefixes for output files")
    parser.add_argument("--enforce_transcript_version", action="store_true", default=False,
                        help="Require transcript version to match in given annotation if a version number is given")
    return parser.parse_args()


def validate_paths(gtf_path, transcript_path):
    errs = False

    try:
        gtf = open(gtf_path)
    except IOError:
        print("Unable to open reference GTF. Ensure path is correct.")
        errs = True

    try:
        transcripts = open(transcript_path)
    except IOError:
        print("Unable to open transcript list. Ensure path is correct")
        errs = True

    if errs:
        sys.exit()
    else:
        return gtf, transcripts


def fetch_field(details, delim, field):

    for item in details.split(delim):
        if field in item:
            return item.split('"')[-2]
    return False


def process_gtf(exon_dict, gtf):

    matched_transcripts = set()
    for line in gtf:
        if line.startswith('#'):
            continue

        line = line.strip().split('\t')

        if line[2] == 'exon':
            tid = fetch_field(line[-1], ';', 'transcript_id')

            if tid in exon_dict or tid.split('.')[0] in exon_dict:

                if tid.split('.')[0] in exon_dict:

                    new_tid = tid.split('.')[0]
                else:
                    new_tid = tid

                if not exon_dict[new_tid]:
                    matched_transcripts.add(new_tid)

                    # GTF
                    gene_name = fetch_field(line[-1], ';', 'gene_name')

                    # GFF
                    if not gene_name:
                        gene_name = fetch_field(line[-1], ';', 'gene ')

                    exon_dict[new_tid] = [tid, gene_name, line[0], line[6], [[line[3], line[4]]]]
                else:
                    exon_dict[new_tid][-1].append([line[3], line[4]])

    gtf.close()

    return matched_transcripts, exon_dict


def print_introns(transcript, opt):

    strand = transcript[3]

    if len(transcript[-1]) == 1:
        return

    for i in range(len(transcript[-1])-1):
        if strand == '+':
            opt.write('\t'.join((transcript[2], transcript[-1][i][1],
                                 transcript[-1][i+1][0], '_'.join((transcript[0], transcript[1])))) + '\n')
        else:
            opt.write('\t'.join((transcript[2], transcript[-1][i+1][1],
                                 transcript[-1][i][0], '_'.join((transcript[0], transcript[1])))) + '\n')


def print_unmatched(unmatched, unmatched_out):
    with open(unmatched_out, 'w') as opt:
        for item in unmatched:
            opt.write(item + '\n')


def main():
    gtf, transcripts = validate_paths(args.annotation_gtf, args.transcripts)

    transcript_list = [x.strip().split('.')[0] if not args.enforce_transcript_version
                       else x.strip() for x in transcripts.readlines()]

    exon_dict = {x: [] for x in transcript_list}
    matched_transcripts, exon_dict = process_gtf(exon_dict, gtf)

    with open(args.output_prefix.rstrip('.') + '.introns.bed', 'w') as opt:
        opt.write('#CHR\tSTART\tEND\tID\n')
        for transcript in exon_dict:
            if exon_dict[transcript]:
                print_introns(exon_dict[transcript], opt)

    print_unmatched([x for x in transcript_list if x not in matched_transcripts],
                    args.output_prefix + '.unmatched_ids.txt')


if __name__ == "__main__":
    args = parse_args()
    main()
