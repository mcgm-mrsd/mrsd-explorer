# Takes set of Cummings outputs and merges

import sys

# argument 1 = file containing list of Cummings outputs to merge, one per line, argument 2 = prefix for output files,
# which must uniquely identify this dataset from any others

with open(sys.argv[1]) as file_list:
    files = [x.strip() for x in file_list.readlines()]

sj_dict = {}
ids = []
chrms = []

for file in files:

    print('Opening ' + file)
    line_count = 0

    with open(file) as sj_ipt:
        for line in sj_ipt:
            
            if line.startswith('Gene'):
                cummings_header = line[:]
                continue

            line_count += 1
#            if line_count % 10000 == 0:
#                print(str(line_count) + ' lines read.')

            line = line.strip().split('\t')

            if '*' in line[3] or '*' in line[4]:
                continue

            if line[2] not in chrms:
                chrms.append(line[2])

            line[5] = int(line[5])
            line[-1] = line[-1].split(',')

            # Built up list of IDs for construction of MRSD output file
            for item in line[-1]:
                    if item.split(':')[0] not in ids:
                        ids.append(item.split(':')[0])

            sj_id = '//'.join((line[2], line[3], line[4]))

            if sj_id not in sj_dict:
                sj_dict[sj_id] = [line[:]]
            else:

                # Check for accidental addition of duplicate entries - will not catch cases where different read counts
                # are supported for the same junction in the same individual in multiple files, however this should be
                # avoidable with correct upstream processing

                duplicates = []
                for gene in sj_dict[sj_id]:
                    for item in line[-1]: 
                        if item in gene[-1]:
                            duplicates.append(item)
                
                line[-1] = [x for x in line[-1] if x not in duplicates]

                if line[0] in [x[0] for x in sj_dict[sj_id]]:
                    idx = [x[0] for x in sj_dict[sj_id]].index(line[0])
                    sj_dict[sj_id][idx][5] = sj_dict[sj_id][idx][5] + sum([int(x.split(':')[1]) for x in line[-1]])
                    sj_dict[sj_id][idx][-1].extend(line[-1])
                else:
                    sj_dict[sj_id].append(line[:])

ids = sorted(ids)
print('Dataset contains the following sample IDs: ' + '\n'.join(ids))

with open(sys.argv[2].rstrip('.') + '.csf.splicing.txt', 'w') as cummings_out:
    with open(sys.argv[2].rstrip('.') + '.junc-counts.txt', 'w') as mrsd_out:

        complete_sjs = []

        cummings_out.write(cummings_header)
        mrsd_out.write('SJID\t' + '\t'.join(ids) + '\n')

        print('Sorting coordinates.')
        is_sorted = False

        key_count = 0

        for sj in sorted(sj_dict.keys(), key=lambda x: (chrms.index(x.split('//')[0]), int(x.split('//')[1]),
                                                        int(x.split('//')[2]))):
            
            if not is_sorted:
                is_sorted = True
                print('Sorting complete.')

            key_count += 1
#            if key_count % 10000 == 0:
#                print(str(key_count) + ' entries processed.')
#                print(len(out_data))

            for entry in sj_dict[sj]:

                out_data = entry[:]
                out_data[-1] = ','.join(sorted(out_data[-1], key=lambda x: int(x.split(':')[1])))
                cummings_out.write('\t'.join([str(x) for x in out_data]) + '\n')

            # For MRSD file, only print splice junctions once regardless of the number of transcripts they are present
            # in, to avoid complications with downstream analysis

            non_redund_sj = '{}-{}-{}'.format(sj.split('//')[0],
                                              sj.split('//')[1],
                                              sj.split('//')[2])

            junc_counts = ['0'] * len(ids)
            for item in sj_dict[sj][0][-1]:
                samp_id, count = item.split(':')
                junc_counts[ids.index(samp_id)] = count
            mrsd_out.write(non_redund_sj + '\t' + '\t'.join([str(x) for x in junc_counts]) + '\n')
            
            del sj_dict[sj]
