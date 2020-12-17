#!/usr/bin/env python

import sys
from seqlib.sam import convert_CIGAR, index_alignment
from seqlib.interval import Interval
import pysam
import mappy as mp
from collections import defaultdict
from copy import deepcopy


def main():
    target = parse_target('./seq_info.txt')
    name = sys.argv[1].split('_')[0]
    bam = pysam.AlignmentFile(sys.argv[1])
    count = defaultdict(int)
    total = 0
    with open(sys.argv[2], 'w') as detail_f, open(sys.argv[3], 'w') as out:
        for read_name, read_info, read_seq in parse_bam(bam):
            total += 1
            if not read_info:
                count['unmap'] += 1
                detail_f.write('{}\tunmap\tunmap\tNone\t{}\n'.format(read_name,
                                                                     read_seq))
                continue
            gDNA_info, target_info = None, None
            for target_name in read_info:
                if target_name.endswith('gDNA'):
                    bond_name = target_name
                    new_read_info = check_gDNA_reads(read_info[target_name],
                                                     read_seq,
                                                     target[target_name]['fa'])
                    gDNA_info = parse_indel(new_read_info, target_name, target,
                                            gDNA=True)
                else:
                    target_info = parse_indel(read_info[target_name],
                                              target_name, target,
                                              bond_name=bond_name)
            result = read_name
            if ((gDNA_info is not None and gDNA_info[1] == 'full') and
                    (target_info is None or gDNA_info[2] < target_info[2])):
                gene = gDNA_info[4].split('_')[0]
                if gDNA_info[2] <= 20 and gDNA_info[3] == 0:
                    result += '\t{}_WT\tWT'.format(gene)
                    count['{}_WT'.format(gene)] += 1
                elif gDNA_info[2] <= 20 and gDNA_info[3] > 0:
                    result += '\t{}_small_indel\tWT_small_indel'.format(gene)
                    count['{}_small_indel'.format(gene)] += 1
                else:
                    result += '\t{}_large_indel\tWT_large_indel'.format(gene)
                    count['{}_large_indel'.format(gene)] += 1
                result += '\t' + '\t'.join(gDNA_info[0])
            elif ((target_info is not None and target_info[1] == 'full') and
                  (gDNA_info is not None and gDNA_info[1] == 'full')):
                gene = target_info[4].split('_')[0]
                if target_info[2] == 0:
                    result += '\t{}_HDR\tHDR'.format(gene)
                    count['{}_HDR'.format(gene)] += 1
                elif target_info[2] <= 20:
                    result += '\t{}_HDR_small_indel'.format(gene)
                    result += '\tHDR_small_indel'
                    count['{}_HDR_small_indel'.format(gene)] += 1
                else:
                    result += '\t{}_HDR_large_indel'.format(gene)
                    result += '\tHDR_large_indel'
                    count['{}_HDR_large_indel'.format(gene)] += 1
                result += '\t' + '\t'.join(target_info[0])
            else:
                indel_lst = []
                target_lst = []
                if target_info:
                    indel_lst.extend(target_info[0])
                    target_lst.append(target_info[4])
                if gDNA_info:
                    indel_lst.extend(gDNA_info[0])
                    target_lst.append(gDNA_info[4])
                indel_type, indel_struc, indel_lst = parse_comindel(indel_lst,
                                                                    read_seq,
                                                                    target_lst,
                                                                    target)
                count[indel_type] += 1
                result += '\t{}\t{}\t{}'.format(indel_type, indel_struc,
                                                '\t'.join(indel_lst))
            result += '\t' + read_seq
            detail_f.write(result + '\n')
        for tag in count:
            percent = count[tag] / total * 100
            out.write('{}\t{}\t{}\t{}\n'.format(name, tag, count[tag],
                                                percent))


def parse_comindel(indel_lst, read_seq, target_lst, target):
    interval = []
    read_len = len(read_seq)
    gene_set = set()
    donor_flag = False
    for indel in indel_lst:
        indel_info = indel.split('|')
        s = int(indel_info[7])
        e = int(indel_info[8])
        seg_type = indel_info[0] + '_' + indel_info[1]
        interval.append([s, e, seg_type, indel])
        gene_name = indel_info[0].split('_')[0]
        target_type = indel_info[0].split('_')[-1]
        if target_type == 'donor':
            donor_flag = True
        gene_set.add(gene_name)
    interval.sort()
    indel_struc = ''
    info = []
    first, last, index = None, None, None
    for n, (s, e, seg, tag) in enumerate(interval):
        if n != 0:
            d = s - index
            if d >= 20:
                seq = read_seq[index: s]
                (align_flag, seg_info,
                 seg_str, left, right) = align_seg(seq, target_lst, target,
                                                   index, read_len - s,
                                                   read_len)
                if align_flag:
                    info.extend(seg_info)
                    indel_struc += '|{}{}|{}'.format(left, seg_str, right)
                else:
                    indel_struc += '|0|{}|0'.format(seq)
            else:
                indel_struc += '|{}'.format(d)
        else:
            first = d = s
            if d >= 20:
                seq = read_seq[:first]
                (align_flag, seg_info,
                 seg_str, left, right) = align_seg(seq, target_lst, target, 0,
                                                   read_len - first, read_len)
                if align_flag:
                    info.extend(seg_info)
                    if left:
                        indel_struc += '{}{}|{}'.format(left, seg_str, right)
                    else:
                        indel_struc += '0|{}|{}'.format(seg_str, right)
                    first -= d - left
                else:
                    indel_struc += '0|{}|0'.format(seq)
            else:
                indel_struc += str(d)
        if last is None or e > last:
            last = e
        index = e
        indel_struc += '|' + seg
        info.append(tag)
    else:
        d = read_len - last
        if d >= 20:
            seq = read_seq[last:]
            align_flag, seg_info, seg_str, left, right = align_seg(seq,
                                                                   target_lst,
                                                                   target,
                                                                   last, 0,
                                                                   read_len)
            if align_flag:
                info.extend(seg_info)
                if right:
                    indel_struc += '|{}{}|{}'.format(left, seg_str, right)
                else:
                    indel_struc += '|{}{}'.format(left, seg_str)
                last += d - right
            else:
                indel_struc += '|0|{}'.format(seq)
        else:
            indel_struc += '|{}'.format(d)
        if read_len - (last - first) <= 40:
            if len(gene_set) == 2:  # two gene
                left_flag, right_flag = False, False
                insert_num = 0
                for n, seg in enumerate(info):
                    gene_target, struc, total_indel, cut_indel = seg.split('|')[:4]
                    total_indel, cut_indel = int(total_indel), int(cut_indel)
                    gene = gene_target.split('_')[0]
                    target_type = gene_target.split('_')[-1]
                    if struc == 'L-LH-C-RH-R':
                        if donor_flag:
                            return '{}_donor_fragment_insert'.format(gene), indel_struc, info
                        else:
                            if total_indel <= 20 and cut_indel == 0:
                                return '{}_WT'.format(gene), 'WT', info
                            elif total_indel <= 20 and cut_indel > 0:
                                return '{}_small_indel'.format(gene), 'WT_small_indel', info
                            else:
                                return '{}_large_indel'.format(gene), 'WT_large_indel', info
                    if n == 0 and target_type == 'gDNA' and struc == 'L-LH-C':
                        left_flag = True
                    if struc == 'LH-I-RH' or struc == 'RH-I-LH':
                        insert_num += 1
                else:
                    if target_type == 'gDNA' and struc == 'C-RH-R':
                        right_flag = True
                if left_flag and right_flag:
                    if insert_num:
                        return '{}_donor_insert({})'.format(gene, insert_num), indel_struc, info
                    elif donor_flag:
                        return '{}_donor_fragment_insert'.format(gene), indel_struc, info
            else:
                insert_num = 0
                left_flag, right_flag = False, False
                out_type = None
                for n, seg in enumerate(info):
                    gene_target, struc, total_indel, cut_indel = seg.split('|')[:4]
                    total_indel, cut_indel = int(total_indel), int(cut_indel)
                    gene = gene_target.split('_')[0]
                    target_type = gene_target.split('_')[-1]
                    if (target_type == 'donor' and total_indel < 10 and
                            (struc == 'LH-I-RH' or struc == 'RH-I-LH')):
                        insert_num += 1
                    elif struc == 'L-LH-C-RH-R':
                        if total_indel <= 20 and cut_indel == 0:
                            out_type = ('{}_WT'.format(gene), 'WT')
                        elif total_indel <= 20 and cut_indel > 0:
                            out_type = ('{}_small_indel'.format(gene), 'WT_small_indel')
                        else:
                            out_type = ('{}_large_indel'.format(gene), 'WT_large_indel')
                    if target_type == 'gDNA':
                        if struc.startswith('L'):
                            left_flag = True
                        if struc.endswith('R'):
                            right_flag = True
                if left_flag and right_flag and insert_num:
                    return '{}_donor_insert({})'.format(gene, insert_num), indel_struc, info
                if out_type:
                    return out_type[0], out_type[1], info
            return 'complicated_structure', indel_struc, info
        else:
            return 'complicated_structure', indel_struc, info


def align_seg(seq, target_lst, target, left_offset, right_offset, read_len):
    seq_len = len(seq)
    flag = False
    read_info = []
    for target_name in target_lst:
        target_fa = target[target_name]['fa']
        if target_name.endswith('gDNA'):
            gDNA = True
        else:
            gDNA = False
        for read in target_fa.map(seq, MD=True):
            left_cigar, right_cigar = '', ''
            if read.q_st > 0 or left_offset > 0:
                left_cigar = '{}S'.format(read.q_st + left_offset)
            if seq_len > read.q_en or right_offset > 0:
                right_cigar = '{}S'.format(seq_len - read.q_en + right_offset)
            if read.strand == 1:  # plus strand
                strand, reverse = '+', False
                cigar = left_cigar + read.cigar_str + right_cigar
            else:  # minus strand
                strand, reverse = '-', True
                cigar = right_cigar + read.cigar_str + left_cigar
            # reverse alignment if possible
            alignment = convert_CIGAR(cigar, read.MD, reverse=reverse)
            aln = ''.join('{}{}'.format(x[0], x[1]) for x in alignment)
            index1, index2 = index_alignment(aln)
            read_info.append([index1, index2, read.r_st, read.r_en, strand,
                              aln, alignment])
            flag = True
        if flag:
            break
    else:
        return False, seq, None, None, None
    read_info.sort()
    info = parse_indel(read_info, target_name, target, bond_name=target_name,
                       gDNA=gDNA)[0]
    seg_structure = ''
    first, last, index = None, None, None
    for n, seg in enumerate(info):
        seg_items = seg.split('|')
        s = int(seg_items[7])
        e = int(seg_items[8])
        seg_info = '_'.join(seg_items[:2])
        if n != 0:
            d = s - index
            seg_structure += '|{}'.format(d)
        else:
            first = s - left_offset
        if last is None or e > last:
            last = e
        index = e
        seg_structure += '|' + seg_info
    return True, info, seg_structure, first, read_len - last - right_offset


def check_gDNA_reads(read_info, read_seq, target_fa):
    read_len = len(read_seq)
    # add padding
    read_lst = deepcopy(read_info)
    for item in read_lst:
        item[1] += 10
    if len(Interval(read_lst)) < len(read_lst):  # realign
        new_read_info = []
        for read in target_fa.map(read_seq, MD=True):
            left_cigar, right_cigar = '', ''
            if read.q_st > 0:
                left_cigar = '{}S'.format(read.q_st)
            if read_len > read.q_en:
                right_cigar = '{}S'.format(read_len - read.q_en)
            if read.strand == 1:  # plus strand
                strand, reverse = '+', False
                cigar = left_cigar + read.cigar_str + right_cigar
            else:  # minus strand
                strand, reverse = '-', True
                cigar = right_cigar + read.cigar_str + left_cigar
            # reverse alignment if possible
            alignment = convert_CIGAR(cigar, read.MD, reverse=reverse)
            aln = ''.join('{}{}'.format(x[0], x[1]) for x in alignment)
            index1, index2 = index_alignment(aln)
            new_read_info.append([index1, index2, read.r_st, read.r_en, strand,
                                  aln, alignment])
        return new_read_info
    else:
        return read_info


def parse_target(fn):
    target = defaultdict(dict)
    with open(fn, 'r') as f:
        for line in f:
            (name, left, donor, right,
                total, cut, in_s, in_e, fa) = line.rstrip().split()
            target_type = 'gDNA' if name.endswith('gDNA') else 'donor'
            left, right = int(left), int(right)
            total, cut = int(total), int(cut)
            in_s, in_e = int(in_s), int(in_e)
            donor = int(donor)
            target[name]['left_bond'] = left
            target[name]['right_bond'] = right
            target[name]['donor'] = donor
            target[name]['total'] = total
            if target_type == 'gDNA':
                target[name]['cut_left'] = cut - 10
                target[name]['cut_right'] = cut + 10
                interval = [[0, left, 'L'], [left, cut - 25, 'LH'],
                            [cut - 25, cut + 25, 'C'],
                            [cut + 25, left + donor, 'RH'],
                            [left + donor, total, 'R']]
                target[name]['fa'] = mp.Aligner(fa, preset='map-pb')
            else:
                interval = [[0, in_s, 'LH'], [in_s, in_e, 'I'], [in_e, total,
                                                                 'RH']]
                target[name]['fa'] = mp.Aligner(fa, preset='sr')
            target[name]['interval'] = interval
            target[name]['rev_interval'] = [[total - x[1], total - x[0], x[2]]
                                            for x in interval][::-1]
    return target


def parse_indel(read_info, target_name, target, bond_name=None, gDNA=False):
    min_indel = None
    indel_info = []
    target_len = target[target_name]['total']
    if not gDNA:  # donor
        left_bond = target[bond_name]['left_bond']
        right_bond = target[bond_name]['right_bond']
        dis = 5
    else:  # gDNA
        left_bond = 0
        right_bond = 0
        cut_left = target[target_name]['cut_left']
        cut_right = target[target_name]['cut_right']
        dis = 20
    for index1, index2, start, end, strand, aln, alignment in read_info:
        if strand == '+':
            interval = target[target_name]['interval']
        else:
            interval = target[target_name]['rev_interval']
            # reverse start and end
            start, end = target_len - end, target_len - start
        aln_info = '{}|{}|{}|{}|{}|{}'.format(strand, start, end,
                                              index1, index2, aln)
        if start >= dis or target_len - end >= dis:
            target_tag = 'partially'
            seg_type = '-'.join(x[2] for x in Interval.mapto([start, end],
                                                             interval))
        else:
            target_tag = 'full'
            if strand == '+':
                seg_type = 'L-LH-C-RH-R' if gDNA else 'LH-I-RH'
            else:
                seg_type = 'R-RH-C-LH-L' if gDNA else 'RH-I-LH'
        cut_indel = 0
        total_indel = 0
        reference_pos = start
        read_tag = 'full'
        for n, (num, t) in enumerate(alignment):
            if t in ('S', 'H'):
                if n == 0:  # first postion
                    if num - left_bond >= dis:
                        read_tag = 'partially'
                else:  # last position
                    if num - right_bond >= dis:
                        read_tag = 'partially'
            if gDNA:
                if t == 'I' and cut_left <= reference_pos <= cut_right:
                    cut_indel += num
                if t in ('D', 'U'):
                    if (cut_left <= reference_pos + num and
                            reference_pos <= cut_right):
                        cut_indel += num
                if t in ('M', 'U', 'D'):
                    reference_pos += num
            else:
                if t in ('D', 'U', 'I'):
                    cut_indel += num
            if t in ('D', 'U', 'I'):
                total_indel += num
        if target_tag == 'partially' or read_tag == 'partially':
            tag = 'partially'
        else:
            tag = 'full'
        if min_indel is None or total_indel < min_indel:
            min_indel = total_indel
            min_cut_indel = cut_indel
            final_tag = tag
        indel_info.append('{}|{}|{}|{}|{}'.format(target_name, seg_type,
                                                  total_indel, cut_indel,
                                                  aln_info))
    return [indel_info, final_tag, min_indel, min_cut_indel, target_name]


def parse_bam(bam):
    read_name = None
    read_info = defaultdict(list)
    seq_flag = False
    read_seq = None
    for read in bam:
        if read.query_name != read_name and read_name is not None:
            yield read_name, read_info, read_seq
            read_info = defaultdict(list)
            seq_flag = False
            read_seq = None
        read_name = read.query_name
        if read.is_unmapped:
            read_seq = read.get_forward_sequence()
            continue
        strand = '-' if read.is_reverse else '+'
        reverse = True if read.is_reverse else False
        # reverse alignment if possible
        alignment = convert_CIGAR(read.cigarstring,
                                  read.get_tag('MD'), reverse=reverse)
        aln = ''.join('{}{}'.format(x[0], x[1]) for x in alignment)
        index1, index2 = index_alignment(aln)
        if not seq_flag and 'H' not in read.cigarstring:
            seq_flag = True
            read_seq = read.get_forward_sequence()
        # not reverse start and end
        info = [index1, index2, read.reference_start, read.reference_end,
                strand, aln, alignment]
        read_info[read.reference_name].append(info)
    else:
        yield read_name, read_info, read_seq


if __name__ == '__main__':
    main()
