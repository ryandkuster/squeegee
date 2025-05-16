#!/usr/bin/env python3

"""
Open the vg deconstruct vcf and parse/summarize the variants

SVTYPE can be SNP, MNP, INS, or DEL (currently no support for BND/DUP/INV)
SVLEN is the length of the alt minus ref (SNP and MNP are REF len)
END is calculated as follows:
    INS: END = POS + (length REF - 1)
    DEL: END = POS + abs(length ALT - length REF)
    SNP: END = POS + (length REF - 1) (i.e., it's just POS)
    MNP: END = POS + (length REF - 1)

The traversal path starts and ends with a shared node. The first node is
included in the actual reference and alternate sequences if all begin
with the same base.

Example:
In this variant, the C and A are stored as individual nodes 100930 and
100932. The T and G are stored in 100929 and 100931.
['CA', 'TG']
['>100928>100930>100932>100933',
 '>100928>100929>100931>100933']

In this variant, the GG is stored as a single node 100964. The CA is
stored in 100963.
['GG', 'CA']
['>100962>100964>100965',
 '>100962>100963>100965']

In this variant, I believe the T is likely stored in 100141, and at
times when there is a single base as an alt this seems to be how it is
represented.
['TAATA', 'T']
['>100141>100142>100143',
 '>100141>100143']

In this variant, the single G is likely stored in 99992, as AA does not
share this common base.
['G', 'GA', 'AA']
['>99991>99992>99995',
 '>99991>99992>99994>99995',
 '>99991>99993>99994>99995']

I belive the first shared A is in 100009, C is in 100013, and A is in
100014. Nodes 100015 and 100016 are included to cover all nodes, much
like how a gvcf covers all loci.
['AC', 'A', 'AA']
['>100009>100013>100014>100015>100016',
 '>100009>100014>100015>100016',
 '>100009>100010>100014>100015>100016']
"""

import argparse
import gzip
import os
import subprocess
import sys


def parse_user_input():
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-v', type=str, required=True,
            help='Path to bcftools norm biallelic vg deconstruct file.')

    parser.add_argument('-o', type=str, required=True,
            help='Path to output vcf.')

    parser.add_argument('-t', type=str, nargs='+', required=False,
            default = ["SNP", "MNP", "INS", "DEL"],
            help='Space separated list of variant types to keep.')

    parser.add_argument('-m', type=int, required=False,
            default = 0,
            help='Minimum SVLEN to keep.')

    parser.add_argument('-M', type=int, required=False,
            default = 1_000_000_0000,
            help='Minimum SVLEN to keep.')

    parser.add_argument('-l', type=int, required=False,
            default = 1_000_000,
            help='Maximum level in snarl tree (LV) to keep.')

    args = parser.parse_args()

    return args


def process_var(var_ls, sv_min, sv_max):
    """
    Split variants by comma (if no comma list length is 1).
    Check variant list for longest variant. If longest variant is over a
    between sv_min and sv_max in either ref or alt, then return True.
    """
    var_min = len(min(var_ls, key=len))
    var_max = len(max(var_ls, key=len))
    var_diff = var_max - var_min

    if sv_min <= var_diff <= sv_max:
        return True

    return False


def last_common_node(var_ls, traversal):
    """
    Search for nodes found in all allele traversal paths. This needs to
    look from the end of the list, as first node will at times be a
    common point in the graph, regardless of if the ref/alt begin there.
    Use the shortest allele path.
    Return index of common final node (default -1 is last node).
    """
    trav_ls = [i[1:].split(">") for i in traversal]
    short_path = len(trav_ls[trav_ls.index(min(trav_ls, key=len))])

    for i in range(-short_path, -1):
        node_ls = []
        for j in trav_ls:
            node_ls += j[i:]
        if len(set(node_ls)) == -i:
            return i

    return -1


def get_variant_type(var_ls):
    """
    Using the ref allele (var_ls[0]), define the alt alleles as
    insertion/deletion and its size.
    """
    ref_len = len(var_ls[0])

    for alt in var_ls[1:]:
        var_diff = len(alt) - ref_len

        if var_diff == 0 and ref_len == 1:
            var_type = "SNP"
        if var_diff == 0 and ref_len != 1:
            var_type = "MNP"
        elif var_diff > 0:
            var_type = "INS"
        elif var_diff < 0:
            var_type = "DEL"

    return var_type


def process_variant(var, var_type):
    """
    Check the var_type and create an END position relative to var[1].
    INS END will be identical to var[1] and DEL END will be var[1] plus
    the absolute value of var_diff.

    END is calculated as follows:
        SNP: END = POS + (length REF - 1) (i.e., it's just POS)
        MNP: END = POS + (length REF - 1)
        INS: END = POS + (length REF - 1)
        DEL: END = POS + (length REF - 1)
    """
    ref_len = len(var[3])
    alt_len = len(var[4])
    var_diff = alt_len - ref_len
    new_var = var[:]
    end = int(var[1]) + ref_len - 1

    if var_type == "SNP":
        new_var[7] = f"SVTYPE={var_type};SVLEN={ref_len};END={end};{var[7]}"
    if var_type == "MNP":
        new_var[7] = f"SVTYPE={var_type};SVLEN={ref_len};END={end};{var[7]}"
    if var_type == "INS":
        new_var[7] = f"SVTYPE={var_type};SVLEN={var_diff};END={end};{var[7]}"
    if var_type == "DEL":
        new_var[7] = f"SVTYPE={var_type};SVLEN={var_diff};END={end};{var[7]}"

    return new_var


def main():
    commit = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode().strip()
    print(commit)
    sys.exit()
    args = parse_user_input()
    assert os.path.abspath(args.v) != os.path.abspath(args.o)
    if not args.o.endswith('.gz'):
        raise argparse.ArgumentTypeError("Output filename must end with '.gz'")

    write_headers = True
    norm_pass = False
    print(f"##{' '.join(sys.argv)}\n")

    with gzip.open(args.v, 'rt') as f, gzip.open(args.o, 'wt') as o:
        for line in f:
            if line.startswith("##contig") and write_headers:
                o.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
                o.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
                o.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of variant">\n')
                write_headers = False

            if line.startswith("##"):
                o.write(line)
                if "bcftools_normCommand" in line:
                    norm_pass = True
            elif line.startswith("#C"):
                if norm_pass is False:
                    print("No bcftools norm command in header!")
                    print("\nExample:\n")
                    print("bcftools norm --output out.vcf --output-type v -m-any --fasta-ref ref.fna --check-ref w in.vcf\n")
                    print("Squeegee relies on bcftools normalized, bi-allelic sites from bcftools norm.")
                    choice = input("Continue anyway? (This is likely a bad choice) : ")
                    if choice in ["y", "Y", "yes", "Yes"]:
                        pass
                    else:
                        sys.exit("See README at https://github.com/ryandkuster/squeegee for help")

                o.write(f"##{' '.join(sys.argv)}\n")
                o.write(line)
            else:
                var = line.rstrip().split()
                ref = var[3].split(",")
                alt = var[4].split(",")
                var_ls = ref + alt
                assert len(var_ls) == 2, "Squeegee relies on biallelic variants only."

                """
                Check if the variant contains any allele between the
                defined sv_min and sv_max.
                """
                if process_var(var_ls, args.m, args.M):
                    level = int(var[7].split(";")[5][3:].split(",")[0])
                else:
                    continue

                if level > args.l:
                    continue

                #traversal = var[7].split(";")[3][3:].split(",")
                #last_node = last_common_node(var_ls, traversal)

                var_type = get_variant_type(var_ls)
                if var_type not in args.t:
                    continue

                """
                Now we process the passing variant.
                """
                new_var = process_variant(var, var_type)
                o.write("\t".join(new_var) + "\n")


if __name__ == "__main__":
    main()



