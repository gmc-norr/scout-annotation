import argparse
import gzip


def mkParser():
    parser = argparse.ArgumentParser(description = "Removes genes from VCF file not present in a defined gene panel")
    parser.add_argument("--vcf",          type = str,    required = True,   help="a file in vcf format")
    parser.add_argument("--panel",        type = str,    required = False,  help="a genepanel file with all genes included listed")
    parser.add_argument("--out",          type = str,    required = True,   help="the name of the output file")

    return parser.parse_args()


def read_panelfile(panel):
    """
    Make a set with all the HGNC ids in the gene panel
    """
    geneset = set()
    with open(panel, 'r') as fh:
        for line in fh:
            if not line.startswith('#'):
                cols = line.strip().split()
                gene = cols[0]
                geneset.add(gene)
    return geneset


def clean_vcf(vcf, zipped, geneset, out):
    """
    Loop over the file and write new one
    """
    if zipped:
        out = out if out[-3:] == '.gz' else out+'.gz'
        fh  = gzip.open(vcf, 'rt')
        out = gzip.open(out, 'wt')
    else:
        fh  = open(vcf, 'r')
        out = open(out, 'w')

    line = fh.readline()
    while line.startswith('##'):
        if line.startswith('##INFO=<ID=CSQ'):
            f_col      = line.split(':')[-1].strip()
            format_col = f_col.strip('\">').split('|')
        out.write(line)
        line = fh.readline()
    out.write(line)

    hgnc = find_ann(format_col, 'HGNC_ID')                  # the index pos of where the ghnc id is listed in the csq annotation

    for line in fh:
        clean = parse(line, geneset, hgnc)
        if clean:
            out.write(clean+'\n')

    fh.close()
    out.close()


def find_ann(cols, annotation):
    """
    Locate index of annotation in a given list
    """
    for ann in cols:
        if ann.startswith(annotation):
            return cols.index(ann)
    print('### ERROR ###    No annotation found, exiting!...')
    exit()


def parse(line, geneset, hgnc):
    """
    Split up line into smaller annotation parts and remove genes not in geneset
    """
    linecol           = line.strip().split()
    ann_col           = linecol[7].split(';')
    csqix             = find_ann(ann_col, 'CSQ')           # find which column contains the csq annotation
    transcripts       = ann_col[csqix].split(',')          # split vep annotation into transcripts
    clean_transcripts = []
    for transcript in transcripts:
        t_cols = transcript.strip().split('|')
        if t_cols[hgnc] in geneset:                        # remove all transcripts not found in geneset
            clean_transcripts.append(transcript)

    if len(clean_transcripts) == 0:                        # only build new line if there are any transcripts left
        return None

    # build new line
    new_transcripts     = ','.join(clean_transcripts)
    if not new_transcripts.startswith('CSQ'):
        new_transcripts = 'CSQ='+new_transcripts
    new_ann_col         = ';'.join(ann_col[:csqix]+[new_transcripts]+ann_col[csqix+1:])
    new_linecol         = '\t'.join(linecol[:7]+[new_ann_col]+linecol[8:])

    return new_linecol


def main():
    args = mkParser()
    print("##  INFO  ###   Running")
    print("##  INFO  ###   Writing new vcf file")
    geneset   = read_panelfile(args.panel)
    zipped    = True if args.vcf[-3:] == '.gz' else False
    clean_vcf(args.vcf, zipped, geneset, args.out)
    print("##  INFO  ###   Done!")


if __name__ == "__main__":
    main()
