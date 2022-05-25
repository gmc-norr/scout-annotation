#!/usr/bin/env python3

import argparse
import gzip

def mkParser():
  parser = argparse.ArgumentParser(description = "Reformats a VCF file to better fit Scout standards. Removes 0/0 gentotypes and adds AD fields if missing")
  parser.add_argument("--vcf",          type = str,    required = True,   help="a file in vcf format")
  parser.add_argument("--out",          type = str,    required = True,   help="the name of the output file")

  return parser.parse_args()


def reformat(vcf, outname, zipped):
  """Opens vcf file, checks lines and writes new vcf file"""
  if zipped:
    infile = gzip.open(vcf, 'rt')
    out    = gzip.open(outname+'.gz', 'wt')
  else:
    infile = open(vcf, 'r')
    out    = open(outname, 'w')

  line = infile.readline()
  while line.startswith('##'):
    out.write(line)
    line = infile.readline()
  out.write(line)

  for line in infile:
    cols     = line.strip().split()

    # check so genotype is not 0/0 or missing
    if not check_geno(cols):
      continue

    # check if AD in format field, otherwise add it
    line = add_AD(cols)
    
    # check for callers and reformat caller info
    line = add_caller_info(cols)

    # write line
    out.write(line)


def check_geno(cols):
  """Returns False if genotype is 0/0 or missing"""
  geno  = cols[9].split(':')
  return True if not geno[0] in ['0/0', '0|0', './.'] else False


def add_AD(cols):
  """Uses AO and RO info if no DP in format fields and make new line"""
  formats = cols[8].strip().split(':')

  if 'AD' in formats:
    line = '\t'.join(cols)
    return line+'\n'

  elif 'AO' in formats and 'RO' in formats:
    geno = cols[9].strip().split(':')
    AD   = "{},{}".format(geno[formats.index('RO')], geno[formats.index('AO')])
    formats.insert(1,'AD')
    geno.insert(1,AD)
    cols[8] = ":".join(formats)
    cols[9] = ":".join(geno)
    line    = "\t".join(cols)
    return line+'\n'

  else:
    print('No DP, AO or RO info available')
    line = "\t".join(cols)
    return line+'\n'


def add_caller_info(cols):
    """ Find callers with pass variants and add info last in format field """
    formats = cols[7].strip().split(';')
    for i in range(len(formats)):
        if 'CALLERS' in formats[i]:
            caller_info = str(formats[i]).split('=')
            callers = str(caller_info[1])
            if "gatk-haplotype" in callers:
                callers = callers.replace("gatk-haplotype","gatk")
            caller_set = callers.replace(',', '-')
            formats.insert(len(formats), 'set=' + caller_set)
            cols[7] = ";".join(formats)
            line = '\t'.join(cols)
            return line + '\n'


def main():
  args = mkParser()
  print("##  INFO  ###   Writing new vcf file")
  zipped  = True if args.vcf[-3:] == '.gz' else False
  outname = args.out if args.out[-3:] != '.gz' else args.out[:-3]
  reformat(args.vcf, outname, zipped)
  print("##  INFO  ###   Done!")

main()
