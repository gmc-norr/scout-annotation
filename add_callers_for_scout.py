#!/usr/bin/env python3.6

import argparse
import gzip

def mkParser():
  parser = argparse.ArgumentParser(description = "Reformats a VCF file to better fit Scout standards. Adds info about caller pass variants ")
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

    # check for callers and reformat caller info
    line = add_caller_info(cols)

    # write line
    out.write(line)

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