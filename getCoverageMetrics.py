#!/usr/bin/env python3

import os, sys

inDir = sys.argv[1]
pattern = sys.argv[2]
wgspattern = sys.argv[3]
def getReads(filestring):

  # Read the data from the input file so we can begin to work with it
  File = open(filestring, 'r')
  Lines = File.readlines()

  # Find specific pieces of data in the lines we collected using different methods
  # Split each line based on tabs; for the line starting with PAIR, take the 6th element
  datalist = []
  reads = 0
  for line in Lines:
    datalist = line.split("\t")
    if datalist[0] == 'PAIR':
      reads = datalist[5]
      # print(reads)
      return reads


def getCoverageFromWGSmetrics(filestring):
  File = open(filestring, 'r')
  Lines = File.readlines()
  # print(Lines)
  lineCount = 0
  for line in Lines:
    datalist = line.split("\t")
    if datalist[0] == 'GENOME_TERRITORY':
      # skip to next line, then take second element
      coverageLine = Lines[lineCount + 1]
      coverage = coverageLine.split("\t")[1]
    lineCount += 1

    
  return coverage

# effectively, need to read in the metrics file, get the number of reads, convert to coverage, then determine the downsample rate (should be a funciton of coverage and desired coverage (30x)

# get the number of reads from the metrics file
# walk through the directory and get all files ending in .alignment_summary_metrics.txt... maybe use os.walk
# with open('samples3.yaml', 'w') as f, open('coverage.txt', 'w') as c:
#   f.write('samples:\n')
with open('coverageValuesFinal.txt', 'w') as c:
  files = [f for f in os.listdir(inDir) if f.endswith(f'{pattern}')]
  print(files)
  # print(files)
  for file in files:
    # if file.startswith('KP') or file.startswith('NPH'):
    #   continue
    filestring = inDir + file
    print(filestring)
    # print(filestring)
    try:
      reads = int(getReads(filestring))
    except:
      print(f'Error in {file}, skipping sample')
      continue
    filestring = filestring.replace(f'{pattern}', f'{wgspattern}')
    if not os.path.exists(filestring):
      print('skipping this sample')
      continue
    print(filestring)
    coverageOld = float(getCoverageFromWGSmetrics(filestring))
    # print(file, reads)
    sample = file.split('/')[-1].split('.')[0]
    # print(reads)
    area_of_coverage = 3e9 / 150 # this is the total area of the genome divided by the read length
    coverage = reads/area_of_coverage # change value as needed
    print(sample, coverage, coverageOld)
    c.write(f'{sample}\t{coverageOld}\t{coverage}\t{reads}\n')
  # c.write(f'{sample}: {coverage}\n')

  # # determine the downsample rate
  # downsampleP = 30/coverage

    # now want to write downsampleP per sample to a file
    # f.write(f'  {sample}: [{str(file)},{str(downsampleP)}]\n')

    