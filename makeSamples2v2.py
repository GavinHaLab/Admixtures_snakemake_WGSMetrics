#!/usr/bin/python

# Import packages as needed
import sys
import getopt

#---------------------------
#--- GET INPUT ARGUMENTS ---
#---------------------------

# Get options passed in as arguments when the script is called
# Once we get file names, we can work with their contents, etc.
# Set defaults in case of an error, etc.
# Note that we already know the output file name and location, but 
# it is provided as an argument for the sake of snakemake's I/O rules
infile = ''
outfile = ''

print('\n')
print('OUTPUT PRINTING FOR makeSamples2.py, LOOK IN logs/samples2yaml.txt FOR OUTPUT')

# Define the main() function that will handle the arguments passed in
def main(argv):

  # Declare global variables to which we will assign values
  global infile
  global outfile

  # Handle trying and exceptions
  try:
    opts, args = getopt.getopt(argv,"hi:o:",["infile=", "outfile="])
  except getopt.GetoptError:
    print('makeSamples2.py -i <infile> -o <outfile>')
    sys.exit(2)

  # Handle argument options
  for opt, arg in opts:
    if opt == '-h':
      print('makeSamples2.py -i <infile> -o <outfile>')
      sys.exit(2)
    elif opt in ("-i", "--infile"):
      infile = arg
    elif opt in ("-o", "--outfile"):
      outfile = arg

  # Print the input arguments to the screen for confirmation
  print('Input file is...  ', infile)
  print('Output file is... ', outfile) 

# Include this to activate the main() function and assign values to variables
if __name__ == "__main__":
	main(sys.argv[1:])

#--------------------------------
#--- PARSE SAMPLES1.YAML FILE ---
#--------------------------------
  
# Read the data from the input file so we can begin to work with it
File = open(infile, 'r')
Lines = File.readlines()

# Start reading out various important values from samples1.yaml file
print('-----------------------')
print('Important values for the process...')

# First, figure out which line numbers have the titles we need
# The data we want will be in between them
# ASSUMPTION: The order of these YAML sections needs to be maintained
# ASSUMPTION: The list order within each YAML section needs to be maintained
print('Find section headings in samples1.yaml:')
puritiesLine = 0
mixturesLine = 0
coverageLine = 0
fractionLine = 0
dataval = ''
count = 0
for line in Lines:
  count += 1
  dataval = line[0:9]
  if dataval == 'purities:':
    print('Purities on line ', count)
    puritiesLine = count
  if dataval == 'mixtures:':
    print('Mixtures on line ', count)
    mixturesLine = count
  if dataval == 'coverage:':
    print('Coverages on line ', count)
    coverageLine = count
  dataval = line[0:10]
  if dataval == 'fractions:':
    print('Fractions on line ', count)
    fractionLine = count

# Record the number of lines in the YAML file
linecount = count
print('Lines in samples1.yaml: ', linecount)

# Now that we know where the purity data is, let's extract sample ID and purity
i = 0
substr = ':'
sampleIDArray = []
purityArray = []
tempArray = []
for i in range(puritiesLine, mixturesLine-1):
  if substr in Lines[i] and Lines[i][0:1] != '#':
    tempArray = Lines[i].split(substr)
    sampleIDArray.append(tempArray[0].strip())
    purityArray.append(tempArray[1].strip())
numSamples = int( len(sampleIDArray) )
print('Number of samples: ', numSamples)

# Now that we know where the coverage data is, let's get mixture ID and coverages
mixIDArray = []
coverArray = []
for i in range(coverageLine, fractionLine - 1):
  if substr in Lines[i] and Lines[i][0:1] != '#':
    tempArray = Lines[i].split(substr)
    mixIDArray.append(tempArray[0].strip())
    coverArray.append(tempArray[1].strip())
numMixtures = int( len(mixIDArray) )
print('Number of mixtures: ', numMixtures)

# We need to get the TF fraction data as well: a TF array for each mixture
# startTF = []
# endTF = []
# incTF = []
for i in range(fractionLine, linecount):
  if substr in Lines[i] and Lines[i][0:1] != '#':
    tempArray = Lines[i].split(',')
    # need to change to allow for input list of TFs instead of range
    tumorFractionList = tempArray
    # need to ensure first element has [ stripped and last element has ] stripped
    tumorFractionList[0] = tumorFractionList[0].split('[')[1]
    tumorFractionList[-1] = tumorFractionList[-1].split(']')[0]

    # startTF.append( tempArray[0].split('[')[1] )
    # endTF.append( tempArray[1].strip() )
    # incTF.append( (tempArray[2].split(']')[0]).strip() )

# Finally, we need to get the mixing partners or pairs of sample IDs
# ASSUMPTION: For each mixture listed, the sample order must be [tumor, normal]
tumorIDArray = []
normalIDArray = []
for i in range(mixturesLine, coverageLine - 1):
  if substr in Lines[i] and Lines[i][0:1] != '#':
    tempArray = Lines[i].split('"')
    tumorIDArray.append(tempArray[1])
    normalIDArray.append(tempArray[3])

#------------------------------------
#--- CREATE LIST OF METRICS FILES ---
#------------------------------------

# First we need the prefix and suffix for each sample metrics file
prefix = './results/sources/'
# suffix = '-chrsOnly.dups_removed.alignment_summary_metrics.txt'
suffix = '-chrsOnly.dups_removed.WGS-metrics.txt'
metricFileArray = []
for i in range(0, numSamples):
  metricFileArray.append( prefix + sampleIDArray[i] + suffix )

#------------------------------------
#--- GET READS FROM METRICS FILES ---
#------------------------------------

# Define a function that will get reads from a location in a file
def getCoverage(filestring):
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

# Now use the function we defined to get reads from metrics files we've already created
sampleCovs = []
for i in range(0, numSamples):
  sampleCovs.append( float(getCoverage(metricFileArray[i])) )
  #print 'Sample # ' + str(i+1) + ', reads = ', SampleCovs[i]

# Now we have everything we need for calculating downsampling probabilities

#----------------------------
#--- START MAIN FOR LOOPS ---
#----------------------------

# The structure is something like this
# For each mixture...
#  Use coverage and partners and purity
#  Find TF start, end, increment for this mixture
#  For each TF in that range...
#    Determine T/N probabilities
#    Write corresponding lines to samples2.yaml

# Create empty arrays to store the new YAML listings we will create below 
newSampleList = []
newMixtureList = []

# Start iterating over the mixture pairings we have
for i in range(0, numMixtures):

  # Gather up the various data we need for this mixture run
  thisCover = float( coverArray[i] )
  thisTumorID = tumorIDArray[i]
  thisNormalID = normalIDArray[i]
  thisMixID = mixIDArray[i]
  print('-----------------------')
  print('Mixture ' + str(i+1) + ': mixID = ' + thisMixID + ', tumor = ' + thisTumorID + ', normal = ' + thisNormalID)  
  print('Coverage = ' + str(thisCover) + 'X')

  # We need to lookup the purity of the tumor partner and reads in that partner
  thisTumorCov = 0
  thisPurity = 0
  for t in range(0,numSamples):
    if sampleIDArray[t] == thisTumorID:
      thisPurity = float( purityArray[t] )
      thisTumorCov = float( sampleCovs[t] )
  print('Tumor sample purity = ' + str(thisPurity))
  print('Tumor sample total reads = ' + str(thisTumorCov))

  # We need to lookup the reads in the normal partner
  thisNormalCov = 0
  for n in range(0,numSamples):
    if sampleIDArray[n] == thisNormalID:
      thisNormalCov = float( sampleCovs[n] )
  print('Normal sample total reads = ' + str(thisNormalCov))

  # Calculate tumor and normal reads for the tumor sample
  tCoverTumor = thisPurity * thisTumorCov
  # tCoverNormal = thisCover - tCoverTumor # need to replace this...
  tCoverNormal = 0 ## replacing with 0 since pdx...
  print()
  # print('Tumor sample has ' + str(tReadsTumor) + ' tumor reads and ' + str(tReadsNormal) + ' normal reads')

  # Start iterating over the TF values
  # for j in range(thisStartTF, thisEndTF+1, thisIncTF):
  for j in tumorFractionList:

    #----------------------------------
    #--- CALCULATE DOWNSAMPLE PROBS ---
    #----------------------------------

    thisTF = round( float(j) / float(100), 2)
    thisTFLabel = str(thisTF)
    digits = len( str(thisTF) )
    if digits == 3:
      thisTFLabel = str(thisTF) + '0'    
    print('Desired TF = ' + thisTFLabel)  

    #We need to solve the following system of equations for tProb and nProb.
    # thisTumorCov * tProb + thisNormalCov * nProb = thisCover
    # (thisTumorCov * thisPurity * tProb) / (thisTumorCov * tProb + thisNormalCov * nProb) = thisTF
    #
    #Solution:
    tProb = (thisTF * thisCover) / (thisTumorCov * thisPurity)
    nProb = (thisCover - thisTumorCov * tProb) / thisNormalCov
   
    #print
    print(f'{thisTumorCov} * tProb + {thisNormalCov} * nProb = {thisCover}')
    print(f'({thisTumorCov} * {thisPurity} * tProb) / ({thisTumorCov} * tProb + {thisNormalCov} * nProb) = {thisTF}')
    print("Solution:")
    print("tProb:", tProb)
    print("nProb:", nProb)

    if tProb < 0 or tProb > 1 or nProb < 0 or nProb > 1:
      print("No possible solution for this mixture.")
      sys.exit()
   
    #----------------------------------
    #--- CREATE TEXT LINES FOR YAML ---
    #----------------------------------
    # Now that we have everything, we can assemble a couple arrays
    # First we need to make sure the probabilities have the right number of digits, i.e. 0.XXXX
    digits = len( str(tProb) )
    tProbLabel = str(tProb) 
    if digits == 3: tProbLabel = str(tProb) + '000'
    elif digits == 4: tProbLabel = str(tProb) + '00'
    elif digits == 5: tProbLabel = str(tProb) + '0'

    digits = len( str(nProb) )
    nProbLabel = str(nProb) 
    if digits == 3: nProbLabel = str(nProb) + '000'
    elif digits == 4: nProbLabel = str(nProb) + '00'
    elif digits == 5: nProbLabel = str(nProb) + '0'

    # Now we can assemble a line of test for the samples
    # The pattern is   sampleID_p0.XXXX: [0.XXXX, results/sources/ID-chrsOnly.dups_removed.bam]
    newTumorLabel = thisTumorID + '_p' + tProbLabel
    newNormalLabel = thisNormalID + '_p' + nProbLabel
    addTumorString = '  ' + newTumorLabel + ': [' + tProbLabel + ', results/sources/' + thisTumorID  + '-chrsOnly.dups_removed.cram]'
    addNormalString = '  ' + newNormalLabel + ': [' + nProbLabel + ', results/sources/' + thisNormalID  + '-chrsOnly.dups_removed.cram]'

    # Now add these lines to our running list, first tumor then normal
    newSampleList.append(addTumorString)
    newSampleList.append(addNormalString)   

    # Now we need to make the line for the new mixture format, tumor then normal
    # The pattern is mixID_TF0.XX: ["newTumorLabel", "newNormalLabel"]
    newMixLabel = '  ' + thisMixID + '_TF' + thisTFLabel + ': ["' + newTumorLabel + '", "' + newNormalLabel + '"]'
    newMixtureList.append(newMixLabel)


    # # Back-convert the current TF value to a decimal with two digits, i.e. 0.XX
    # thisTF = round( float(j) / float(100), 2)
    # thisTFLabel = str(thisTF)
    # digits = len( str(thisTF) )
    # if digits == 3:
    #   thisTFLabel = str(thisTF) + '0'    
    # print('   TF = ' + thisTFLabel)

    # # Now we need to find the downsample probabilities, if they exist, for our pair
    # # First, get the overall numbers we need to satisfy for the target
    # targetReads = int( round( thisCover * int(3e7) ) )
    # targetReadsTumor = int( round( thisTF * targetReads ) )
    # targetReadsNormal = targetReads - targetReadsTumor
    # print('   Target reads: total = ' + str(targetReads) + ', tumor = ' + str(targetReadsTumor) + ', normal = ' + str(targetReadsNormal))

    # # ASSUMPTION: Note that we are truncating probabilities to 0.XXXX, which might not be precise enough
    # # Calculate tumor reads, check if we have enough for what we need
    # tProb = 0
    # if tReadsTumor < targetReadsTumor:
    #   print('   Insufficient tumor reads in tumor sample. Prob = 0')
    # else:
    #   tProb = round( float(targetReadsTumor) / float(tReadsTumor), 4)
    #   print('   Tumor downsample probability = ' + str(tProb))

    # # Now we need to find how many reads we need from the normal sample, if any
    # # This is like the total we need minus normal reads in the downsampled tumor file
    # needNormal = targetReadsNormal - int( round( tProb * tReadsNormal ) )
     
    # # Let's see if we have enough reads in the normal sample; if so, find the probability
    # nProb = 0
    # if needNormal > 0:
    #   if thisNormalCov < needNormal:
    #     print('   Insufficient reads in the normal sample. Prob = 0')
    #   else:
    #     nProb = round( float(needNormal) / float(thisNormalCov), 4)
    #     print('   Normal downsample probability = ' + str(nProb))
    # else:
    #   print('   Not possible, too many normal reads in tumor sample...')

    # #----------------------------------
    # #--- CREATE TEXT LINES FOR YAML ---
    # #----------------------------------

    # # Now that we have everything, we can assemble a couple arrays
    # # First we need to make sure the probabilities have the right number of digits, i.e. 0.XXXX
    # digits = len( str(tProb) )
    # tProbLabel = str(tProb) 
    # if digits == 3: tProbLabel = str(tProb) + '000'
    # elif digits == 4: tProbLabel = str(tProb) + '00'
    # elif digits == 5: tProbLabel = str(tProb) + '0'

    # digits = len( str(nProb) )
    # nProbLabel = str(nProb) 
    # if digits == 3: nProbLabel = str(nProb) + '000'
    # elif digits == 4: nProbLabel = str(nProb) + '00'
    # elif digits == 5: nProbLabel = str(nProb) + '0'

    # # Now we can assemble a line of test for the samples
    # # The pattern is   sampleID_p0.XXXX: [0.XXXX, results/sources/ID-chrsOnly.dups_removed.bam]
    # newTumorLabel = thisTumorID + '_p' + tProbLabel
    # newNormalLabel = thisNormalID + '_p' + nProbLabel
    # addTumorString = '  ' + newTumorLabel + ': [' + tProbLabel + ', results/sources/' + thisTumorID  + '-chrsOnly.dups_removed.bam]'
    # addNormalString = '  ' + newNormalLabel + ': [' + nProbLabel + ', results/sources/' + thisNormalID  + '-chrsOnly.dups_removed.bam]'

    # # Now add these lines to our running list, first tumor then normal
    # newSampleList.append(addTumorString)
    # newSampleList.append(addNormalString)   

    # # Now we need to make the line for the new mixture format, tumor then normal
    # # The pattern is mixID_TF0.XX: ["newTumorLabel", "newNormalLabel"]
    # newMixLabel = '  ' + thisMixID + '_TF' + thisTFLabel + ': ["' + newTumorLabel + '", "' + newNormalLabel + '"]'
    # newMixtureList.append(newMixLabel)

# Add a space to print out
print('\n')

#-----------------------------
#--- WRITE THE OUTPUT FILE ---
#-----------------------------

# Open the output file, config/samples2.yaml, for writing data
ofile = open(outfile, 'w')

# Write some data in either tumor or normal sample case
ofile.writelines('# SAMPLES2.YAML, created by python script: makeSamples2.py\n')
ofile.writelines('# Individual downsamples are listed below\n')
ofile.writelines('# Format is the following... sampleID : [downsample probability, file path]\n')
ofile.writelines('samples:\n')
for i in range(0, len(newSampleList) ):
  ofile.writelines(newSampleList[i] + '\n')

# Write the line break between dictionary groups
ofile.writelines('\n')

# Now write the mixture dictionary
ofile.writelines('# New mixture listings are listed below\n')
ofile.writelines('# Format is the following... mixtureID: [tumor partner ID, normal partner ID]\n')
ofile.writelines('mixtures:\n')
for i in range(0, len(newMixtureList) ):
  ofile.writelines(newMixtureList[i] + '\n')

# Now close the output file, done writing
ofile.close()
