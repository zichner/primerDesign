import csv
import xlrd
import xlwt
import subprocess
import ConfigParser
import sys
import os
import string
import random
from datetime import datetime


FILE_CHRCHECKFIRST = ""
PARAM_NUMBER_N = 100
UNIQUE_PRIMERS = set("")
NON_UNIQUE_PRIMERS = set("")
FORBIDDEN_PRIMERS = set("")

def main():
    """Perform primer design"""
    config = ConfigParser.ConfigParser()
    config.optionxform = str
    config.read(sys.argv[1])
    maxSeqDist = config.getint("PrimerDesign", "max_dist_for_sequencing")
    maxPriDist = config.getint("PrimerDesign", "max_pcr_product_size")
    global FILE_CHRCHECKFIRST
    FILE_CHRCHECKFIRST = "/tmp/primerDesign_chrToCheckFirst_" + str(os.getpid()) + "_" + ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5)) + ".txt"
    file = open(FILE_CHRCHECKFIRST, "w")
    file.write(config.get("General", "check_first_chr"))
    file.close()

    varList = read_variant_list_from_file(sys.argv[2])
    varListTmp = varList
    primerList = [[-1, -1, -1, -1] for i in range(len(varList))]
    for paramSet in [(int(maxSeqDist/3), 25), (int(maxSeqDist/3), 100), (maxSeqDist,50), (maxSeqDist,200), (maxSeqDist,2000), (int(maxPriDist/3),50), (int(maxPriDist/3),200), (int(maxPriDist/3),1000), (maxPriDist,50), (maxPriDist,200), (maxPriDist,1000), (maxPriDist,5000)]:
        if len(varListTmp) == 0: break
        primerTargetSize, primerNum = paramSet
        primerPairs = determine_primer_pairs(varListTmp, primerTargetSize, primerNum, config)
        select_primer_pairs(primerList, primerPairs, varListTmp, config.getint("PrimerDesign", "max_dist_for_sequencing"))
        varListTmp = []
        for i in range(len(primerList)):
            if ((primerTargetSize <= maxSeqDist) and (primerList[i][0] <= 1 or primerList[i][1] <= 1)) or ((primerTargetSize > maxSeqDist) and (primerList[i][0] < 1 or primerList[i][1] < 1)):
                varListTmp.append(varList[i])
    resultList1 = createResultList(varList, primerList, config.getint("PrimerDesign", "max_dist_for_sequencing"), 0)
    resultList2 = createResultList(varList, primerList, config.getint("PrimerDesign", "max_dist_for_sequencing"), 1)
    if len(sys.argv) == 4:
       write_primer_list_to_file(sys.argv[3], resultList1, resultList2)
    else:
        print_array(resultList1)
        print ""
        print_array(resultList2)
    os.remove(FILE_CHRCHECKFIRST)

def read_variant_list_from_file(fileName):
    """Read in a list of variants from file"""
    varRows = []
    if fileName.split(".")[-1] == "xls":
        workbook = xlrd.open_workbook(fileName)
        worksheet = workbook.sheet_by_index(0)
        for i in range(worksheet.nrows):
            varRows.append([str(v) for v in worksheet.row_values(i)])
    else:
        varListFH = open(fileName, 'rb')
        varRows = list(csv.reader(varListFH, delimiter='\t', quoting=csv.QUOTE_NONE))
        varListFH.close()
        
    varList = []
    iter = 1
    for row in varRows:
        svStartChr, svStartPos, svEndChr, svEndPos, svType, svComment = row
        svStartPos = int(float(svStartPos))
        svEndPos = int(float(svEndPos))
        if svType in ["del", "dup", "inv3to3", "inv5to5", "trans3to3", "trans3to5", "trans5to3", "trans5to5", "snv"]:
            varList.append([iter, svStartChr, svStartPos, svEndChr, svEndPos, svType, svComment])
            iter += 1
        elif svType in ["inv",  "invRef"]:
            varList.append([iter, svStartChr, svStartPos, svEndChr, svEndPos, "invRefA", svComment])
            iter += 1
            varList.append([iter, svStartChr, svStartPos, svEndChr, svEndPos, "invRefB", svComment])
            iter += 1
        elif svType == "invAlt":
            varList.append([iter, svStartChr, svStartPos, svEndChr, svEndPos, "invAltA", svComment])
            iter += 1
            varList.append([iter, svStartChr, svStartPos, svEndChr, svEndPos, "invAltB", svComment])
            iter += 1
        else:
            sys.stderr.write(svType+" is not an accepted svType!\n")
    return varList


def write_primer_list_to_file(fileName, resultList1, resultList2):
    """Write list of primers to file"""
    if fileName.split(".")[-1] == "xls":
        workbook = xlwt.Workbook()
        worksheet = workbook.add_sheet('AtLeastOneUniquePrimer')
        maxCharacters = [0]*len(resultList1[0])
        for i in range(len(resultList1)):
            primerResult = resultList1[i]
            for j in range(len(primerResult)):
                if is_number(primerResult[j]) and not (j in [9, 10, 15, 16]): primerResult[j] = float(primerResult[j])
                worksheet.write(i,j,primerResult[j])
                maxCharacters[j] = max(maxCharacters[j], len(str(primerResult[j])))
        for i in range(len(maxCharacters)):
            worksheet.col(i).width = 330 * (maxCharacters[i]+1)

        worksheet = workbook.add_sheet('BothPrimersUnique')
        maxCharacters = [0]*len(resultList1[0])
        for i in range(len(resultList2)):
            primerResult = resultList2[i]
            for j in range(len(primerResult)):
                if is_number(primerResult[j]) and not (j in [9, 10, 15, 16]): primerResult[j] = float(primerResult[j])
                worksheet.write(i,j,primerResult[j])
                maxCharacters[j] = max(maxCharacters[j], len(str(primerResult[j])))
        for i in range(len(maxCharacters)):
            worksheet.col(i).width = 330 * (maxCharacters[i]+1)
        workbook.save(fileName)
    else:
        primerWriterFH = open(fileName, 'wb')
        primerWriter = csv.writer(primerWriterFH, delimiter='\t')
        primerWriter.writerows(resultList1)
        primerWriter.writerow([])
        primerWriter.writerows(resultList2)
        primerWriterFH.close()


def determine_primer_pairs(varList, primerTargetSize, primerNum, config):
    """Determine a set of primer pairs for a list of variants"""
    minPrimerOffset = config.getint("PrimerDesign","min_primer_offset")
    targetSequences = []
    sys.stderr.write("Number of variants left: "+str(len(varList))+"\n")
    for variant in varList:
        primerTargetSeqFwd, primerTargetSeqRev = get_primer_target_sequence(*variant, primerTargetSize=primerTargetSize, primerOffset=minPrimerOffset, samtools=config.get("Programs","samtools"), genomeFile=config.get("General","reference_genome"))
        targetSequences.append([variant[0], primerTargetSeqFwd,  primerTargetSeqRev])
    primer3Input = generate_primer3_input(targetSequences, primerNum, config.getint("PrimerDesign","max_pcr_product_size")-(2*config.getint("PrimerDesign","min_primer_offset")), config.items("Primer3"))
    primer3Output, primer3OutputTable, primerSeqList = call_primer3(primer3Input, config.get("Programs","primer3"), config.getfloat("PrimerDesign","primer_qual_cutoff"))
    uniqPrimers = determine_unique_primers(set(primerSeqList), config.get("Programs", "blast"), config.get("General", "reference_genome"), config.getint("PrimerDesign","min_mismatches"), config.getint("PrimerDesign","min_mismatches_close3p"), config.getint("PrimerDesign","min_dist3p"), config.getint("PrimerDesign","max_mismatches_forbidden"),  config.get("Programs","num_cpus"), config.get("PrimerDesign","word_size"))
    sys.stderr.write("unique primers: "+str(len(UNIQUE_PRIMERS))+"\n")
    sys.stderr.write("non-unique primers: "+str(len(NON_UNIQUE_PRIMERS))+"\n")
    sys.stderr.write("forbidden primers: "+str(len(FORBIDDEN_PRIMERS))+"\n")
    iter = 0
    for primer3OutputSubTable in primer3OutputTable:
        for primerPair in primer3OutputSubTable:
            if ((primerPair[2] in FORBIDDEN_PRIMERS) or (primerPair[3] in FORBIDDEN_PRIMERS)):
                primerPair.extend([False, False])
            else:
                primerPair.append(primerPair[2] in uniqPrimers)
                primerPair.append(primerPair[3] in uniqPrimers)
            primerPair.append(len(targetSequences[iter][1]) - int(primerPair[4].split(",")[0]) + minPrimerOffset - 1)
            primerPair.append(int(primerPair[5].split(",")[0]) - len(targetSequences[iter][1]) - PARAM_NUMBER_N + minPrimerOffset)
        iter += 1
    return primer3OutputTable


def select_primer_pairs(primerList, primerPairs, varList, maxDistForSequencing):
    """Select suitable primer pairs and add to list of variants"""
    iter = -1
    for primerPairSet in primerPairs:
        iter += 1
        if len(primerPairSet) == 0:
           continue
        numSeqBp1 = primerList[int(primerPairSet[0][0])-1][0]
        numSeqBp2 = primerList[int(primerPairSet[0][0])-1][1]
        primerScore1, primerScore2 = (1000.0, 1000.0)
        if numSeqBp1 > -1:
            primerScore1 = float(primerList[int(primerPairSet[0][0])-1][2][1])
        if numSeqBp2 > -1:
            primerScore2 = float(primerList[int(primerPairSet[0][0])-1][3][1])
        for primerPair in primerPairSet:
            if (not primerPair[8]) and (not primerPair[9]):
                continue
            if (varList[iter][5] == "dup" and (primerPair[10]+primerPair[11] > varList[iter][4]-varList[iter][2])):
                continue
            numSeqBp = 0
            if primerPair[10] < maxDistForSequencing:
                numSeqBp += 1
            if primerPair[11] < maxDistForSequencing:
                numSeqBp += 1
            if (primerPair[8] or primerPair[9]) and (numSeqBp > numSeqBp1):
                primerList[int(primerPair[0])-1][2] = primerPair[:]
                primerList[int(primerPair[0])-1][0] = numSeqBp
                primerScore1 = float(primerPair[1])
                numSeqBp1 = numSeqBp
            if (primerPair[8] and primerPair[9]) and (numSeqBp > numSeqBp2):
                primerList[int(primerPair[0])-1][3] = primerPair[:]
                primerList[int(primerPair[0])-1][1] = numSeqBp
                primerScore2 = float(primerPair[1])
                numSeqBp2 = numSeqBp


def createResultList(varList, primerList, maxDistForSequencing, primerSetId=0):
    """Combines the variant and the primer list and computes expected band sizes"""
    naString = "."
    resultList = []
    for i in range(len(varList)):
        resultList.append(varList[i][1:])
        if primerList[i][primerSetId] > -1:
            resultList[i].extend(primerList[i][primerSetId+2][1:4])
            resultList[i].extend(primerList[i][primerSetId+2][6:12])
        else:
            resultList[i].extend([naString]*13)
        
    for i in range(len(varList)):
        if primerList[i][primerSetId] > -1:
            if resultList[i][4] == "del" or resultList[i][4] == "snv":
                resultList[i].extend([resultList[i][1]-resultList[i][13], resultList[i][3]+resultList[i][14]])
                resultList[i].extend([resultList[i][13]+resultList[i][14], resultList[i][16]-resultList[i][15]])
            elif resultList[i][4] == "dup":
                resultList[i].extend([resultList[i][1]+resultList[i][13], resultList[i][3]-resultList[i][14], resultList[i][13]+resultList[i][14]])
                if (resultList[i][16]<resultList[i][15]):
                    resultList[i].append(resultList[i][15]-resultList[i][16])
                else:
                    resultList[i].append(naString)
            elif resultList[i][4] == "inv3to3":
                resultList[i].extend([resultList[i][1]-resultList[i][13], resultList[i][3]-resultList[i][14], resultList[i][13]+resultList[i][14], naString])
            elif resultList[i][4] == "inv5to5":
                resultList[i].extend([resultList[i][1]+resultList[i][13], resultList[i][3]+resultList[i][14], resultList[i][13]+resultList[i][14], naString])
            elif resultList[i][4] == "trans3to3":
                resultList[i].extend([resultList[i][1]-resultList[i][13], resultList[i][3]-resultList[i][14], resultList[i][13]+resultList[i][14], naString])
            elif resultList[i][4] == "trans3to5":
                resultList[i].extend([resultList[i][1]-resultList[i][13], resultList[i][3]+resultList[i][14], resultList[i][13]+resultList[i][14], naString])
            elif resultList[i][4] == "trans5to3":
                resultList[i].extend([resultList[i][1]+resultList[i][13], resultList[i][3]-resultList[i][14], resultList[i][13]+resultList[i][14], naString])
            elif resultList[i][4] == "trans5to5":
                resultList[i].extend([resultList[i][1]+resultList[i][13], resultList[i][3]+resultList[i][14], resultList[i][13]+resultList[i][14], naString])
            elif resultList[i][4] == "invRefA":
                if resultList[i+1][13] == naString:
                    resultList[i] = varList[i][1:]
                    resultList[i].extend([naString]*13)
                else:
                    resultList[i].extend([resultList[i][1]-resultList[i][13], resultList[i][1]+resultList[i][14], resultList[i][13]+resultList[i][14], resultList[i][13]+resultList[i+1][13]])
            elif resultList[i][4] == "invRefB":
                if resultList[i-1][13] == naString:
                    resultList[i] = varList[i][1:]
                    resultList[i].extend([naString]*13)
                else:
                    resultList[i].extend([resultList[i][3]-resultList[i][13], resultList[i][3]+resultList[i][14], resultList[i][13]+resultList[i][14], resultList[i-1][14]+resultList[i][14]])
            elif resultList[i][4] == "invAltA":
                if resultList[i+1][13] == naString:
                    resultList[i] = varList[i][1:]
                    resultList[i].extend([naString]*13)
                else:
                    resultList[i].extend([resultList[i][1]-resultList[i][13], resultList[i][3]-resultList[i][14], resultList[i][13]+resultList[i][14], resultList[i][13]+resultList[i+1][13]])
            elif resultList[i][4] == "invAltB":
                if resultList[i-1][13] == naString:
                    resultList[i] = varList[i][1:]
                    resultList[i].extend([naString]*13)
                else:
                    resultList[i].extend([resultList[i][1]+resultList[i][13], resultList[i][3]+resultList[i][14], resultList[i][13]+resultList[i][14], resultList[i-1][14]+resultList[i][14]])

    for i in range(len(varList)):
        row = resultList[i][0:9]
        row.extend(resultList[i][11:13])
        row.extend(resultList[i][9:11])
        row.extend(resultList[i][13:15])
        if resultList[i][13] != naString:
            row.extend([d<=maxDistForSequencing for d in resultList[i][13:15]])
        else:
            row.extend([naString]*2)
        row.extend(resultList[i][15:])
        resultList[i] = row
    return resultList


def generate_primer3_input(targetSequences, primerNum, maxPrimerDist, primer3Parameters):
    """Generate the input for primer3_core"""
    numberNs = PARAM_NUMBER_N
    primer3Input = ""
    for target in targetSequences:
        primer3Input += "SEQUENCE_ID=" + str(target[0]) + "\n"
        primer3Input += "SEQUENCE_TEMPLATE=" + target[1] + 'N'*numberNs + target[2] + "\n"
        primer3Input += "SEQUENCE_TARGET=" + str(len(target[1])) + "," + str(numberNs) + "\n"
        primer3Input += "PRIMER_PRODUCT_SIZE_RANGE=" + str(numberNs) + "-" + str(min(len(target[1]) + len(target[2]), maxPrimerDist) + numberNs) + "\n"
        primer3Input += "PRIMER_TASK=generic\n"
        primer3Input += "PRIMER_PICK_LEFT_PRIMER=1\n"
        primer3Input += "PRIMER_PICK_INTERNAL_OLIGO=0\n"
        primer3Input += "PRIMER_PICK_RIGHT_PRIMER=1\n"
        primer3Input += "PRIMER_MAX_NS_ACCEPTED=1\n"
        primer3Input += "PRIMER_NUM_RETURN="+ str(primerNum) + "\n"
        primer3Input += "PRIMER_EXPLAIN_FLAG=1\n"
        for param in primer3Parameters:
            primer3Input += param[0] + "=" + param[1] + "\n"
        primer3Input += "=\n"
    return primer3Input


def call_primer3(primer3Input, primer3Program, maxPrimerPenalty):
    """Call primer3_core"""
    p1 = subprocess.Popen(primer3Program, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    primer3Output = p1.communicate(input=primer3Input)[0]
    primer3OutputTable = []
    primer3OutputSubTable = []
    primerSeqList = []
    numValues = 0
    for line in primer3Output.splitlines():
        value = line.split("=")[-1]
        if ("SEQUENCE_ID=" in line):
            variantID = value
        if ("PRIMER_PAIR_" in line) and ("_PENALTY" in line):
            primerPenalty = value
            numValues += 1
        elif ("PRIMER_LEFT_" in line) and ("_SEQUENCE" in line):
            primerSeqLeft = value
            primerSeqList.append(value)
            numValues += 1
        elif ("PRIMER_RIGHT_" in line) and ("_SEQUENCE" in line):
            primerSeqRight = value
            primerSeqList.append(value)
            numValues += 1
        elif ("PRIMER_LEFT_" in line) and ("," in line) and not ("EXPLAIN" in line):
            primerPosLeft = value
            numValues += 1
        elif ("PRIMER_RIGHT_" in line) and ("," in line) and not ("EXPLAIN" in line):
            primerPosRight = value
            numValues += 1
        elif ("PRIMER_LEFT_" in line) and ("TM" in line):
            primerTmLeft = value
            numValues += 1
        elif ("PRIMER_RIGHT_" in line) and ("TM" in line):
            primerTmRight = value
            numValues += 1

        if numValues == 7:
            if float(primerPenalty) < maxPrimerPenalty:
                primer3OutputSubTable.append([variantID, primerPenalty, primerSeqLeft, primerSeqRight, primerPosLeft, primerPosRight, primerTmLeft, primerTmRight])
            numValues = 0
        if line == "=":
            primer3OutputTable.append(primer3OutputSubTable)
            primer3OutputSubTable = []

    return primer3Output, primer3OutputTable, primerSeqList


def run_blast(sequences, eValue, blast, genomeFile, numThreads, wordSize, seqidListFile=""):
    """Run blast for a set of sequences"""
    blastInput = ""
    iter = 0
    for primerSeq in sequences:
        iter += 1
        blastInput += ">" + primerSeq + "\n"
        blastInput += primerSeq + "\n"
    if seqidListFile == "":
        p1 = subprocess.Popen([blast, "-task", "blastn", "-db", genomeFile, "-evalue", str(eValue), "-num_threads", numThreads, "-outfmt", "6 std gaps nident", "-dust", "no", "-gapopen", "4", "-gapextend", "2", "-penalty", "-2", "-reward", "2", "-word_size", wordSize, "-max_target_seqs", "500"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    else:
        p1 = subprocess.Popen([blast, "-task", "blastn", "-db", genomeFile, "-seqidlist", seqidListFile, "-evalue", str(eValue), "-num_threads", numThreads, "-outfmt", "6 std gaps nident", "-dust", "no", "-gapopen", "4", "-gapextend", "2", "-penalty", "-2", "-reward", "2", "-word_size", wordSize, "-max_target_seqs", "2"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    return p1.communicate(input=blastInput)[0]


def determine_nonunique_sequences_from_blastoutput(blastOutput, minMismatches, minMismatchesClose3p, minDist3p, maxMismatchesToBeForbidden):
    """Determine non-unique sequences from blast output"""
    prevSeq = ""
    skipSeq = ""
    nonUniqueSequences = []
    completelyNonUniqueSequences = []
    for blastHit in blastOutput.splitlines():
        seq, chr, a, b, c, d, e, f, pos1, pos2, score1, score2, g, h = blastHit.split("\t")
        a = float(a)
        b = int(b)
        c = int(c)
        d = int(d)
        e = int(e)
        f = int(f)
        g = int(g)
        h = int(h)
        length = len(seq)
        if seq != skipSeq and ((h>length-minMismatchesClose3p) or (f<=length-minDist3p and h>length-minMismatches)):
            if seq==prevSeq:
                nonUniqueSequences.append(seq)
                if h>length-maxMismatchesToBeForbidden: completelyNonUniqueSequences.append(seq)
                skipSeq = seq
            else:
                prevSeq = seq
    return (nonUniqueSequences, completelyNonUniqueSequences)
    
    
def determine_unique_primers(primerList, blast, genomeFile, minMismatches, minMismatchesClose3p, minDist3p, maxMismatchesToBeForbidden, numThreads, wordSize):
    """Check primer uniqueness using blast"""
    sys.stderr.write(str(datetime.now())+" Primers to test: "+str(len(primerList))+"\n")
    global UNIQUE_PRIMERS
    global NON_UNIQUE_PRIMERS
    global FORBIDDEN_PRIMERS
    potUniquePrimers = primerList - NON_UNIQUE_PRIMERS - UNIQUE_PRIMERS - FORBIDDEN_PRIMERS
    for seqidListFile in [FILE_CHRCHECKFIRST, ""]:
        for eValue in [0.01, 0.1, 1, 10, 100]:
            if (seqidListFile == "" and eValue < 5):
                continue
            blastOutput = run_blast(potUniquePrimers, eValue, blast, genomeFile, numThreads, wordSize, seqidListFile)
            nonUniquePrimers, forbiddenPrimers = determine_nonunique_sequences_from_blastoutput(blastOutput, minMismatches, minMismatchesClose3p, minDist3p, maxMismatchesToBeForbidden)
            NON_UNIQUE_PRIMERS = NON_UNIQUE_PRIMERS | set(nonUniquePrimers)
            potUniquePrimers = potUniquePrimers - set(nonUniquePrimers)
            FORBIDDEN_PRIMERS = FORBIDDEN_PRIMERS | set(forbiddenPrimers)
            sys.stderr.write(str(datetime.now())+" potUnique "+str(eValue)+" "+str(len(potUniquePrimers))+"\n")
    UNIQUE_PRIMERS = UNIQUE_PRIMERS | potUniquePrimers
    return (potUniquePrimers | UNIQUE_PRIMERS)


def get_DNA_sequence(chr, start, end, samtools, genomeFile):
    """Get the DNA sequence of a given genomic region."""
    chrRegionString = create_chr_region_str(chr, start, end)
    p1 = subprocess.Popen([samtools, "faidx", genomeFile, chrRegionString], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["grep", "-v", ">"], stdin=p1.stdout, stdout=subprocess.PIPE)
    p3 = subprocess.Popen(["tr", "-d", "'\n'"], stdin=p2.stdout, stdout=subprocess.PIPE)
    sequence = p3.communicate()[0]
    return sequence


def create_chr_region_str(chr, start, end):
    """Create chromosomal region string"""
    return chr + ":" + str(start) + "-" + str(end)
    

def get_primer_target_sequence(id, svStartChr, svStartPos, svEndChr, svEndPos, svType, svComment, primerTargetSize, primerOffset, samtools, genomeFile):
    """Get the sequences in which primers will be placed"""
    if svType in ["del", "inv3to3", "trans3to3", "trans3to5", "snv", "invRefA", "invAltA"]:
        targetSeq1Start = svStartPos - primerOffset - primerTargetSize
        targetSeq1End = svStartPos - primerOffset
        targetSeq1 = get_DNA_sequence(svStartChr, targetSeq1Start, targetSeq1End, samtools, genomeFile).upper()
    elif svType in ["invRefB"]:
        targetSeq1Start = max(svEndPos - primerOffset - primerTargetSize, svStartPos + primerOffset)
        targetSeq1End = svEndPos - primerOffset
        targetSeq1 = get_DNA_sequence(svStartChr, targetSeq1Start, targetSeq1End, samtools, genomeFile).upper()
    elif svType in ["trans5to3", "trans5to5"]:
        targetSeq1Start = svStartPos + primerOffset
        targetSeq1End = svStartPos + primerOffset + primerTargetSize
        targetSeq1 = reverseComplementSequence(get_DNA_sequence(svStartChr, targetSeq1Start, targetSeq1End, samtools, genomeFile).upper())
    elif svType in ["dup", "inv5to5", "invAltB"]:
        targetSeq1Start = svStartPos + primerOffset
        targetSeq1End = min(svStartPos + primerOffset + primerTargetSize, svEndPos - primerOffset)
        targetSeq1 = reverseComplementSequence(get_DNA_sequence(svStartChr, targetSeq1Start, targetSeq1End, samtools, genomeFile).upper())

    if svType in ["del", "inv5to5", "snv", "invRefB", "invAltB"]:
        targetSeq2Start = svEndPos + primerOffset
        targetSeq2End = svEndPos + primerOffset + primerTargetSize
        targetSeq2 = get_DNA_sequence(svStartChr, targetSeq2Start, targetSeq2End, samtools, genomeFile).upper()
    elif svType in ["invRefA"]:
        targetSeq2Start = svStartPos + primerOffset
        targetSeq2End = min(svStartPos + primerOffset + primerTargetSize, svEndPos - primerOffset)
        targetSeq2 = get_DNA_sequence(svStartChr, targetSeq2Start, targetSeq2End, samtools, genomeFile).upper()
    elif svType in ["dup", "inv3to3", "invAltA"]:
        targetSeq2Start = max(svEndPos - primerTargetSize - primerOffset, svStartPos + primerOffset)
        targetSeq2End = svEndPos - primerOffset
        targetSeq2 = reverseComplementSequence(get_DNA_sequence(svStartChr, targetSeq2Start, targetSeq2End, samtools, genomeFile).upper())
    elif svType in ["trans3to5", "trans5to5"]:
        targetSeq2Start = svEndPos + primerOffset
        targetSeq2End = svEndPos + primerOffset + primerTargetSize
        targetSeq2 = get_DNA_sequence(svEndChr, targetSeq2Start, targetSeq2End, samtools, genomeFile).upper()
    elif svType in ["trans3to3", "trans5to3"]:
        targetSeq2Start = svEndPos - primerTargetSize - primerOffset
        targetSeq2End = svEndPos - primerOffset
        targetSeq2 = reverseComplementSequence(get_DNA_sequence(svEndChr, targetSeq2Start, targetSeq2End, samtools, genomeFile).upper())
    return (targetSeq1, targetSeq2)


def reverseComplementSequence(sequence):
    '''Reverse complements a DNA sequence'''
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'N': 'N', 'n': 'n'} 
    letters = list(sequence)
    letters = [basecomplement[base] for base in letters]
    complement = ''.join(letters)
    return complement[::-1]


def print_array(arrayVariable):
    """Print 2D array to terminal"""
    print "\n".join(["\t".join(map(str, r)) for r in arrayVariable])


def is_number(s):
    """Check whether string is a number"""
    try:
        float(s)
        return True
    except ValueError:
        return False


sys.stderr.write(str(datetime.now())+"\n")
main()
sys.stderr.write(str(datetime.now())+"\n")
