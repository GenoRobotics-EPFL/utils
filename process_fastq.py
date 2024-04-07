import glob, gzip, math, random
import numpy as np
from Bio import SeqIO, Align
import os.path as ospath

import os
import gzip
import shutil
import sys 
import matplotlib.pyplot as plt


#----------------------------
#DEAL WITH FASTQ FILES
#----------------------------

def extract_gz(src_file, dst_file):
    """
    Extracts a gzipped file.

    Args:
        src_file (str): Path to the source gzipped file.
        dst_file (str): Path to the destination file.

    Returns:
        None
    """
    with gzip.open(src_file, 'rb') as f_in:
        with open(dst_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)



def read_fastq(fastq_filepath=None):
    """
    Read a fastq file and return the reads.
    :param fastq_filepath: filepath of the .fastq file [str]
    :return: reads from the fastq file [list of Bio.SeqRecord.SeqRecord]
    """
    if fastq_filepath is None: fastq_filepath = "data/rbcL_Qiagen_tomato.fastq" # default path (example)
    if fastq_filepath.lower().endswith('.gz'):
        f = gzip.open(fastq_filepath, 'rt')
    else:
        f = open(fastq_filepath, 'rt')
    reads = []
    for read in SeqIO.parse(f, "fastq"):
        reads.append(read)
    return reads



def split_fastq(input_path, output_dir, base_name, percentile=20):
    """
    Split a FASTQ file into two separate files based on sequence length.

    Args:
        input_path (str): Path to the input FASTQ file.
        output_dir (str): Directory where the output files will be saved.
        base_name (str): Base name for the output files.
        percentile (int, optional): The percentile value to split the sequences. Defaults to 20.

    Returns:
        tuple: A tuple containing the paths to the top percentile sequences file and the remaining sequences file.
    """
    # Read sequences and sort by length
    sequences = list(SeqIO.parse(input_path, "fastq"))
    sequences.sort(key=lambda x: len(x), reverse=True)

    # Split into top percentile and remaining sequences
    split_index = len(sequences) * percentile // 100
    top_sequences = sequences[:split_index]
    remaining_sequences = sequences[split_index:]

    # Write to separate files
    top_sequences_path = os.path.join(output_dir, f"{base_name}_top{percentile}.fastq")
    remaining_sequences_path = os.path.join(output_dir, f"{base_name}_remaining{100-percentile}.fastq")
    SeqIO.write(top_sequences, top_sequences_path, "fastq")
    SeqIO.write(remaining_sequences, remaining_sequences_path, "fastq")

    return top_sequences_path, remaining_sequences_path


def concatenate_fastq(src_folder=None, dst=None):
    """
    Concatenate all .fastq from a folder into a single .fastq file.
    Args:
        folder (str): folder containing the .fastq files (usually fastq_pass folder)
        dst (str): destination file (.fastq)
    Returns:
        None
    """
    if src_folder is None: src_folder = "fastq_pass" # default folder (example)
    if dst is None: dst = f"{src_folder}/concatenation.fastq" # default destination (example)
    if dst[-6:] != ".fastq": dst = f"{dst}.fastq"
    def get_file_iterator(filename):
        if filename.lower().endswith('.gz'):
            f = gzip.open(filename, 'rt')
        else:
            f = open(filename, 'rt')
        return f
    fastq_list = [fastq for fastq in glob.glob(f"{src_folder}/*.fastq*")]
    fastq_iterators = [SeqIO.parse(get_file_iterator(fastq), "fastq") for fastq in fastq_list]
    while True:
        for fq in fastq_iterators:
            try:
                SeqIO.write(next(fq), open(dst,"at"), "fastq")
            except StopIteration:
                fastq_iterators.remove(fq)
        if len(fastq_iterators) == 0:
            break


def extract_fastq(main_dir, barcode_nb, dst):
    """
    Extract all fastq files under a certain barcode/sample number from an expedition results folder.

    Args:
        main_dir(str): expedition folder path
        sample_nb(int): sample number for which to extract fastq. ie 6 to extract from folder barcode06
        dst(str): destination file path
    Returns:
        None
    """
    search_dirs=[]
    for root, dirs, files in os.walk(main_dir):
        if "fastq_pass" in dirs:
            search_dirs.append(ospath.join(root, "fastq_pass"))
        if "barcoding" in dirs:
            search_dirs.append(ospath.join(root, "barcoding"))
        
    for search_dir in search_dirs:
        for root, dirs, files in os.walk(search_dir):
            if f"barcode{barcode_nb}" in dirs:
                fastq_dir = ospath.join(root, f"barcode{barcode_nb}")
                print(fastq_dir)
                concatenate_fastq(fastq_dir, dst)
            elif f"barcode0{barcode_nb}" in dirs:
                fastq_dir = ospath.join(root, f"barcode0{barcode_nb}")
                print(fastq_dir)
                concatenate_fastq(fastq_dir, dst)

#----------------------------------
#GET AND VISUALIZE QUALITY SCORES
#----------------------------------

def getContentFile(pathFastqFile):
    """
    Open file and get all the reads.

    Args:
        pathFastqFile(str): fastq file to read

    Returns:
        allReadsAsSring(str): all fastq faile reads
    """
    # open file and get all the reads
    allReadsAsString = []

    with open(pathFastqFile) as fastq:
        allReadsAsString = fastq.readlines()

    return allReadsAsString


def extractReadQuality(fileContent):
    """
    Extract the quality of the file reads. Note that the reads are concatenated in one list.

    Args:
        fileContent(str): fastq file content as a string

    Returns:
        qualityOfReads(list of str): quality scores as a string list
    """
    # the fileContent must contain multiple reads, where each read is 4 lines:
    if not (len(fileContent) % 4 == 0):
        raise ValueError("fastq file must have reads composed of 4 lines")

    qualityOfReads = []

    # the quality of each read is always the 4th element
    for i in range(3, len(fileContent), 4):

        line = fileContent[i]

        # remove the newline character at the end
        line = line[:-1]

        qualityOfReads.append(line)

    return qualityOfReads


def convertQuality(allReadsQuality):
    """
    Converts the qualities from fastq format to int.

    Args:
        allReadsQuality(list of str): quality scores in the original format

    Returns:
        qualityOfReads(list of int): quality scores converted in int
    """
    # note that here we can't distinguish the reads from each other anymore
    # they are all in one list
    convertedQualities = []

    for readQuality in allReadsQuality:

        for rawQuality in readQuality:

            # transform the character into the int quality
            score = ord(rawQuality) - 33

            convertedQualities.append(score)

    return convertedQualities


def visualiseQualities(pathFastqFile, readQualityConverted):
    """
    Creates a graph representing the frequency of each score value.

    Args:
        pathFastqFile(str): path of the fastq file in directory
        readQualityConverted(list of int): quality scores converted in int

    Returns:
        None
    """
    plt.figure(figsize=(10, 6))

    plt.hist(readQualityConverted, bins=range(0, 95))

    plt.ylabel("Frequency")

    plt.xlabel("Base quality score 0-93 (higher is better)")
    plt.xlim(0, 95)

    nameFile = os.path.basename(pathFastqFile)
    plt.title(nameFile)

    nameOutput = nameFile + ".png"
    plt.savefig(nameOutput)  # saves in current directory


#--------------------------------
#CREATE AND DAMAGE SEQUENCES
#--------------------------------
    
def create_random_sequence(length=500, seed=None):
    """
    Generate a random sequence.
    Args:
    length(int): length of the sequence generated
    seed(int): random seed 
    Returns:
        random sequence(str)
    """
    if seed is not None: random.seed(seed)
    return "".join(["ATGC"[random.randint(0,3)] for i in range(length)])

def damage_sequence(sequence, mutation_rate=0.05, deletion_rate=0.05, insertion_rate=0.05):
    """
    Damage a sequence randomly with mutation (substitution), deletion and insertion.
    Warning: to avoid infinite loop, don't set the insertion_rate to 1.0 (or near).

    Args:
        sequence(str): original sequence to damage
        mutation_rate(float): mutation rate (between 0.0 and 1.0)
        deletion_rate(float): deletion_rate (between 0.0 and 1.0)
        insertion_rate(float): insertion_rate (between 0.0 and 1.0)

    Returns:
        sequence(str): damaged sequence
    """
    if mutation_rate < 0 or mutation_rate > 1:
        raise Exception("[damage_sequence] mutation_rate is incorrect (must be between 0.0 and 1.0)")
    if deletion_rate < 0 or deletion_rate > 1:
        raise Exception("[damage_sequence] deletion_rate is incorrect (must be between 0.0 and 1.0)")
    if insertion_rate < 0 or insertion_rate > 1:
        raise Exception("[damage_sequence] insertion_rate is incorrect (must be between 0.0 and 1.0)")
    sequence = "".join(["ATGC".replace(b, '')[random.randint(0, 2)] if random.random() < mutation_rate else b for b in sequence]) # mutation / substitution
    sequence = "".join(['' if random.random() < deletion_rate else b for b in sequence]) # deletion
    insertion_extension_rate = insertion_rate  # can be changed if the extension rate is different
    def get_insert(extension_rate=insertion_extension_rate):
        insert = "ATGC"[random.randint(0,3)]
        while random.random() < extension_rate: insert += "ATGC"[random.randint(0,3)] # extension (after insertion)
        return insert
    sequence = "".join([sequence[i:i+1] + get_insert() if random.random() < insertion_rate else sequence[i:i+1] for i in range(len(sequence)+1)]) # insertion
    return sequence

