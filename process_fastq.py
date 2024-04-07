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
    :param folder: folder containing the .fastq files (usually fastq_pass folder) [str]
    :param dst: destination file (.fastq) [str]
    :return: None
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


def extract_fastq(main_dir, sample_nb, dst):
    """
    extract all fastq files under a certain barcode/sample number from an expedition results folder.

    Parameters
    ----------
    main_dir: expedition folder path
    sample_nb: sample number for which to extract fastq. ie 6 to extract from folder barcode06
    dst: destination file path
    """
    for root, dirs, files in os.walk(main_dir):
        if "fastq_pass" in dirs:
            pass_dir = ospath.join(root, "fastq_pass")
            break
    for root, dirs, files in os.walk(pass_dir):
        if f"barcode{sample_nb}" in dirs:
            fastq_dir = ospath.join(root, f"barcode{sample_nb}")
        elif f"barcode0{sample_nb}" in dirs:
            fastq_dir = ospath.join(root, f"barcode0{sample_nb}")
    concatenate_fastq(fastq_dir, dst)

#----------------------------------
#GET AND VISUALIZE QUALITY SCORES
#----------------------------------

def getContentFile(pathFastqFile: str):

    # open file and get all the reads
    allReadsAsString = []

    with open(pathFastqFile) as fastq:
        allReadsAsString = fastq.readlines()

    return allReadsAsString


def extractReadQuality(fileContent):

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
    :param length: length of the sequence generated [int]
    :param seed: random seed [int]
    :return: random sequence [str]
    """
    if seed is not None: random.seed(seed)
    return "".join(["ATGC"[random.randint(0,3)] for i in range(length)])

def damage_sequence(sequence, mutation_rate=0.05, deletion_rate=0.05, insertion_rate=0.05):
    """
    Damage a sequence randomly with mutation (substitution), deletion and insertion
    Warning: to avoid infinite loop, don't set the insertion_rate to 1.0 (or near)
    :param sequence: original sequence to damage [str]
    :param mutation_rate: mutation rate (between 0.0 and 1.0) [float]
    :param deletion_rate: deletion_rate (between 0.0 and 1.0) [float]
    :param insertion_rate: insertion_rate (between 0.0 and 1.0) [float]
    :return: damaged sequence [str]
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

