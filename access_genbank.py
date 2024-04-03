from Bio import SeqIO, Entrez
import pandas as pd
import os.path as ospath

Entrez.email = 'contact@genorobotics.org'

#-------------------------------
#extract information from ncbi
#-------------------------------

def build_search_term(gene_name, min_len=0, max_len=-1):
    """ generate a search query for NCBI's nucleotide database (GenBank). Bio.Entrez does not provide filters such as sequence length and gene name
    so a search query keywords must be used instead.
    Used to make database queries in extract_database and download_database
    
    Parameters
    ----------
        gene_name: (str) name of the gene of interest
        min_len: (int, optional) lower bound of the sequence length filter, by default 0
        max_len: (int, optional) upper bound of the sequence length filter, by default -1 means no limit
    
    Returns
    ----------
        term: (str) search query
    """
    term = f"{gene_name}[Gene Name]"
    if max_len == -1:
        term += f" NOT 0:{min_len}[Sequence Length]"
    else:
        term += f" AND {min_len}:{max_len}[Sequence Length]"
    return term


def extract_database(gene_name, seqlen_start=0, seqlen_stop=-1):
    '''Request data on GenBank from NCBI (with Entrez module), with filtering options

    Parameters
    ----------
        gene_name: (str) name of the gene to extract from GenBank
        seqlen_start: (int, optional) lower bound of the sequence length filter, by default 0
        seqlen_stop: (int, optional) upper bound of the sequence length filter, by default -1 means no limit

    Returns
    ----------
        seq_record: (str) raw data of all sequences
    '''
    handle = Entrez.esearch(db= "nucleotide", term= build_search_term(gene_name, seqlen_start, seqlen_stop), retmax = 10000)
    record = Entrez.read(handle)
    matching_requests = record["IdList"]
    request_counter = 0
    while len(record["IdList"]) == 10000:
        request_counter += 1
        handle = Entrez.esearch(db= "nucleotide", term= build_search_term(gene_name, seqlen_start, seqlen_stop)
                                , retstart = request_counter*10000,  retmax = 10000)
        record = Entrez.read(handle)
        matching_requests += record["IdList"]
    seq_handle = Entrez.efetch(db="nucleotide", id=matching_requests, retmode = "fasta", rettype = "fasta")
    seq_record = seq_handle.read()
    return seq_record


def download_database(gene_name, seqlen_start=0, seqlen_stop=-1):
    '''Download data to a .fasta file, with filtering options

    Args:
        gene_name: (str) name of the gene to extract from GenBank
        seqlen_start: (int) lower bound of the sequence length filter
        seqlen_stop: (int) upper bound of the sequence length filter

    Returns:
        data_path: (str) path to downloaded data .fasta file
    '''
    data_path = f"{gene_name}[{seqlen_start},{seqlen_stop}].fasta"
    with open(data_path, mode="w") as file:
        file.write(extract_database(gene_name, seqlen_start, seqlen_stop))
    return data_path

def download_sequence(species, gene_name, dst, start_length=None, stop_length= None, id = None, permissive_search = True):
    """
    download sequence from GenBank through the Entrez database. 

    Parameters:
    ----------
    species(str): name of species
    gene_name(str): name of gene
    dst(str,Path-like): destination file path
    start_length(int): minimum length of sequence
    stop_length(int): maximum length of sequence
    id(list): list of NCBi ids of sequences to download. If provided, overrides gene_name and species.
    permissive_search(bool, default = True): when True, if Advanced NCBI query returns nothing, replace it with a less precise general query.
    """
    
    if id == None:
        search_term = f"{gene_name}[Gene Name] AND {species}[Organism]"
        if start_length!= None or stop_length!= None:
                search_term += f" {start_length}:{stop_length}[Sequence Length]"
        handle = Entrez.esearch(db = "nucleotide", term = search_term, retmax = 1, idtype = "acc")
        search_result = Entrez.read(handle)
        handle.close()
        id = search_result["IdList"]
        n=0
        for i in id:
            n+=1
        if n==0 and permissive_search:
            search_term = f"{gene_name} {species} {start_length}:{stop_length}[Sequence Length]"
            handle = Entrez.esearch(db = "nucleotide", term = search_term, retmax = 1, idtype = "acc")
            search_result = Entrez.read(handle)
            handle.close()
            id = search_result["IdList"]
    n=0
    for i in id:
        n+=1
    if n==1:
        handle = Entrez.efetch(db="nucleotide", id=id, retmode = "fasta", rettype = "fasta")
        sequence = handle.read()
        handle.close()
        with open(dst, mode="a") as writer:
            writer.write(sequence)