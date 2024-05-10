from Bio import SeqIO, Entrez
import pandas as pd
import os.path as ospath

Entrez.email = 'contact@genorobotics.org'

#-------------------------------
#extract information from ncbi
#-------------------------------

def build_search_term(gene_name, min_len=0, max_len=-1):
    """ 
    Generate a search query for NCBI's nucleotide database (GenBank). Bio.Entrez does not provide filters such as sequence length and gene name
    so a search query keywords must be used instead.
    Used to make database queries in extract_database and download_database
    
    Parameters: 
        gene_name(str): name of the gene of interest
        min_len(int, optional): lower bound of the sequence length filter, by default 0
        max_len(int, optional): upper bound of the sequence length filter, by default -1 means no limit
    
    Returns:
        term(str): search query
    """
    term = f"{gene_name}[Gene Name]"
    if max_len == -1:
        term += f" NOT 0:{min_len}[Sequence Length]"
    else:
        term += f" AND {min_len}:{max_len}[Sequence Length]"
    return term


def extract_database(gene_name, seqlen_start=0, seqlen_stop=-1):
    '''
    Request data on GenBank from NCBI (with Entrez module), with filtering options

    Parameters:
        gene_name(str): name of the gene to extract from GenBank
        seqlen_start(int, optional): lower bound of the sequence length filter, by default 0
        seqlen_stop(int, optional): upper bound of the sequence length filter, by default -1 means no limit

    Returns:
        seq_record(str): raw data of all sequences
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
    '''
    Download data to a .fasta file, with filtering options

    Args:
        gene_name(str): name of the gene to extract from GenBank
        seqlen_start(int): lower bound of the sequence length filter
        seqlen_stop(int): upper bound of the sequence length filter

    Returns:
        data_path(str): path to downloaded data .fasta file
    '''
    data_path = f"{gene_name}[{seqlen_start},{seqlen_stop}].fasta"
    with open(data_path, mode="w") as file:
        file.write(extract_database(gene_name, seqlen_start, seqlen_stop))
    return data_path

def download_sequence_from_id(id: list, dst = None):
    """
    download sequence from GenBank through the Entrez database using unique identifiers.

    Args:
        id (list): list of IDs of sequences to download from GenBank
        dst(str): destination file path

    Returns:
        Bio.SeqIO.FastaIO.FastaIterator: iterator over the records corresponding to the input ids
    
    """
    handle = Entrez.efetch(db="nucleotide", id=id, retmode = "fasta", rettype = "fasta")
    record_parser = SeqIO.parse(handle, "fasta")
    if dst != None:
        with open(dst, mode="a") as writer:
            SeqIO.write(record_parser,writer,"fasta")
        return SeqIO.parse(dst,"fasta")
    return record_parser

def get_species_from_ids(ids:list):
    """
    Return the species from GenBank entries with certain ids.

    Parameters:
        ids: list of ids to search

    Returns:
        dictionary that matches ids to organism names
    
    """
    species_dict= dict()
    for id in ids:
        handle = Entrez.efetch(db="nucleotide", id=id, retmode = "gb", rettype = "gb")
        records = SeqIO.parse(handle,"genbank")
        one_record = next(records)
        species_dict[id]=one_record.annotations["organism"]
    return species_dict

def download_sequence(species, gene_name, dst = None, start_length=None, stop_length= None, permissive_search = True):
    """
    Download sequence from GenBank through the Entrez database. 

    Parameters:
        species(str): name of species
        gene_name(str): name of gene
        dst(str,Path-like): destination file path, if None does not write to a file, only returns 
        start_length(int): minimum length of sequence
        stop_length(int): maximum length of sequence
        id(list): list of NCBi ids of sequences to download. If provided, overrides gene_name and species.
        permissive_search(bool, default = True): when True, if Advanced NCBI query returns nothing, replace it with a less precise general query.
    
    Returns:
        Bio.SeqRecord.SeqRecord: sequence found from GenBank

    """
    

    search_term = f"{gene_name}[Gene Name] AND {species}[Organism]"
    if start_length!= None or stop_length!= None:
            search_term += f" {start_length}:{stop_length}[Sequence Length]"
    handle = Entrez.esearch(db = "nucleotide", term = search_term, retmax = 1, idtype = "acc")
    search_result = Entrez.read(handle)
    handle.close()
    id = search_result["IdList"]
    if len(id)==0 and permissive_search:
        search_term = f"{gene_name} {species} {start_length}:{stop_length}[Sequence Length]"
        handle = Entrez.esearch(db = "nucleotide", term = search_term, retmax = 1, idtype = "acc")
        search_result = Entrez.read(handle)
        handle.close()
        id = search_result["IdList"]
    if len(id)==1:
        return next(download_sequence_from_id(id,dst))
    else: 
         raise ValueError("Your search did not return any results")
