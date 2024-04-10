import process_fastq
from Bio import SeqIO, Align

#In order to run the tests, enter in terminal : pytest test_utils.py
#Warning : the 3 last tests have an extremely small chance to fail due to functions based on randomness, running the tests several times can confirm if it is the reason or not

def test_read_fastq():
    filePath = "test_material/rbcL_Qiagen_tomato_5000_test_file.fastq"
    reads = process_fastq.read_fastq(filePath)
    assert len(reads) == 2
    assert str(reads[0].seq) == "GTATGCTTCGTTCATTCAAATTTGGGTGTTTAGCAATTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTGGATCGGCAGGTGTTTAACCGTTTTCGCATTTATCGTGGCTTAACACCAGCGTTTTTCATTGTAACTTAAAATTATATAAGAGAAGAAGAATCTTTGATTTTTTTTTTTGAAAAAGGTAACCGAGCTTCTTTAGTAATAAGACTATTCAAATTACAATATTCGTGGAAAATCGTAATAAATATTGAAGGCATCTTTTAATAGCGAAGTTTGAACAAAATTTCCAA"
    assert str(reads[1].seq) == "GTTGTACTTCGTTCCAGTTATCAGATGTTGGGTGTTTAGCCGTTTTCGCATTTATCATTGAAACAACCGCGTTTTCGTGCGCCGCTTCACCTACAATGGAAGTAAACATATTGGTAACAGAACCTTTATGTAAAGGTCTAAAGGAGTGGCTACATAAGCAATATATTGATCTTTTCTCCAACAACGCGCTCGATGCGGTAGCATCGCACCTTTGTAACGATCAAGACTGGTAGAGTCCATCGGTCCATACAATTGTCCATGTACCAGTAGAAGATTCGGCTGCGGCCCCTGCTACTGGGTGGAACTCCAGGTTGAGGTTACTCGGAATGCTAATATATCAGTATCCTTGGTTTGGTACTCAGGAGTATAATAAGTCAATTTGTACTCTTTAACACCAGCTTTGAATCAACACTTTAGTCTCTG"
    
def test_split_fastq():
    inputPath="test_material/rbcL_Qiagen_tomato_5000_test_file.fastq"
    outputDir="test_material"
    output_name="outputTest"
    assert process_fastq.split_fastq(inputPath,outputDir,output_name) == ('test_material\\outputTest_top20.fastq', 'test_material\\outputTest_remaining80.fastq')

def test_quality_scores():
    filePath = "test_material/rbcL_Qiagen_tomato_5000_test_file.fastq"
    reads = process_fastq.getContentFile(filePath)
    assert reads==['@632f5334-2c26-4f91-a611-bbcc85d83ea1 runid=84dc383f16e4874203b6120bf223ea75d3815a54 read=14 ch=417 start_time=2021-11-26T13:56:48Z flow_cell_id=FAR79813 protocol_group_id=rbcL_Qiagen_tomato sample_id=no_sample\n', 'GTATGCTTCGTTCATTCAAATTTGGGTGTTTAGCAATTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTGGATCGGCAGGTGTTTAACCGTTTTCGCATTTATCGTGGCTTAACACCAGCGTTTTTCATTGTAACTTAAAATTATATAAGAGAAGAAGAATCTTTGATTTTTTTTTTTGAAAAAGGTAACCGAGCTTCTTTAGTAATAAGACTATTCAAATTACAATATTCGTGGAAAATCGTAATAAATATTGAAGGCATCTTTTAATAGCGAAGTTTGAACAAAATTTCCAA\n', '+\n', "$$%&$$%00*268%%)($$#####%++]<;7-$$#$#$257166>>8;7820.(''''')668.0+,-0.,(()''(,,++55]<<*')10024</8&.479;00023)(())&%%$%#$#%'02UFF=)((&$#(,,*%&20577C]]/-.]]-970,+05*+(..1049;:9853453+*/..-.*&*%%&%'/*7?<9=;&&&*/10-.1/1:?A:9:71008::B]]>=//.*,*(),/)+//00347%&'%#**&%()'&'(((()'***+*33373356*,+')),0,)&\n", '@c4124d21-64a3-4b8e-ac5a-771ea99ff56f runid=84dc383f16e4874203b6120bf223ea75d3815a54 read=704 ch=203 start_time=2021-11-26T14:00:32Z flow_cell_id=FAR79813 protocol_group_id=rbcL_Qiagen_tomato sample_id=no_sample\n', 'GTTGTACTTCGTTCCAGTTATCAGATGTTGGGTGTTTAGCCGTTTTCGCATTTATCATTGAAACAACCGCGTTTTCGTGCGCCGCTTCACCTACAATGGAAGTAAACATATTGGTAACAGAACCTTTATGTAAAGGTCTAAAGGAGTGGCTACATAAGCAATATATTGATCTTTTCTCCAACAACGCGCTCGATGCGGTAGCATCGCACCTTTGTAACGATCAAGACTGGTAGAGTCCATCGGTCCATACAATTGTCCATGTACCAGTAGAAGATTCGGCTGCGGCCCCTGCTACTGGGTGGAACTCCAGGTTGAGGTTACTCGGAATGCTAATATATCAGTATCCTTGGTTTGGTACTCAGGAGTATAATAAGTCAATTTGTACTCTTTAACACCAGCTTTGAATCAACACTTTAGTCTCTG\n', '+\n', "%84695532.,-.:+&''+'&(%'('&%$&&..]2233-,-58>DB]479;>9899)))..-)((-0'(7879<?<,+')]67767=:/.//-./2243632221/'&',.,&&((%&'().,+')'$'$##,04434-62060-('&$$())+7:7967<L]](*(>=5;8?A>:7@67=<663/.-3-3'.5=;2++;9:;99&&()'%*%%+++,,/>@9?@HABA334+)*&)+-.0/25:750**)(()23581100/0-/44=8]2@3333101334214778743.+'(%%.45<D<8:5/*6;9942*6062069(58762226-2/]]=>''':6:.,344095-+,*)())()(&&)--,.**,+.&%/7@>;2012441+*+))68123?=660-=&56')&522688<]1+"]
    quality = process_fastq.extractReadQuality(reads)
    assert quality==["$$%&$$%00*268%%)($$#####%++]<;7-$$#$#$257166>>8;7820.(''''')668.0+,-0.,(()''(,,++55]<<*')10024</8&.479;00023)(())&%%$%#$#%'02UFF=)((&$#(,,*%&20577C]]/-.]]-970,+05*+(..1049;:9853453+*/..-.*&*%%&%'/*7?<9=;&&&*/10-.1/1:?A:9:71008::B]]>=//.*,*(),/)+//00347%&'%#**&%()'&'(((()'***+*33373356*,+')),0,)&", "%84695532.,-.:+&''+'&(%'('&%$&&..]2233-,-58>DB]479;>9899)))..-)((-0'(7879<?<,+')]67767=:/.//-./2243632221/'&',.,&&((%&'().,+')'$'$##,04434-62060-('&$$())+7:7967<L]](*(>=5;8?A>:7@67=<663/.-3-3'.5=;2++;9:;99&&()'%*%%+++,,/>@9?@HABA334+)*&)+-.0/25:750**)(()23581100/0-/44=8]2@3333101334214778743.+'(%%.45<D<8:5/*6;9942*6062069(58762226-2/]]=>''':6:.,344095-+,*)())()(&&)--,.**,+.&%/7@>;2012441+*+))68123?=660-=&56')&522688<]1"]
    converted_quality = process_fastq.convertQuality(quality)
    assert converted_quality==[3, 3, 4, 5, 3, 3, 4, 15, 15, 9, 17, 21, 23, 4, 4, 8, 7, 3, 3, 2, 2, 2, 2, 2, 4, 10, 10, 60, 27, 26, 22, 12, 3, 3, 2, 3, 2, 3, 17, 20, 22, 16, 21, 21, 29, 29, 23, 26, 22, 23, 17, 15, 13, 7, 6, 6, 6, 6, 6, 8, 21, 21, 23, 13, 15, 10, 11, 12, 15, 13, 11, 7, 7, 8, 6, 6, 7, 11, 11, 10, 10, 20, 20, 60, 27, 27, 9, 6, 8, 16, 15, 15, 17, 19, 27, 14, 23, 5, 13, 19, 22, 24, 26, 15, 15, 15, 17, 18, 8, 7, 7, 8, 8, 5, 4, 4, 3, 4, 2, 3, 2, 4, 6, 15, 17, 52, 37, 37, 28, 8, 7, 7, 5, 3, 2, 7, 11, 11, 9, 4, 5, 17, 15, 20, 22, 22, 34, 60, 60, 14, 12, 13, 60, 60, 12, 24, 22, 15, 11, 10, 15, 20, 9, 10, 7, 13, 13, 16, 15, 19, 24, 26, 25, 24, 23, 20, 18, 19, 20, 18, 10, 9, 14, 13, 13, 12, 13, 9, 5, 9, 4, 4, 5, 4, 6, 14, 9, 22, 30, 27, 24, 28, 26, 5, 5, 5, 9, 14, 16, 15, 12, 13, 16, 14, 16, 25, 30, 32, 25, 24, 25, 22, 16, 15, 15, 23, 25, 25, 33, 60, 60, 29, 28, 14, 14, 13, 9, 11, 9, 7, 8, 11, 14, 8, 10, 14, 14, 15, 15, 18, 19, 22, 4, 5, 6, 4, 2, 9, 9, 5, 4, 7, 8, 6, 5, 6, 7, 7, 7, 7, 8, 6, 9, 9, 9, 10, 9, 18, 18, 18, 22, 18, 18, 20, 21, 9, 11, 10, 6, 8, 8, 11, 15, 11, 8, 5, 4, 23, 19, 21, 24, 20, 20, 18, 17, 13, 11, 12, 13, 25, 10, 5, 6, 6, 10, 6, 5, 7, 4, 6, 7, 6, 5, 4, 3, 5, 5, 13, 13, 60, 17, 17, 18, 18, 12, 11, 12, 20, 23, 29, 35, 33, 60, 19, 22, 24, 26, 29, 24, 23, 24, 24, 8, 8, 8, 13, 13, 12, 8, 7, 7, 12, 15, 6, 7, 22, 23, 22, 24, 27, 30, 27, 11, 10, 6, 8, 60, 21, 22, 22, 21, 22, 28, 25, 14, 13, 14, 14, 12, 13, 14, 17, 17, 19, 18, 21, 18, 17, 17, 17, 16, 14, 6, 5, 6, 11, 13, 11, 5, 5, 7, 7, 4, 5, 6, 7, 8, 13, 11, 10, 6, 8, 6, 3, 6, 3, 2, 2, 11, 15, 19, 19, 18, 19, 12, 21, 17, 15, 21, 15, 12, 7, 6, 5, 3, 3, 7, 8, 8, 10, 22, 25, 22, 24, 21, 22, 27, 43, 60, 60, 7, 9, 7, 29, 28, 20, 26, 23, 30, 32, 29, 25, 22, 31, 21, 22, 28, 27, 21, 21, 18, 14, 13, 12, 18, 12, 18, 6, 13, 20, 28, 26, 17, 10, 10, 26, 24, 25, 26, 24, 24, 5, 5, 7, 8, 6, 4, 9, 4, 4, 10, 10, 10, 11, 11, 14, 29, 31, 24, 30, 31, 39, 32, 33, 32, 18, 18, 19, 10, 8, 9, 5, 8, 10, 12, 13, 15, 14, 17, 20, 25, 22, 20, 15, 9, 9, 8, 7, 7, 8, 17, 18, 20, 23, 16, 16, 15, 15, 14, 15, 12, 14, 19, 19, 28, 23, 60, 17, 31, 18, 18, 18, 18, 16, 15, 16, 18, 18, 19, 17, 16, 19, 22, 22, 23, 22, 19, 18, 13, 10, 6, 7, 4, 4, 13, 19, 20, 27, 35, 27, 23, 25, 20, 14, 9, 21, 26, 24, 24, 19, 17, 9, 21, 15, 21, 17, 15, 21, 24, 7, 20, 23, 22, 21, 17, 17, 17, 21, 12, 17, 14, 60, 60, 28, 29, 6, 6, 6, 25, 21, 25, 13, 11, 18, 19, 19, 15, 24, 20, 12, 10, 11, 9, 8, 7, 8, 8, 7, 8, 7, 5, 5, 8, 12, 12, 11, 13, 9, 9, 11, 10, 13, 5, 4, 14, 22, 31, 29, 26, 17, 15, 16, 17, 19, 19, 16, 10, 9, 10, 8, 8, 21, 23, 16, 17, 18, 30, 28, 21, 21, 15, 12, 28, 5, 20, 21, 6, 8, 5, 20, 17, 17, 21, 23, 23, 27, 60, 16]

def test_create_random_sequence():
    seed=1
    random_seq = process_fastq.create_random_sequence(50,seed)
    assert random_seq == "TAGACCCCTACACCACGTAGAAAACTCATCCTGTTCGACATGAGCTGGCC"

def test_mutation():
    sequence = "ATCG"
    mutated_sequence = process_fastq.damage_sequence(sequence, mutation_rate=1.0)
    assert mutated_sequence != sequence

def test_deletion():
    sequence = "ATCG"
    deleted_sequence = process_fastq.damage_sequence(sequence, deletion_rate=1.0)
    assert deleted_sequence != sequence

def test_insertion():
    sequence = "ATCG"
    inserted_sequence = process_fastq.damage_sequence(sequence, insertion_rate=0.98)
    assert inserted_sequence != sequence