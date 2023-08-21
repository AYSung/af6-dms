from Bio.Seq import Seq

F_PRIMER = 'GCATGCGGCGG'
R_PRIMER = 'CCACAGTCCCCGAGA'
VECTOR_BACKBONE = 'GAATTCTCAGTGGGCAGAGCGCACATCGC'
DMS_REGION = 'CTGCCCGGGCCGGAGGTGTCTGGGCGGAGCGTGGCTGCGGCCAGCGGACCGGGCGCCTGGGGCACTGACCACTACTGCCTGGAGCTGCTGCGGAAACGGGATTATGAAGGTTATTTATGCTCCCTGCTGCTCCCTGCAGAATCCCGAAGCTCTGTTTTTGCACTGAGGGCCTTTAATGTGGAACTGGCTCAGGTTAAAGACTCAGTCTCTGAGAAAACAATTGGACTGATGCGAATGCAGTTTTGGAAAAAAACTGTGGAAGATATATACTGTGACAATCCACCACATCAGCCTGTGGCCATTGAACTATGGAAGGCTGTTAAAAGACATAATCTGACTAAAAGATGGCTTATGAAAATCGTCGATGAAAGAGAAAAAAATCTGGATGACAAAGCATATCGTAATATCAAGGAACTGGAAAATTATGCTGAAAACACACAGAGCTCTCTTCTTTACTTAACACTAGAAATATTGGGTATAAAGGATCTTCATGCAGATCATGCTGCAAGTCATATTGGAAAAGCACAAGGCATTGTCACTTGCTTGAGAGCAACACCATATCATGGGAGCAGAAGAAAGGTGTTCCTTCCCATGGATATTTGTATGCTGCATGGTGTTTCACAAGAGGACTTTCTACGGAGGAACCAAGATAAAAATGTGAGAGATGTAATATATGACATTGCCAGTCAAGCACACTTGCACCTAAAGCATGCTAGGTCCTTTCACAAAACTGTTCCTGTGAAAGCATTTCCTGCTTTTCTTCAGACGGTTTCTCTAGAGGACTTTCTAAAGAAAATTCAGCGAGTGGATTTTGATATATTCCACCCATCTTTACAGCAGAAGAATACATTACTTCCATTATATTTGTATATTCAGTCATGGAGAAAAACATATTGA'


INDEXED_WT = f'{DMS_REGION[-9:]}CT{VECTOR_BACKBONE[:9]}'.upper()
INDEXED_WT_RC = str(Seq(INDEXED_WT).reverse_complement())
WT_SEQUENCE = f'{F_PRIMER}{DMS_REGION}{VECTOR_BACKBONE}{R_PRIMER}'
WT_SEQUENCE_RC = str(Seq(WT_SEQUENCE).reverse_complement())

REFERENCE = (
    f'{WT_SEQUENCE}{WT_SEQUENCE}{WT_SEQUENCE_RC}{WT_SEQUENCE_RC}{WT_SEQUENCE[:300]}'
)
REFERENCE_SEQUENCE = Seq(REFERENCE)
