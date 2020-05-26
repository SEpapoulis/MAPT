#from .. import pna_designer
import datetime
import sys
sys.path.append("..")

import pna_designer


#gathering sequence ids
#taxid = sil.db.get_taxon_byname('Fungi')
taxid = 6516 #Fungi
#taxid = 25326 #Microcystis
silva = pna_designer.silva_manager("silva_db_compiled/")

fungi_taxpath = silva.find_taxpath("Fungi")

F='GTGYCAGCMGCCGCGGTAA'
R='CCGYCAATTYMTTTRAGTTT'


pd=pna_designer.PNA_Designer(result_file='test.csv',target_silva_accession='GAFF01033989.4391.6188',sequence_silva_path=fungi_taxpath[0],
    primer_F=F,primer_R=R,kmer_range=(9,14),database_dir="silva_db_compiled/")

