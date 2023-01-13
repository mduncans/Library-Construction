# -*- coding: utf-8 -*-
# @Author: charl
# @Date:   2018-02-03 14:42:51
# @Last Modified by:   charl
# @Last Modified time: 2019-03-20 15:52:27
import os
from tkinter import Tk
from tkinter.filedialog import askopenfilename

#UPDATE library_generator DIRECTORY HERE:
construction_dir = '/Users/mattsmith/Dropbox (University of Michigan)/from_box/My_Stuff/Python_Codes/mds_pkgs/Library Construction'

cwd = os.getcwd()
os.chdir(construction_dir)

amino_acids = [
                'G',
                'A',
                'L',
                'I',
                'V',
                'M',
                'H',
                'P',
                'F',
                'Y',
                'W',
                'S',
                'T',
                'C',
                'R',
                'K',
                'N',
                'Q',
                'D',
                'E',
                'STOP'
                ]
#this dictionary is used by codon_construction.py to create every possible
#combination of degenerate codons
degenerate_nucleotides = {
						'A':['A'],
						'C':['C'],
						'T':['T'],
						'G':['G'],
						'N':['A','C','T','G'],
						'B':['C','T','G'],
						'D':['A','T','G'],
						'H':['A','C','T'],
						'V':['A','C','G'],
						'K':['G','T'],
						'M':['A','C'],
						'R':['A','G'],
						'S':['C','G'],
						'W':['A','T'],
						'Y':['C','T']
					}
#this dictionary is used to identify the CDR regions in the framework
#based on kabat numbering scheme
# can customize cdr regions using alternate numbering schemes but update
#numbers to reflect kabat numbers
cdr_positions = {
				'HCDR1': ('26','35B'),
				'HCDR2': ('50','65'),
				'HCDR3': ('93','102'),
				'LCDR1': ('24','34'),
				'LCDR2': ('50','56'),
				'LCDR3': ('89','97')
}
#this dictionary is used to generate primers for libary construction
#it is specific for the pd208 Fab framework
flanking_regions = {
					'hcdr_upstream': 'GGTGGTACTGCTACCGCTGG',
					'hcdr_downstream': 'CCAGATGTAGATTTAGAGGAAGGAGCC',
					'lcdr_upstream': 'TTGGATAAAAGAGAAGCTGCTAGCG',
					'lcdr_downstream': 'GGAAGATGAAGACAGATGGTGCAGC'
					}

compliment = {
			'A':'T',
			'T':'A',
			'C':'G',
			'G':'C'
			}

#these excel files contain information necessary for the construction of 
#any library based on this scheme
#translation of 64 naturally occuring codons
natural_codons = 'codon_table.xlsx'
#translation of degenerate codons
translated_codons = 'translated_codons.xlsx'
#variable region frameworks aligned to kabat numbering
framework_file = 'framework_alignment.xlsx'
#frequency of occurence of each amino acid at each variable region site
root = Tk()
root.withdraw()
diversity_file = askopenfilename()
diversity_file = diversity_file.split('/')[-1]
root.update() 
#sequence of fab variable region gene blocks
nucleotide_file = 'fab_geneblocks.xlsx' 


#Version 2 Rules
diversification = {
	1: (['Y', 'D', 'G', 'S', 'C', 'V', 'F', 'A'], [('Y', 2), ('A', 2), ('D', 2)], ['Y', 'D', 'A']),
	2: (['Y', 'D', 'G', 'S', 'C', 'V', 'F', 'A'], [('Y', 2), ('A', 2)], ['Y', 'A']),
	3: (['Y', 'D', 'G', 'S', 'C', 'V', 'F', 'A'], [('Y', 2)], ['Y']),
	4: (['Y', 'D', 'G', 'S', 'C', 'V', 'F', 'A'], [('A', 2)], ['A']),
	5: (['Y', 'D', 'G', 'S', 'C', 'V', 'F', 'A'], [('D', 2)], ['D']),
	6: (['N', 'T'], [('Y', 2), ('A', 2), ('D', 2)], ['Y', 'A', 'D']),
	7: (['N', 'T'], [('Y', 2), ('A', 2)], ['Y', 'A']),
	8: (['N', 'T'], [('Y', 2)], ['Y']),
	9: (['N', 'T'], [('A', 2)], ['A']),
	10: (['N', 'T'], [('D', 2)], ['D']),
	11: (['P', 'L'], [('A', 2)], ['A']),
	12: (['M'], [('A', 2)], ['A']),
	13: (['Q', 'E'], [('A', 2), ('E', 2)], ['A', 'E']),
	14: (['Q', 'E'], [('A', 2)], ['A']),
	15: (['Q', 'E'], [('E', 2)], ['E']),
	16: (['I'], [('Y', 2)], ['Y']),
	17: (['I'], [('A', 2), ('D', 2)], ['A', 'D']),
	18: (['I'], [('A', 2)], ['A']),
	19: (['I'], [('D', 2)], ['D']),
}


#Natural Diversity No Bias
diversification_no_bias = {
	1: (
		['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
		'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'],
		[('A', 0)], 
		[]
		)
}

homologs = {
	'A': 'KCT', 
	'C': 'TSC', 
	'D': 'GAM',
	'E': 'GAM',
	'F': 'TWC',
	'G': 'GST',
	'H': 'MAC',
	'I': 'RTT',
	'K': 'ARG',
	'L': 'MTC',
	'M': 'MTG',
	'N': 'RAC',
	'P': 'SCA',
	'Q': 'SAA',
	'R': 'ARG',
	'S': 'KCC',
	'T': 'ASC',
	'V': 'RTT',
	'W': 'TKG',
	'Y': 'TWC'
}

AA_homologs = {
	'A': 'S',
	'C': 'S',
	'D': 'E',
	'E': 'D',
	'F': 'Y',
	'G': 'A',
	'H': 'N',
	'I': 'V',
	'K': 'R',
	'L': 'I',
	'M': 'L',
	'N': 'D',
	'P': 'A',
	'Q': 'E',
	'R': 'K',
	'S': 'A',
	'T': 'S',
	'V': 'I',
	'W': 'L',
	'Y': 'F',
}

def soft(codon):
	soft_conversion = {
		'A': '(70101010)',
		'G': '(10107010)',
		'C': '(10701010)',
		'T': '(10101070)'
	}
	soft = ''
	for nt in codon:
		soft += soft_conversion[nt]

	return soft

soft_codons = {
	'A':	'677',
	'C':	'868',
	'D':	'657',
	'E':	'656',
	'F':	'888',
	'G':	'668',
	'H':	'758',
	'I':	'588',
	'K':	'556',
	'L':	'787',
	'M':	'586',
	'N':	'557',
	'P':	'777',
	'Q':	'756',
	'R':	'768',
	'S':	'878',
	'T':	'577',
	'V':	'688',
	'W':	'866',
	'Y':	'857',
}
soft_conversion = {
	'5': '(70101010)',
	'6': '(10107010)',
	'7': '(10701010)',
	'8': '(10101070)'
}
soft_codons_numbers = {
	k: ''.join([soft_conversion[n] for n in v]) for k, v in soft_codons.items()
}

os.chdir(cwd)