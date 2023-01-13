# -*- coding: utf-8 -*-
# @Author: charl
# @Date:   2018-05-21 20:28:04
# @Last Modified by:   charl
# @Last Modified time: 2019-03-20 15:55:42

import pandas as pd 
import numpy as np 
import itertools
import constants


class Codons(object):

	def __init__(self):
		#load and parse naturally occuring codon file
		self.natural_codon_file = pd.ExcelFile(constants.natural_codons)
		self.natural_codons = self.natural_codon_file.parse(
							sheet_name = 'codons', index_col = 'Codon'
							)
		self.degenerate_nucleotides = constants.degenerate_nucleotides
		self.amino_acids = constants.amino_acids
		#create all potential combinations of naturally occuring codons
		self.all_codons = self.construct_codons()
		#check if excel file exists already containing codon translations
		#if not, tranlate and create this file for future use
		try:
			translated_file = pd.ExcelFile(constants.translated_codons)
			self.translated_codons = translated_file.parse(
				sheet_name = 'translated_codons', index_col = 0
				)
		except:
			self.translated_codons = self.translate_codons()

	def construct_codons(self):
		#creates all possible codons using degenerate nucleotides
		codons_iterator = itertools.product(
			self.degenerate_nucleotides.keys(), repeat = 3
			)
		codons = [''.join(i) for i in codons_iterator]
		return codons

	def translate_codons(self):
		#determine amino acid complement for each degenerate codon
		codon_frame = pd.DataFrame(
			0, index = self.all_codons, columns = self.amino_acids
			)
		for degenerate_codon in codon_frame.index:
			unpacked_codon = self.unpack_degenerates(degenerate_codon)
			for codon in unpacked_codon:
				amino_acid = self.natural_codons.loc[codon, 'Amino acid']
				if not codon_frame.loc[degenerate_codon, amino_acid]:
					codon_frame.loc[degenerate_codon, amino_acid] = 1
		writer = pd.ExcelWriter('translated_codons.xlsx')
		codon_frame.to_excel(writer, sheet_name = 'translated_codons')
		writer.save()
		return codon_frame

	def unpack_degenerates(self, codon):
		#create all possible codons from natural nucleotides
		codon_holder = ['']

		for nucleotide in codon:
			possible_nucleotides = self.degenerate_nucleotides[nucleotide]
			codon_holder *= len(possible_nucleotides)
			for index, possible_nucleotide in enumerate(possible_nucleotides):
				start = (index*len(codon_holder))//len(possible_nucleotides)
				stop = (index + 1)*len(codon_holder)//len(possible_nucleotides)
				for position in range(start, stop):
					codon_holder[position] += possible_nucleotide
		return codon_holder
