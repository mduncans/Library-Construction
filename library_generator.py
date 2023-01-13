# -*- coding: utf-8 -*-
# @Author: charl
# @Date:   2018-05-22 11:07:06
# @Last Modified by:   charl
# @Last Modified time: 2019-04-17 15:49:52

import pandas as pd 
import numpy as np 
import itertools
import constants
import codon_construction
import copy
import os
from openpyxl import load_workbook
from tabulate import tabulate

class LibraryGenerator(object):

	def __init__(self, codon_list, framework, library_regions, excluded_residues, size_limit, residue_limit, cysteine_allowance, Homolog_Library = False, no_bias = False):

		self.codon_list = codon_list
		self.framework = framework
		self.excluded_residues = excluded_residues
		self.size_limit = size_limit
		self.residue_limit = residue_limit
		self.cysteine_allowance = cysteine_allowance
		self.Homolog_Library = Homolog_Library

		if no_bias:
			self.diversification = constants.diversification_no_bias
		else:
			self.diversification = constants.diversification

		self.natural_diversity = self.load_diversity()
		self.heavy_frame, self.light_frame = self.load_framework(framework)
		self.variable_heavy, self.variable_light = self.trim_frameworks(library_regions, framework)	
		self.framework_codon_table(framework)	
		print('Frameworks Trimmed')
		
		self.residue_catalog = self.identify_diversification()
		print('Diversification identified')
		self.residue_catalog = self.generate_library()
		print('Library Generated')
		self.best_codon()
		print('Best codons chosen')
		self.heavy_out, self.light_out = self.format_ouput(framework, library_regions, Homolog_Library, no_bias)
		print('Output created')

	def load_diversity(self):
		#load natural diversity profile for each site in each variable region
		#into a pandas dataframe
		div_file = pd.ExcelFile(constants.diversity_file)
		h = div_file.parse(sheetname = 'heavy', index_col = 'Amino acids')
		l = div_file.parse(sheetname = 'light', index_col = 'Amino acids')
		
		h.columns = h.columns.astype(str)
		l.columns = l.columns.astype(str)

		return {'heavy': h, 'light': l}

	def load_framework(self, framework):
		#load the framework aminio acid sequence for all frameworks in database
		#into pandas dataframe
		frame_file = pd.ExcelFile(constants.framework_file)
		heavy_frame = frame_file.parse(sheetname = 'heavy', index_col = 'Name')
		light_frame = frame_file.parse(sheetname = 'light', index_col = 'Name')
		#remove all frameworks from dataframe except specified framework
		heavy_frame.drop(heavy_frame[heavy_frame.index != framework].index, inplace = True)
		light_frame.drop(light_frame[light_frame.index != framework].index, inplace = True)

		heavy_frame.columns = heavy_frame.columns.astype(str)
		light_frame.columns = light_frame.columns.astype(str)

		return (heavy_frame, light_frame)

	def framework_codon_table(self, framework):
		nucleotide_file = pd.ExcelFile(constants.nucleotide_file)
		nucleotide_table = nucleotide_file.parse(sheetname = 'frameworks', index_col = 'Framework')

		self.heavy_framework_codons = pd.DataFrame(columns = ['Residue', 'Codon'])
		j = 0
		for i, pos in enumerate(self.heavy_frame.T.index):
			if self.heavy_frame.T.loc[pos][0] != '-':
				self.heavy_framework_codons.loc[pos] = [self.heavy_frame.T.loc[pos][0], nucleotide_table.loc[framework, 'VH sequence'][j*3 : 3*j + 3]]
				j += 1
			else:
				self.heavy_framework_codons.loc[pos] = ['-', '-']

		self.light_framework_codons = pd.DataFrame(columns = ['Residue', 'Codon'])
		j = 0
		for i, pos in enumerate(self.light_frame.T.index):
			if self.light_frame.T.loc[pos][0] != '-':
				self.light_framework_codons.loc[pos] = [self.light_frame.T.loc[pos][0], nucleotide_table.loc[framework, 'VL sequence'][j*3 : 3*j + 3]]
				j += 1
			else:
				self.light_framework_codons.loc[pos] = ['-', '-']

	def trim_frameworks(self, library_regions, framework):
		#isolate potential sites for library variability using CDR definitions
		variable_heavy = pd.Series()
		variable_light = pd.Series()
		#select variable sites in heavy chain
		for region in library_regions:
			if region[0] == 'H':
				variable_heavy = variable_heavy.append(self.heavy_frame.loc[framework, constants.cdr_positions[region][0]:constants.cdr_positions[region][1]])
		#select variable sites in light chain
			else:
				variable_light = variable_light.append(self.light_frame.loc[framework, constants.cdr_positions[region][0]:constants.cdr_positions[region][1]])

		return (variable_heavy, variable_light)

	def identify_diversification(self):
		#Identify all sites amenable to diversification based on DY design paradigm
		#could make residues and cutoffs a dictionary to alter design strategy

		residue_catalog = pd.DataFrame(columns = ['Chain', 'Wild type', 'Tier', 'Initial Diversity', 'Variability', 'Variable residues', 'Codon', 'Diversity Coverage', 'Homolog Codon', 'Soft Codon'])
		#Iterate through potential variable sites for heavy chain to identify teir 1, 2, and 3 sites. 
		for residue in self.variable_heavy.index:
			#Calculate natural diversity of the site in natural antibodies
			variability = self.natural_diversity['heavy'][residue].max() #highly conserved residues have variability near 1
			if self.variable_heavy[residue] != '-':
				init_diversity = self.natural_diversity['heavy'][residue][self.variable_heavy[residue]]
			else:
				init_diversity = ''

			for tier, requirements in self.diversification.items():
				if self.variable_heavy[residue] in requirements[0]:
					diversity_count = 0
					for nd in requirements[1]:
						if self.natural_diversity['heavy'].loc[nd[0], residue] >= nd[1]:
							diversity_count	+=1
						else:
							break
					if diversity_count == len(requirements[1]):
						residue_catalog.loc[residue] = ['heavy', self.variable_heavy[residue], tier, init_diversity, variability, '','','', '', '']
						break

		print(residue_catalog)
		#Iterate through potential variable sites in light chain
		for residue in self.variable_light.index:
			#Calculate natural diversity of site in natural antibodies
			variability = self.natural_diversity['light'][residue].max()
			if self.variable_light[residue] != '-':
				init_diversity = self.natural_diversity['light'][residue][self.variable_light[residue]]
			else:
				init_diversity = ''
			
			for tier, requirements in self.diversification.items():
				if self.variable_light[residue] in requirements[0]:
					diversity_count = 0
					for nd in requirements[1]:
						if self.natural_diversity['light'].loc[nd[0], residue] >= nd[1]:
							diversity_count	+=1
						else:
							break
					if diversity_count == len(requirements[1]):
						residue_catalog.loc[residue] = ['light', self.variable_light[residue], tier, init_diversity, variability, '','','', '', '']

		#sort residue catalog to prioritize variable sites with the most natural diversity
		print(residue_catalog)
		residue_catalog.sort_values('Variability', axis = 0, inplace = True)
		
		#Drop residues with wt > 50 natural diversity.
		if self.Homolog_Library == False:
			for residue in residue_catalog.index:
				if self.natural_diversity[residue_catalog.loc[residue, 'Chain']][residue][residue_catalog.loc[residue, 'Wild type']] > 50:
					residue_catalog = residue_catalog.drop(residue, axis = 0)
			print(tabulate(residue_catalog, headers = residue_catalog.columns, tablefmt = 'github'))
		
		#Unless this is homolog library
		else:
			for residue in residue_catalog.index:
				if self.natural_diversity[residue_catalog.loc[residue, 'Chain']][residue][residue_catalog.loc[residue, 'Wild type']] >= 80:
					residue_catalog = residue_catalog.drop(residue, axis = 0)

				elif self.natural_diversity[residue_catalog.loc[residue, 'Chain']][residue][residue_catalog.loc[residue, 'Wild type']] >= 50:
					if self.natural_diversity[residue_catalog.loc[residue, 'Chain']][residue][constants.AA_homologs[residue_catalog.loc[residue, 'Wild type']]] < 5:
						residue_catalog = residue_catalog.drop(residue, axis = 0)

			print(tabulate(residue_catalog, headers = residue_catalog.columns, tablefmt = 'github'))

		return residue_catalog

	def generate_library(self):
		#a library with no variable sites has size 1, this counter keeps track of size
		library_size = 1

		#iterate through potential variable sites to choose codons
		for residue in self.residue_catalog.index:
			#Permit excluded residue if it is the wildtype residue
			wild_type = self.residue_catalog.loc[residue, 'Wild type']
			if wild_type in self.excluded_residues:
				temp_excluded = self.excluded_residues.copy()
				temp_excluded.remove(wild_type)
			else:
				temp_excluded = self.excluded_residues.copy()

			#choose residues, codon. Report sampled diversity
			new_residue, new_codon, new_diversity = self.select_codon(
				residue, wild_type, self.residue_catalog.loc[residue, 'Chain'], self.residue_catalog.loc[residue, 'Tier'], temp_excluded, self.residue_limit, 0
				)
			
			#record results of codon selection in residue catalog
			self.residue_catalog.loc[residue, 'Variable residues'] = ''.join(new_residue)
			self.residue_catalog.loc[residue, 'Codon'] = ','.join(new_codon)
			self.residue_catalog.loc[residue, 'Diversity'] = new_diversity
			self.residue_catalog.loc[residue, 'Homolog Codon'] = constants.homologs[wild_type]
			
			if self.residue_catalog.loc[residue, 'Chain'] == 'heavy':
				self.residue_catalog.loc[residue, 'Soft Codon'] = constants.soft_codons_numbers[wild_type]
			elif self.residue_catalog.loc[residue, 'Chain'] == 'light':
				self.residue_catalog.loc[residue, 'Soft Codon'] = constants.soft_codons_numbers[wild_type]

			#adjust library size based on number of residues 
			library_size *= len(new_residue)
			#exit residue selection process when size limit is exceeded
			if library_size > self.size_limit:
				break
		
		#drop unvaried residues from the residue catalog
		self.residue_catalog.dropna(axis = 0, inplace = True)
		
		#call function to remove cysteine residues such that total libary
		#variants with cysteine is less than residue limit
		self.fix_cysteine()
		'''
		self.fix_asparagine(
			excluded_residues, residue_limit
			)
		'''
		return self.residue_catalog

	def select_codon(self, position, wild_type, chain, tier, excluded_residues, residue_limit, cysteine_flag):
		#make deep copy of codon list so that it can be altered
		available_codons = self.codon_list.translated_codons.copy(deep = True)

		#establish amino acid residues that must be included in selected codon
		fixed_residues = self.diversification[tier][2].copy()
		if wild_type not in fixed_residues:
			fixed_residues.append(wild_type)

		#remove tyrosine from fixed residues if we are fixing cysteine
		if wild_type == 'G' and cysteine_flag == 1 and 'Y' in fixed_residues:
			fixed_residues.remove('Y')
		
		#remove all codons that code for an excluded residue
		if excluded_residues == None:
			print(position)
		
		for residue in excluded_residues:
			available_codons.drop(available_codons[available_codons[residue] == 1].index, inplace = True)

		#remove codons that do not code for a fixed residue
		for residue in fixed_residues:
			available_codons.drop(available_codons[available_codons[residue] == 0].index, inplace = True)
			
		#remove codons that code for more residues than the residue limit
		available_codons.drop(available_codons[available_codons.sum(axis = 1) > residue_limit].index, inplace = True)

		#convert binary value to natural diversity for each residue at the variable site
		for residue in available_codons.columns:
			available_codons[residue] *= (self.natural_diversity[chain].loc[residue, position])

		#calculate diversity sampled and sort
		available_codons['Diversity Sum'] = available_codons.sum(axis = 1)
		available_codons.sort_values('Diversity Sum', axis = 0, inplace = True, ascending = False)
		
		#select codon that provides most natural diversity
		try:
			codon = available_codons.index[0]
		except:
			print(position)
			print(excluded_residues)

		residues = self.codon_list.translated_codons.loc[codon,:][self.codon_list.translated_codons.loc[codon,:] == 1].index.tolist()
		codons = [codon]
		#select additional potential codons that code for same residues
		for replicate in available_codons.index[1:]:
			if available_codons.loc[replicate,:].equals(available_codons.loc[codon,:]):
				codons.append(replicate)
		diversity = available_codons.loc[codon, 'Diversity Sum']


		return (residues, codons, diversity)

	def fix_cysteine(self):
		#calculate cysteine frequency in unmodified library
		cysteine_frequency = self.calculate_cysfreq()
		#perform cysteine reduction if cysteine frequency exceeds limit
		if cysteine_frequency <= self.cysteine_allowance:
			return
		#iterate through twice to avoid removing cysteine from G wild-type
		for iteration in range(0,2,1):
			if cysteine_frequency <= self.cysteine_allowance:
				break
			#iterate through variable residues starting with least naturally diverse site
			for residue in self.residue_catalog.index[::-1]:
				wild_type = self.residue_catalog.loc[residue, 'Wild type']
				if wild_type == 'G' and iteration == 0:
					continue
				#check for cysteine in variable residues
				if 'C' in self.residue_catalog.loc[residue, 'Variable residues']:
					#don't alter residue if cysteine is wild-type
					if 'C' == wild_type:
						continue
					else:
						#add cysteine to excluded residues list
						temp_excluded = self.excluded_residues + ['C']
						#select a new codon based on cysteine exclusion
						new_residue, new_codon, new_diversity = self.select_codon(
							residue, wild_type, self.residue_catalog.loc[residue, 'Chain'], self.residue_catalog.loc[residue, 'Tier'], temp_excluded, self.residue_limit, 1
							)
						self.residue_catalog.loc[residue, 'Variable residues'] = ''.join(new_residue)
						self.residue_catalog.loc[residue, 'Codon'] = ','.join(new_codon)
						self.residue_catalog.loc[residue, 'Diversity'] = new_diversity
						self.residue_catalog.loc[residue, 'Cysteine removed?'] = 'Yes'

				cysteine_frequency = self.calculate_cysfreq()
				if cysteine_frequency <= self.cysteine_allowance:
					break

		return

	def calculate_cysfreq(self):
		#calculate the number of library variants with at least one cysteine
		site_probabilities = []
		for amino_acids in self.residue_catalog['Variable residues']:
			if 'C' in amino_acids:
				site_probabilities.append((len(amino_acids)-1.0)/len(amino_acids))
			else:
				site_probabilities.append(1)

		cysteine_frequency = 1 - np.prod(site_probabilities)

		return cysteine_frequency
	
	'''
	def fix_asparagine(self, excluded_residues, residue_limit):
		#identify and modify sites with potential for N-linked glycosylation
		for residue in self.residue_catalog.index:
			if self.residue_catalog.loc[residue, 'Chain'] == 'heavy':
				framework = self.heavy_frame
			else:
				framework = self.light_frame
			wild_type = self.residue_catalog.loc[residue, 'Wild type']
			if wild_type in ['N', 'S', 'T']:!!!!!!!!!!!!!!!!
				continue
			if 'N' in self.residue_catalog.loc[residue, 'Variable residues']:
				residue_index = framework.columns.get_loc(residue)
				residue_index += 2
				if framework.iloc[0, residue_index] in ['S', 'T']:
					if framework.columns[residue_index] not in self.residue_catalog.index:
						print self.framework
						temp_excluded = excluded_residues + ['N']
						new_residue, new_codon, new_diversity = self.select_codon(
							residue, wild_type, 
							self.residue_catalog.loc[residue, 'Chain'],
							self.residue_catalog.loc[residue, 'Tier'], 
							temp_excluded, residue_limit, 1
							)
						print residue
						print new_residue
						print new_codon
						print new_diversity
						self.residue_catalog.loc[residue, 'Variable residues'] = ''.join(new_residue)
						self.residue_catalog.loc[residue, 'Codon'] = ','.join(new_codon)
						self.residue_catalog.loc[residue, 'Diversity'] = new_diversity
						self.residue_catalog.loc[residue, 'Asp fix?'] = 'Yes'
			for amino_acid in ['S', 'T']:
				if amino_acid in self.residue_catalog.loc[residue, 'Variable residues']:
					residue_index = framework.columns.get_loc(residue)
					residue_index -= 2
					if framework.iloc[0, residue_index] == 'N':
						if framework.columns[residue_index] not in self.residue_catalog.index:
							print self.framework
							temp_excluded = excluded_residues + ['S', 'T']
							new_residue, new_codon, new_diversity = self.select_codon(
								residue, wild_type, 
								self.residue_catalog.loc[residue, 'Chain'],
								self.residue_catalog.loc[residue, 'Tier'], 
								temp_excluded, residue_limit, 1
								)
							self.residue_catalog.loc[residue, 'Variable residues'] = ''.join(new_residue)
							self.residue_catalog.loc[residue, 'Codon'] = ','.join(new_codon)
							self.residue_catalog.loc[residue, 'Diversity'] = new_diversity
							self.residue_catalog.loc[residue, 'Asp fix?'] = 'Yes'
		return
	'''
	def best_codon(self):

		for residue in self.residue_catalog.index:
			if self.residue_catalog.loc[residue, 'Codon']:
				potential_codons = self.residue_catalog.loc[residue, 'Codon'].split(',')
				codon_series = pd.Series(0, index = potential_codons)
				for codon in potential_codons:
					unpacked_codons = self.codon_list.unpack_degenerates(codon)
					total_usage = 0.0
					for single in unpacked_codons:
						total_usage += ((self.codon_list.natural_codons.loc[single, 'Human usage'] + self.codon_list.natural_codons.loc[single, 'Yeast usage'])/2.0)
					total_usage /= len(unpacked_codons)
					codon_series.loc[codon] = total_usage
				codon_series.sort_values(ascending = False, inplace = True)
				self.residue_catalog.loc[residue, 'Best codon'] = codon_series.index[0]

	def format_ouput(self, framework, library_regions, Homolog_Library, no_bias):

		heavy_out = pd.DataFrame(
			'', index = self.variable_heavy.index, 
			columns = ['Wild type', 'Tier', 'Initial Diversity', 'Variable residues', 'Codon', 
			'Best codon', 'Variability', 'Diversity', 'Cysteine removed?', 'Homolog Codon', 'Soft Codon'] + 
			self.natural_diversity['heavy'].index.tolist())
		
		light_out = pd.DataFrame(
			'', index = self.variable_light.index, 
			columns = ['Wild type', 'Tier', 'Initial Diversity', 'Variable residues', 'Codon', 
			'Best codon', 'Variability', 'Diversity', 'Cysteine removed?', 'Homolog Codon', 'Soft Codon'] + 
			self.natural_diversity['light'].index.tolist())

		for residue in self.variable_heavy.index:
			if (residue in self.residue_catalog.index and self.residue_catalog.loc[residue, 'Chain'] == 'heavy'):
				columns = ['Wild type', 'Tier', 'Initial Diversity', 'Variable residues', 'Codon', 'Best codon', 'Variability', 'Diversity', 'Cysteine removed?', 'Homolog Codon', 'Soft Codon']
				heavy_out.loc[residue, columns] = self.residue_catalog.loc[residue, columns]
			else:
				heavy_out.loc[residue, 'Homolog Codon'] = self.heavy_framework_codons.loc[residue][1]
				heavy_out.loc[residue, 'Soft Codon'] = self.heavy_framework_codons.loc[residue][1]
				
			columns = self.natural_diversity['heavy'].index
			heavy_out.loc[residue, columns] = self.natural_diversity['heavy'].loc[columns, residue]
			heavy_out.loc[residue, 'Wild type'] = self.variable_heavy[residue]
			
		
		for residue in self.variable_light.index:
			if (residue in self.residue_catalog.index and self.residue_catalog.loc[residue, 'Chain'] == 'light'):
				columns = ['Wild type', 'Tier', 'Initial Diversity', 'Variable residues', 'Codon', 'Best codon', 'Variability', 'Diversity', 'Cysteine removed?', 'Homolog Codon', 'Soft Codon']
				light_out.loc[residue, columns] = self.residue_catalog.loc[residue, columns]
			else:
				light_out.loc[residue, 'Homolog Codon'] = self.light_framework_codons.loc[residue][1]
				light_out.loc[residue, 'Soft Codon'] = self.light_framework_codons.loc[residue][1]
			
			columns = self.natural_diversity['light'].index
			light_out.loc[residue, columns] = self.natural_diversity['light'].loc[columns, residue]
			light_out.loc[residue, 'Wild type'] = self.variable_light[residue]
		
		if not Homolog_Library:
			regions = ''
			for r in library_regions:
				regions += r + '_'
			regions = regions[0:-1]

			if no_bias:
				regions += '_No_Bias'
		else:
			regions = 'Homolog'


		mydir = os.path.dirname(os.path.realpath(__file__))
		writer = pd.ExcelWriter(os.path.join(mydir, 'V3 Libraries/' + framework + '_' + regions + '.xlsx'))
		heavy_out.to_excel(writer, sheet_name = 'heavy')
		light_out.to_excel(writer, sheet_name = 'light')
		writer.save()

		return [heavy_out, light_out]

class PrimerGenerator(object):

	def __init__(self, codons, framework, library, library_regions):

		self.codons = codons
		self.library = library
		print(framework)
		self.aligned_heavy, self.aligned_light = self.align_nucleotides(framework)
		print('nucleotides aligned')
		self.variable_primers = self.diversity_primers(framework, library_regions)
		print('Primers created')
		self.terminal_primers = self.terminal_selection(library_regions)
		print('Final Primers selected')
		self.primers_out(framework, library_regions)
		print('Primers output')

	def align_nucleotides(self, framework):
		nucleotide_file = pd.ExcelFile(constants.nucleotide_file)
		nucleotide_table = nucleotide_file.parse(sheetname = 'frameworks', index_col = 'Framework')

		aligned_heavy = pd.DataFrame(index = self.library.heavy_frame.columns, columns = ['Amino acid', 'Codon'])
		aligned_light = pd.DataFrame(index = self.library.light_frame.columns, columns = ['Amino acid', 'Codon'])

		alignments = [aligned_heavy, aligned_light]
		aa_sequences = [self.library.heavy_frame, self.library.light_frame]
		domains = ['VH sequence', 'VL sequence']

		for alignment, aa_sequence, domain in zip(alignments, aa_sequences, domains):

			alignment['Amino acid'] = aa_sequence.loc[framework, :]
			start = 0
			end = 3
			for position in alignment.index:
				if alignment.loc[position, 'Amino acid'] != '-':
					alignment.loc[position, 'Codon'] = (nucleotide_table.loc[framework, domain][start:end])
					start += 3
					end += 3

		return alignments

	def diversity_primers(self, framework, library_regions):
		primers = []

		for region in library_regions:
			if region[0] == 'H':
				lib_residues = self.library.heavy_out
				alignment = self.aligned_heavy.drop(self.aligned_heavy[self.aligned_heavy['Amino acid'] == '-'].index)
			else:
				lib_residues = self.library.light_out
				alignment = self.aligned_light.drop(self.aligned_light[self.aligned_light['Amino acid'] == '-'].index)

			border = constants.cdr_positions[region]
			potential_bubble = lib_residues.loc[border[0]:border[1],:]
			true_bubble = self.find_bubble(potential_bubble)
			
			if true_bubble is not None:
				bubble_sequence = self.generate_bubble(true_bubble, region)
			else:	
				primers.append(['','',''])
				continue
			
			prebubble_alignment = alignment.loc[:true_bubble.index[0], :]
			postbubble_alignment = alignment.loc[true_bubble.index[-1]:, :]

			prebubble_sequence = ''.join(prebubble_alignment.iloc[-11:-1, 1].tolist())
			postbubble_sequence = ''.join(postbubble_alignment.iloc[1:11, 1].tolist())

			degenerate_primer = (prebubble_sequence + bubble_sequence + postbubble_sequence)
			reverse_overlap = self.reverse_compliment(prebubble_sequence)
			primers.append([degenerate_primer, reverse_overlap, postbubble_sequence])

		return primers

	def find_bubble(self, potential_bubble):
		
		start_bubble = ''
		end_bubble = ''
		for residue in potential_bubble.index:
			if potential_bubble.loc[residue, 'Variable residues']:
				start_bubble = residue
				break

		for residue in potential_bubble.index[::-1]:
			if potential_bubble.loc[residue, 'Variable residues']:
				end_bubble = residue
				break

		if not start_bubble:
			return None
		else:
			return potential_bubble.loc[start_bubble:end_bubble, :]

	def generate_bubble(self, true_bubble, region):

		bubble_sequence = ''
		if region[0] == 'H':
			alignment = self.aligned_heavy
		else:
			alignment = self.aligned_light

		for residue in true_bubble.index:
			if true_bubble.loc[residue, 'Wild type'] == '-':
				continue
			elif true_bubble.loc[residue, 'Variable residues']:
				bubble_sequence += true_bubble.loc[residue, 'Best codon']
			else:
				block_codon = alignment.loc[residue, 'Codon']
				best_mismatch = self.generate_mismatch(block_codon)
				bubble_sequence += best_mismatch

		return bubble_sequence

	def generate_mismatch(self, block_codon):

		encoded_aa = self.codons.natural_codons.loc[block_codon, 'Amino acid']
		if encoded_aa in ['W', 'M']:
			return block_codon
		viable_choices = self.codons.natural_codons[self.codons.natural_codons.loc[:,'Amino acid'] == encoded_aa]

		viable_choices.sort_values('Average', inplace = True, ascending = False)
		viable_choices.drop(viable_choices.index[-1], inplace = True)

		for choice in viable_choices.index:
			viable_choices.loc[choice, 'Mismatch'] = 0
			for index, nucleotide in enumerate(block_codon):
				if nucleotide != choice[index]:
					viable_choices.loc[choice, 'Mismatch'] += 1
		
		viable_choices.sort_values('Mismatch', inplace = True, ascending = False, kind = 'mergesort')

		best_mismatch = viable_choices.index[0]
		return best_mismatch

	def reverse_compliment(self, sequence):

		reverse_compliment = ''

		for residue in sequence[::-1]:
			reverse_compliment += constants.compliment[residue]

		return reverse_compliment

	def terminal_selection(self, library_regions):

		terminal_primers = []
		heavy_flag = False
		light_flag = False

		for region in library_regions:
			if region[0] == 'H':
				heavy_flag = True
			else:
				light_flag = True

		if heavy_flag == True and light_flag == True:
			terminal_primers.append(constants.flanking_regions['lcdr_upstream'])
			terminal_primers.append(constants.flanking_regions['hcdr_downstream'])
		elif heavy_flag == True and light_flag == False:
			terminal_primers.append(constants.flanking_regions['hcdr_upstream'])
			terminal_primers.append(constants.flanking_regions['hcdr_downstream'])
		elif heavy_flag == False and light_flag == True:
			terminal_primers.append(constants.flanking_regions['lcdr_upstream'])
			terminal_primers.append(constants.flanking_regions['lcdr_downstream'])

		return terminal_primers

	def primers_out(self, framework, library_regions):

		title_catalog = ['Terminal forward', 'Terminal reverse']
		primer_catalog = [self.terminal_primers[0], self.terminal_primers[1]]


		for primer, region in zip(self.variable_primers, library_regions):

			title_catalog.extend([region + ' degenerate', region + ' reverse overlap', region + ' forward helper'])
			primer_catalog.extend([primer[0], primer[1], primer[2]])

		output = pd.Series(primer_catalog, index = title_catalog)
		mydir = os.path.dirname(os.path.realpath(__file__))
		path = os.path.join(mydir, 'Nanobody/' + framework + '.xlsx')
		book = load_workbook(path)
		writer = pd.ExcelWriter(path, engine = 'openpyxl')
		writer.book = book
		writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
		output.to_excel(writer, sheet_name = 'primers')
		writer.save()


def main(framework, library_regions, excluded_residues, size_limit, residue_limit, cysteine_allowance, Homolog_Library = False, no_bias = False):
	cwd = os.getcwd()
	os.chdir(constants.construction_dir)
	codons = codon_construction.Codons()
	
	for scaffold in framework:
		print('designing {} library now'.format(scaffold))
		library = LibraryGenerator(codons, scaffold, library_regions, excluded_residues, size_limit, residue_limit, cysteine_allowance, Homolog_Library, no_bias)

		#print(f'Designing {scaffold} Primers now.')
		#PrimerGenerator(codons, scaffold, library, library_regions)
	os.chdir(cwd)
	return library

def order_sheet(scaffolds, seed):

	mydir = os.path.dirname(os.path.realpath(__file__))
	order_frame = pd.DataFrame(0, index = [], columns = ['Sequence', 'Amt', 'Prep', 'Description'])
	for scaffold in scaffolds:
		path = os.path.join(mydir, 'library_designs_2/' + scaffold + '.xlsx')
		excel_file = pd.ExcelFile(path)
		excel_table = excel_file.parse(sheetname = 'primers')
		for primer in excel_table.index:
			if 'Terminal' in primer:
				continue
			if not pd.isnull(excel_table.loc[primer, 0]):
				name = 'AD pr' + str(seed)
				seed += 1
				if len(excel_table.loc[primer, 0]) >= 60:
					size = '4nmU'
				else:
					size = '25nm'
				order_frame.at[name, :] = [excel_table.loc[primer, 0], size, 'STD', scaffold + ' ' + primer]
	path = os.path.join(mydir, 'Nanobody/' + 'order.xlsx')

	writer = pd.ExcelWriter(path)
	order_frame.to_excel(writer)
	writer.save()

def find_codon(aa_list): 
	codons = pd.read_excel('translated_codons.xlsx', index_col = 0)
	codons.insert(21, 'Total', [0 for _ in range(len(codons.index))])
	
	for c in codons.index: 
		codons.loc[c, 'Total'] = sum(codons.loc[c, codons.columns != 'Total'])

	pos_codons = [] 
	for c in codons.index: 
		if codons.loc[c, 'Total'] == len(aa_list): 
			count = 0 
			for aa in aa_list: 
				if codons.loc[c, aa] == 1: 
					count += 1 
			if count == len(aa_list): 
				pos_codons.append(c) 
	return pos_codons

if __name__ == '__main__':
	main(['aSyn5'], ['HCDR1', 'HCDR2'], ['R', 'K', 'H', 'STOP'], 6**8, 6, 0.33, False)
	