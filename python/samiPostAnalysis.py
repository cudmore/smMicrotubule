"""
Robert Cudmore
20200418
"""

import os
import itertools
from collections import OrderedDict

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt

import scipy.stats

import seaborn as sns

class samiPostAnalysis:
	"""
	"""
	def __init__(self):
		wtfPath = '../analysis/wt-female_results.csv'
		kofPath = '../analysis/ko-female_results.csv'
		wtmPath = '../analysis/wt-male_results.csv'
		komPath = '../analysis/ko-male_results.csv'

		if not os.path.isfile(wtfPath):
			wtfPath = 'analysis/wt-female_results.csv'
			kofPath = 'analysis/ko-female_results.csv'
			wtmPath = 'analysis/wt-male_results.csv'
			komPath = 'analysis/ko-male_results.csv'
		
		print('samiPostAnalysis() is loading analysis files')
		# load as dataframe
		print('    loading:', wtfPath)
		wtf = pd.read_csv(wtfPath)
		print('    loading:', kofPath)
		kof = pd.read_csv(kofPath)
		print('    loading:', wtmPath)
		wtm = pd.read_csv(wtmPath)
		print('    loading:', komPath)
		kom = pd.read_csv(komPath)

		# strip out 'inf', todo: do htis in samiAnalysis
		wtf[wtf['tortuosity']==np.inf] = np.nan
		kof[kof['tortuosity']==np.inf] = np.nan
		wtm[wtm['tortuosity']==np.inf] = np.nan
		kom[kom['tortuosity']==np.inf] = np.nan

		dfList = [wtf, kof, wtm, kom]

		# append them all together
		self.df = pd.concat(dfList, axis=0, ignore_index=True)

		self.removeDictList = []
		
		self.genotypes = ['wt', 'ko']
		self.sexes = ['female', 'male']

		# this works ... IS NOT USED
		# a = list(itertools.starmap(_getCondition, itertools.product(self.genotypes, self.sexes)))
		self.condList = [] # ['wtf', 'wtm', kof', 'kom']
		for genotype in self.genotypes:
			for sex in self.sexes:
				cond = self._getCondition(genotype, sex)
				self.condList.append(cond)

		#print(self.df.head())
		
	def _getCondition(self, genotype, sex):
		"""
		return: ('wtf', 'wtm', 'kof', 'kom')
		"""
		if sex == 'female':
			sex = 'f'
		elif sex == 'male':
			sex = 'm'
		condition = genotype + sex
		return condition

	# not implemented
	def setRemove(self, removeDictList):
		"""
		removeDictList = [
							{
							'genotype': ('wt', 'ko'),
							'sex': ('female', 'male')
							'filename': fn
							}
						]
		
		There might be overlapping cell with same filename?
		Maybe use myCellNumber
		"""
		for removeDict in removeDictList:
			print(removeDict)
	
	def getPruned(self, pruneDict):
		"""
		return a copy of self.df after pruning
		
		This is our main analysis engine !!!
		"""
		
		df = self.df
		
		if pruneDict is None:
			return df
		
		statName = pruneDict['statName']

		# reduce by 'branchType'
		branchTypeList = pruneDict['branchType']
		if not isinstance(branchTypeList,list): branchTypeList = [branchTypeList]
		if branchTypeList:
			df = df[df['branchType'].isin(branchTypeList)]
		
		# reduce by genotype
		genotypeList = pruneDict['genotype']
		if not isinstance(genotypeList,list): genotypeList = [genotypeList]
		if genotypeList:
			df = df[df['genotype'].isin(genotypeList)]
		
		# reduce by sex
		sexList = pruneDict['sex']
		if not isinstance(sexList,list): sexList = [sexList]
		if sexList:
			df = df[df['sex'].isin(sexList)]
		
		# reduce by value
		minValue = pruneDict['minValue']
		if minValue is not None:
			df = df[df[statName]>=minValue]

		#
		return df
		
	def getCounts(self, pruneDict, asDict=False):
		"""
		return a dataframe with raw and cell count for each (genotype, sex)
		"""
		df = self.getPruned(pruneDict)
		
		rowList = []
		
		if asDict:
			retDict = OrderedDict()
		
		# this will yield 4x groups (wtf, wtm, kof, kom)
		for genotype, sex in itertools.product(self.genotypes, self.sexes):
			genotypeList = [genotype]
			sexList = [sex]
			
			df2 = df[df['genotype'].isin(genotypeList) & df['sex'].isin(sexList)]
			
			rawList = df2['len3d'].values # values within this group
			nRaw = rawList.shape[0]		

			cellList = df2['myCellNumber'].unique() # cell within this group
			nCell = len(cellList)		

			if asDict:
				condKey = self._getCondition(genotype, sex)
				retDict[condKey] = OrderedDict()
				retDict[condKey]['nRaw'] = nRaw
				retDict[condKey]['nCell'] = nCell
			else:
				rowDict = OrderedDict({'genotype':genotype, 'sex':sex, 'nRaw':nRaw, 'nCell':nCell})
				rowList.append(rowDict)
			
		if asDict:
			theRet = retDict
		else:
			dfRet = pd.DataFrame(rowList)     
			theRet = dfRet
			
		#
		return theRet
			
	def getCellMean(self, pruneDict, verbose=False):
		"""
		return a df with mean across cells
		"""
		
		statName = pruneDict['statName']
		
		# todo: for now, pretty sure this only works for one stat
		theseColumns = pruneDict['statName']
		statFuncList = [np.mean, np.std, scipy.stats.sem, 'count']
		
		# new
		'''
		statFuncList = {
			statName+'_': ['mean']
			}
		'''
		
		# todo: replace this by self.getPruned()
		df = self.getPruned(pruneDict)
		
		df = df.groupby(['genotype', 'sex', 'myCellNumber'])[theseColumns].agg(statFuncList)
		
		# never used this
		#df.columns = ["".join(x) for x in df.columns.ravel()]
		
		# this was working but columns become (mean,std,sem,count), e.g. we loose stat name
		df.columns = df.columns.map(''.join)
		df = df.reset_index()

		if verbose:
			print('getCellMean()', pruneDict)
			print(df)
			
		return df
		
	def myStatTest(self, pruneDict, v1, v2):
		"""
		Kruskal-Wallis H-test tests the null hypothesis that the
		population median of all of the groups are equal.
		It is a non-parametric version of ANOVA. The test works on 2 or more independent samples,
		which may have different sizes.
		
		see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kruskal.html
			h,p = scipy.stats.kruskal(v1, v2, nan_policy='omit')
		
		Mann-Whitney rank test on samples x and y.

		see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html#scipy.stats.mannwhitneyu
			u, p = scipy.stats.mannwhitneyu(v1, v2, alternative='two-sided') # {None, ‘two-sided’, ‘less’, ‘greater’}
		
		"""
		v1 = v1[np.isfinite(v1)]
		v2 = v2[np.isfinite(v2)]
		
		statTest = pruneDict['statTest']
		
		t = None
		prob = None
		
		try:
			if statTest == 'T-Test':
				t, prob = scipy.stats.ttest_ind(v1, v2) # this is a two-tailed test
			elif statTest == 'Mann-Whitney':
				t,prob = scipy.stats.mannwhitneyu(v1, v2, alternative='two-sided')
			elif statTest == 'Kruskal-Wallis':
				t,prob = scipy.stats.kruskal(v1, v2, nan_policy='omit')
			else:
				print('warning: myStatTest() defaulting to T-Test, unknown statTest:', statTest)
				t, prob = scipy.stats.ttest_ind(v1, v2) # this is a two-tailed test
		except (ValueError) as e:
			pass
						
		return t, prob
		
	def getPairwiseCellComparison(self, pruneDict):
		"""
		todo: just use getPairwiseGroupComparison() ???
		
		Assuming one (genotype,sex), make statistical comparison between each cell mean
		
		For example:	
			pruneDict['genotype'] = ['wt']
			pruneDict['sex'] = ['male']

		The reason I am doing this is I am worried
		    most cells within a given (genotype, sex) will have significant differences?
		
		todo: perform one-tailed, scipy.stats.ttest_ind is two-tailed
		todo: pre-allocate each cells raw values so we do not duplicate fetching them in nested loop
		"""
		
		statName = pruneDict['statName']
		
		# reduce/prune down to just one group, e.g. (wt, m)
		df = self.getPruned(pruneDict)
		
		cellList = df['myCellNumber'].unique() # cell within this group
		n = len(cellList)

		rowList = []
		
		for i in range(n):
			iValues = df[df['myCellNumber'] == i][statName].values
			rowDict = {'cells':i} # first column is cell in a n*n+1 matrix
			for j in range(n):
				if j>=i:
					jValues = df[df['myCellNumber'] == j][statName].values
					# run t-test
					#t, prob = scipy.stats.ttest_ind(iValues, jValues) # this is a two-tailed test
					t, prob = self.myStatTest(pruneDict, iValues, jValues)
					
					rowDict[j] = prob
				else:
					rowDict[j] = ''
					
			# after j loop, inside i loop
			rowList.append(rowDict)
								
		dfRet = pd.DataFrame(rowList)     
		#
		return dfRet
		
	def _getStatStr(self, pruneDict, name1, values1, name2=None, values2=None, asDict=True):

		rowStr = ''
		rowDictList = []
		
		values1 = values1[np.isfinite(values1)]
		if values2 is not None:
			values2 = values2[np.isfinite(values2)]
		
		mean = float(np.mean(values1))
		std = float(np.std(values1))
		sem = float(scipy.stats.sem(values1))
		n = int(values1.shape[0])
		if asDict:
			rowDict = OrderedDict()
			rowDict['name'] = name1
			rowDict['mean'] = mean	
			rowDict['std'] = std	
			rowDict['sem'] = sem	
			rowDict['n'] = n
			rowDict['f'] = ''
			rowDict['p'] = ''
			rowDictList.append(rowDict)
			
		else:
			rowStr += '{0},{1},{2},{3},{4}'.format(name1, mean, std, sem, n)
		
		if name2 is not None and values2 is not None:
		
			#t, prob = scipy.stats.ttest_ind(values1, values2) # this is a two-tailed test
			t, prob = self.myStatTest(pruneDict, values1, values2)
			
			mean = float(np.mean(values2))
			std = float(np.std(values2))
			sem = float(scipy.stats.sem(values2))
			n = int(values2.shape[0])

			if asDict:
				rowDict = OrderedDict()
				rowDict['name'] = name2
				rowDict['mean'] = mean	
				rowDict['std'] = std	
				rowDict['sem'] = sem	
				rowDict['n'] = n	
				rowDict['f'] = t	
				rowDict['p'] = prob	
				rowDictList.append(rowDict)
			else:
				rowStr += '\n'
				rowStr += '{0},{1},{2},{3},{4},f={5},p={6}'.format(name2, mean, std, sem, n, t, prob)

		if asDict:
			return rowDictList
		else:
			return rowStr
		
	def getPairwiseGroupComparison(self, pruneDict, doCellMean=False):
		"""
		Get all pairwise comparisons between our different 'groups'.
		For a 4x4 matrix, we get the upper right triangle including diagonal
		
			[('wt', 'f'), ('wt', 'm'), ('ko', 'f'), ('ko', 'm')]
		
		This pools all measurements across all cells (does not group by cell)
		
		todo: rewrite this to use per cell mean
		
		IMPORTANT, if I chase pruneDict, it propogates back to caller????
		"""
		
		localPruneDict = pruneDict.copy()
		
		if doCellMean:
			myStat = 'mean'
		else:
			myStat = localPruneDict['statName']

		# there are 4 groups
		groupsList = list(itertools.product(self.genotypes, self.sexes))
		n = len(groupsList)

		rowList = []
		rowStrList = [] # to output comma separated (genotype,sex), m, sd,se, n, ...p
		
		for iIdx, (iGenotype,iSex) in enumerate(itertools.product(self.genotypes, self.sexes)): 
			localPruneDict['genotype'] = iGenotype
			localPruneDict['sex'] = iSex
			if doCellMean:
				dfi = self.getCellMean(localPruneDict, verbose=False) # need to use 'mean' as stat
			else:
				dfi = self.getPruned(localPruneDict)
			iValues = dfi[myStat].values

			#
			rowDict = {'groups':(iGenotype,iSex)} # first column is group name, making a 4x4 df
			
			for jIdx, (jGenotype,jSex) in enumerate(itertools.product(self.genotypes, self.sexes)): 
				if jIdx >= iIdx:
					localPruneDict['genotype'] = jGenotype
					localPruneDict['sex'] = jSex
					if doCellMean:
						dfj = self.getCellMean(localPruneDict, verbose=False) # need to use 'mean' as stat
					else:
						dfj = self.getPruned(localPruneDict)
					jValues = dfj[myStat].values
					#
					#t, prob = scipy.stats.ttest_ind(iValues, jValues) # this is a two-tailed test
					'''
					print('getPairwiseGroupComparison()', iIdx, jIdx, 'i', iValues.shape, 'j', jValues.shape)
					if iIdx==0 and jIdx==2:
						print('iValues:')
						print(iValues)
						print('jValues:')
						print(jValues)
					'''
					t, prob = self.myStatTest(localPruneDict, iValues, jValues)
					#
					rowDict[(jGenotype,jSex)] = prob

					rowDict2 = self._getStatStr(localPruneDict, (iGenotype,iSex), iValues, (jGenotype,jSex), jValues, asDict=True)
					rowStrList += rowDict2
				else:
					rowDict[(jGenotype,jSex)] = ''
			# after j loop, inside i loop
			rowList.append(rowDict)
								
		dfRet = pd.DataFrame(rowList)     
		dfRet2 = pd.DataFrame(rowStrList)     
		
		#
		return dfRet, dfRet2
		
	def getDefaultPruneDict(self):
		pruneDict  = OrderedDict()
		pruneDict['genotype'] = [] # [] means all
		pruneDict['sex'] = [] # [] means all
		pruneDict['branchType'] = [2] # None means all
		pruneDict['useRemove'] = True # todo: implement this to exclude individual cells
		pruneDict['statName'] = 'len3d' # corresponds to numeric column in self.df
		pruneDict['minValue'] = 1 # if specified will only accept values >= ['minValue']
		pruneDict['doCellMean'] = False # if true, prunng will return mean for each cell, rather than raw branch data
		pruneDict['statTest'] = 'Mann-Whitney' # one of ('T-Test', 'Mann-Whitney', 'Kruskal-Wallis')
		return pruneDict
	
	def _plotCondMeanLegend(self, pruneDict, ax, doCellMean):
		"""
		getting sloppy here
		
		see: https://stackoverflow.com/questions/58267394/seaborn-barplot-add-xticks-for-hue
		"""
		nDict = self.getCounts(pruneDict, asDict=True)
		if doCellMean:
			wtf_n = nDict['wtf']['nCell']
			wtm_n = nDict['wtm']['nCell']
			kof_n = nDict['kof']['nCell']
			kom_n = nDict['kom']['nCell']
		else:
			wtf_n = nDict['wtf']['nRaw']
			wtm_n = nDict['wtm']['nRaw']
			kof_n = nDict['kof']['nRaw']
			kom_n = nDict['kom']['nRaw']

		xAxisLabels = ['wt f\nn='+str(wtf_n), 'wt m\nn='+str(wtm_n), 'ko f\nn='+str(kof_n), 'ko m\nn='+str(kom_n)]
		ax.set_xticks([-0.2,0.2, 0.8,1.2])
		ax.set_xticklabels(xAxisLabels)
		
	def plotCondMean(self,pruneDict=None, doCellMean=None, ax=None):
		"""
		all branch length pooled, todo: do this per cell
		violon plot across genotype, sex is color (x) and pruneDict['statname'] (y)
		
		TODO: make new version of this to plot violin over mean per cell !!!
			see jupyter notebook for code !!!
		"""
				
		if doCellMean:
			myStatName = 'mean'
			df = self.getCellMean(pruneDict, verbose=False) # need to use 'mean' as stat
		else:
			myStatName = pruneDict['statName']
			df = self.getPruned(pruneDict)
		
		xOrder = self.genotypes
		hueOrder = self.sexes
		hue = 'sex'
				
		if ax is None:
			fig, ax = plt.subplots(1, 1, figsize=(10,5))
			
		#
		# violin
		g = sns.violinplot(x="genotype", y=myStatName, order=xOrder, hue_order=hueOrder, hue=hue, data=df, ax=ax)
		
		xlabels = ['{:,.1f}'.format(x) for x in g.get_xticks()]
		g.set_xticklabels(xlabels)

		#
		# strip
		# when kind="strip"  we also need dodge=True to seperate hue (e.g. sex)
		g = sns.stripplot(x="genotype", y=myStatName,
			order=xOrder, hue_order=hueOrder, hue=hue,
			dodge=True, palette=['#91bfdb','#fc8d59'], data=df, ax=ax)
		
		myFontSize = 16

		g.legend(bbox_to_anchor=(1, 1), ncol=1)
		# i want to remove duplicate legends???
		# this removes the entire legend?
		# see: diff b/w/ ._legend and .legend_ https://stackoverflow.com/questions/54781243/hide-legend-from-seaborn-pairplot
		#g.legend_.remove()
		
		
		# title
		statName = pruneDict['statName']
		if doCellMean:
			ax.set_title('Mean {0} across mean of each cell'.format(statName), fontsize=myFontSize)
		else:
			ax.set_title('Mean {0} across all cells'.format(statName), fontsize=myFontSize)
				
		g.set_xlabel(g.get_xlabel(), size = myFontSize)
		g.set_ylabel(g.get_ylabel(), size = myFontSize)
		#
		g.set_xticklabels(g.get_xticks(), size = myFontSize)
		g.set_yticklabels(g.get_yticks(), size = myFontSize)

		# yticklabels were sometimes blowing up into 1.000000000001
		ylabels = [str(round(y,1)) for y in g.get_yticks()] # used below
		g.set_yticklabels(ylabels, size = myFontSize) # using from above

		self._plotCondMeanLegend(pruneDict, ax, doCellMean)

		#
		return g

	def plotHist(self, pruneDict=None, ax=None):
		"""
		plot a histogram of ['statName']

		use lists for pruneDict['genotype'] and ['sex'] to plot multiple conditions
		
		todo: make per cell histogram
		"""
		
		localPruneDict = pruneDict.copy()

		genotypes = localPruneDict['genotype']
		sexes = localPruneDict['sex']
		statName = localPruneDict['statName']
		
		defaultColorList = ['k', 'r', 'g', 'b']
		
		# step through genotype/sex and make multiple raw datasets
		#df = self.getPruned(localPruneDict)
		valuesList = []
		colors = []
		labels = [] # for legend
		idx = 0
		for genotype in genotypes:
			for sex in sexes:
				# pull values
				localPruneDict['genotype'] = [genotype]
				localPruneDict['sex'] = [sex]
				df = self.getPruned(localPruneDict)
				values = df[statName].values
				#
				valuesList.append(values)
				colors.append(defaultColorList[idx])
				labelStr = genotype + ' ' + sex # for legend
				labels.append(labelStr)
				# increment
				idx += 1
		
		if ax is None:
			fig, ax = plt.subplots(1, 1, figsize=(10,5))

		ax.hist(valuesList, bins='auto', density=True, histtype='bar', color=colors, label=labels) #, alpha=0.7, rwidth=0.85)
		ax.legend(prop={'size': 10})
		
		return ax # call plt.show() to show !!!
		
if __name__ == '__main__':
	spa = samiPostAnalysis()
	
	if 0:
		pruneDict = spa.getDefaultPruneDict()
		pruneDict['sex'] = ['female']
		pruneDict['statName'] = 'len3d'
		cellMean = spa.getCellMean(pruneDict, verbose=True)
		#print(cellMean)
	
	'''
	pruneDict = spa.getDefaultPruneDict()
	pruneDict['branchType'] = [2]
	cellMean = spa.getCellMean(pruneDict)
	print(cellMean)
	'''
	
	# see for multiple seaborn plots
	# http://seaborn.pydata.org/introduction.html#figure-level-and-axes-level-functions
	if 0:
		fig, axs = plt.subplots(1, 2, figsize=(10,5))
		axs = np.ravel(axs) # flatten

		pruneDict = spa.getDefaultPruneDict()
		pruneDict['statName'] = 'len3d'
		pruneDict['minValue'] = 2
		pruneDict['branchType'] = [2]

		# plot violin across all segment lengths
		g = spa.plotCondMean(pruneDict, doCellMean=False, ax=axs[0])

		# per cell violin
		g = spa.plotCondMean(pruneDict, doCellMean=True, ax=axs[1])
		plt.show()
			
	# for a given (genotype, sex), make
	# (1) tables of mean for each cell
	# (2) p-value matrix holding all pair-wise t-test(s) between cells
	if 0:
		pruneDict = spa.getDefaultPruneDict()
		pruneDict['genotype'] = ['wt']
		pruneDict['sex'] = ['male']
		pruneDict['statName'] = 'len3d'
		pruneDict['minValue'] = 2
		pruneDict['branchType'] = [2]

		'''
		# cell mean table
		cellMean = spa.getCellMean(pruneDict)
		print(cellMean)
		'''
		
		# p-value matrix
		dfResult = spa.getPairwiseCellComparison(pruneDict)
		print(dfResult)
		#with np.printoptions(precision=5, suppress=True): # force a pretty-print
		#	print(result)
	
	if 1:
		pruneDict = spa.getDefaultPruneDict()
		print('pruneDict:', pruneDict)
		pruneDict['branchType'] = [2]
		df1, df2 = spa.getPairwiseGroupComparison(pruneDict) # this was changing prune dict???
		print('pruneDict:', pruneDict)
		print('df1:')
		print(df1)
		print('df2:')
		print(df2)
		
	if 0:
		pruneDict = spa.getDefaultPruneDict()
		pruneDict['minValue'] = 1
		pruneDict['branchType'] = [2]
		df = spa.getCounts(pruneDict, asDict=False)
		print('getCounts():')
		print(df)
		
	if 0:
		pruneDict = spa.getDefaultPruneDict()
		pruneDict['statName'] = 'len3d'
		pruneDict['minValue'] = 1
		pruneDict['branchType'] = [2]
		pruneDict['statTest'] = 'T-Test'
		df1, df2 = spa.getPairwiseGroupComparison(pruneDict, doCellMean=True)
		