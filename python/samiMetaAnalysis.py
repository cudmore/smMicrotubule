import numpy as np

from scipy.stats import ttest_ind 

def myStats(a, b, name='', acceptLengthGreater=None):
	"""
	a: length
	b: tortuosity
	
	this has been modified so much I have lost track
	todo: rewrite
	"""
	decimalPlaces = 3
	c = a[~np.isnan(a)] # length
	d = b[~np.isnan(a)] # tort
	n_total = len(c)
	if acceptLengthGreater is not None:
		d = d[c>acceptLengthGreater] # tortuosity
		c = c[c>acceptLengthGreater] # length
	'''
	tortMean = np.mean(d[np.isfinite(d)])
	tortStd = np.std(d[np.isfinite(d)])
	'''
	
	# tort has inf when euclid=0 in len/euclid --->>> they will have different n
	d = d[np.isfinite(d)]
	tortMean = np.nanmean(d)
	tortStd = np.nanstd(d)
	print(',', name, ',', round(np.sum(c),decimalPlaces), ',', round(np.mean(c),decimalPlaces), ',', round(np.std(c),decimalPlaces), ',', len(c), ',', n_total, ',', round(tortMean,decimalPlaces), ',', round(tortStd,decimalPlaces), ',', len(d))
	
print('\n========= Branch Length =========')

#
# branch length
#

wtf = np.load('/Users/cudmore/box/data/sami/analysis/wt-female_branchLength.npy')
wtm = np.load('/Users/cudmore/box/data/sami/analysis/wt-male_branchLength.npy')
kof = np.load('/Users/cudmore/box/data/sami/analysis/ko-female_branchLength.npy')
kom = np.load('/Users/cudmore/box/data/sami/analysis/ko-male_branchLength.npy')

wtf_tortuosity = np.load('/Users/cudmore/box/data/sami/analysis/wt-female_tortuosity.npy')
wtm_tortuosity = np.load('/Users/cudmore/box/data/sami/analysis/wt-male_tortuosity.npy')
kof_tortuosity = np.load('/Users/cudmore/box/data/sami/analysis/ko-female_tortuosity.npy')
kom_tortuosity = np.load('/Users/cudmore/box/data/sami/analysis/ko-male_tortuosity.npy')

wtf_branchType = np.load('/Users/cudmore/box/data/sami/analysis/wt-female_branchType.npy')
wtm_branchType = np.load('/Users/cudmore/box/data/sami/analysis/wt-male_branchType.npy')
kof_branchType = np.load('/Users/cudmore/box/data/sami/analysis/ko-female_branchType.npy')
kom_branchType = np.load('/Users/cudmore/box/data/sami/analysis/ko-male_branchType.npy')


condStrList = ['wtf', 'wtm', 'kof', 'kom']
condList = [wtf, wtm, kof, kom]
tortuosityList = [wtf_tortuosity, wtm_tortuosity, kof_tortuosity, kom_tortuosity]
branchTypeList = [wtf_branchType, wtm_branchType, kof_branchType, kom_branchType]

# only accept branch length greater than acceptLengthGreater
# 20200407, this is now in um (I am saving length/tortuosity scaled)
acceptLengthGreater = 1 # 0.4 # use None to have to pruning
print('acceptLengthGreater:', acceptLengthGreater)

# junction-to-junction is 2
acceptBranchType = 2 # use None to have no pruning
print('acceptBranchType:', acceptBranchType)

# change between wt and knockout
#pChange_m = wtm

#
# print the mean/std/n
print(',', 'cond, len tot, len mean, len std, n, n_total, tort mean, tort std, tort n') # n_total is total before pairing down with acceptLengthGreater
for condIdx, cond in enumerate(condList):
	name = condStrList[condIdx]
	for j in range(cond.shape[1]):
		# this is going to be hard to read 2-3 months from now
		myStats(cond[:,j], tortuosityList[condIdx][:,j], name=name, acceptLengthGreater=acceptLengthGreater)
	print('')

#
# strip down by acceptLengthGreater and acceptBranchType
for condIdx, cond in enumerate(condList):
	# todo: do the same for tort, use np.isfinit()
	condStr = condStrList[condIdx]
	condTort = tortuosityList[condIdx]
	if acceptBranchType is not None:
		cond = cond[branchTypeList[condIdx]==acceptBranchType]
		condTort = condTort[branchTypeList[condIdx]==acceptBranchType]
	cond = cond[~np.isnan(cond)] # this flattens 2d to 1d !!!
	if acceptLengthGreater is not None:
		condTort = condTort[cond>acceptLengthGreater] # do this first
		cond = cond[cond>acceptLengthGreater]
	condList[condIdx] = cond 
	tortuosityList[condIdx] = condTort 
	# this is sloppy, I have 2x copies of each wtf..., one as a variable and another in a list?
	if condStr == 'wtf':
		wtf = cond
		wtf_tortuosity = condTort
	elif condStr == 'wtm':
		wtm = cond
		wtm_tortuosity = condTort
	elif condStr == 'kof':
		kof = cond
		kof_tortuosity = condTort
	elif condStr == 'kom':
		kom = cond
		kom_tortuosity = condTort
			
# 4x comparisons
# 1) wt male vs female
# 2) wt male vs knockout male
# 3) wt female vs knockout female
# 4) ko male vs ko female
decimalPlaces = 3
condStrList2 = [('wtm','wtf'), ('wtm','kom'), ('wtf','kof'), ('kom','kof')]
condList2 = [(wtm,wtf), (wtm,kom), (wtf,kof), (kom,kof)]
# todo: put variables in a dict accsed by wtm/wtf/... keys loadedTort['wtm']
condList3 = [(wtm_tortuosity,wtf_tortuosity), (wtm_tortuosity,kom_tortuosity), (wtf_tortuosity,kof_tortuosity), (kom_tortuosity,kof_tortuosity)]
for compareIdx, (a, b) in enumerate(condList2):
	compareStrTuple = condStrList2[compareIdx]
	fStat, pValue = ttest_ind(a, b)
	fStat = round(fStat,decimalPlaces)
	#pValue = round(pValue,decimalPlaces)
	print(condStrList2[compareIdx][0], 'vs', condStrList2[compareIdx][1])
	
	print(',', compareStrTuple[0], 'len', ',', round(np.mean(a),decimalPlaces), ',', round(np.std(a),decimalPlaces), ',', a.shape[0]) 
	print(',', compareStrTuple[1], 'len', ',', round(np.mean(b),decimalPlaces), ',', round(np.std(b),decimalPlaces), ',', b.shape[0], ',', 'f=', ',', fStat, ',', 'p=', ',', pValue) 
	#print('f=', round(fStat,decimalPlaces), ',', 'p=', pValue)

	# tort
	theseTort = condList3[compareIdx]
	a_ = theseTort[0]
	a_ = a_[np.isfinite(a_)]
	b_ = theseTort[1]
	b_ = b_[np.isfinite(b_)]

	fStat, pValue = ttest_ind(a_, b_)
	fStat = round(fStat,decimalPlaces)
	#pValue = round(pValue,decimalPlaces)
	print(',', compareStrTuple[0], 'tort', ',', round(np.mean(a_),decimalPlaces), ',', round(np.std(a_),decimalPlaces), ',', a_.shape[0]) 
	print(',', compareStrTuple[1], 'tort', ',', round(np.mean(b_),decimalPlaces), ',', round(np.std(b_),decimalPlaces), ',', b_.shape[0], ',', 'f=', ',', fStat, ',', 'p=', ',', pValue) 
	