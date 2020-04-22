	#
	# compile raw data into different pandas dataframes (dfList_*) and save them as .xlsx below
	#d = {}
	#dfList = [] # not used?
	dfList_branchLength = []
	dfList_euclideanDistance = []
	dfList_branchType = []
	dfList_tortuosity = []
	maxNumBranchLength = 0
	maxNumEuclideanDistance = 0
	maxNumBranchType = 0
	maxNumTortuosity = 0
	for idx, result in enumerate(resultsList): # resultsList can be across files
		# scale
		'''
		xVoxel = resultsList[idx]['xVoxel']
		yVoxel = resultsList[idx]['yVoxel']
		zVoxel = resultsList[idx]['zVoxel']
		'''
		
		# code is becoming sloppey, I use these to get tortuosity in k == 'euclideanDistance'
		thisBranchLength = None
		thisEuclideanDistance = None
		
		for k,v in result['data'].items():
			colName = k + '_' + str(idx)
			#print(idx, 'colName:', colName, type(v))
			#d[colName] = v
			#
			'''
			dfDict = {colName: v}
			df = pd.DataFrame(dfDict, index=range(len(v)))
			dfList.append(df)
			'''
			#
			if k == 'branchLength':
				#print('branchlength shape:', v.shape)
				if len(v) > maxNumBranchLength:
					maxNumBranchLength = len(v)
				#v = v * xVoxel # assuming x/y voxel is the same
				#vScaled = np.multiply(v, xVoxel)
				vScaled = v
				thisBranchLength = vScaled # used in k == 'euclideanDistance' to get tortuosity
				dfDict = {colName: vScaled}
				df = pd.DataFrame(dfDict, index=range(len(v)))
				dfList_branchLength.append(df)
			if k == 'euclideanDistance':
				#print('branchlength shape:', v.shape)
				if len(v) > maxNumEuclideanDistance:
					maxNumEuclideanDistance = len(v)
				#v = v * xVoxel # assuming x/y voxel is the same
				#vScaled = np.multiply(v, xVoxel)
				vScaled = v
				thisEuclideanDistance = vScaled # used in k == 'euclideanDistance' to get tortuosity
				dfDict = {colName: vScaled}
				df = pd.DataFrame(dfDict, index=range(len(v)))
				dfList_euclideanDistance.append(df)
		
			if k == 'tortuosity':
				# assuming we have both branchLength and euclideanDistance, calculate the scale tortuosity
				# this is REALLY bad style !!!
				if len(v) > maxNumTortuosity:
					maxNumTortuosity = len(v)
				'''
				# this will print 'divide by zero encountered in true_divide' and value will become inf
				vScaled = np.divide(thisBranchLength, thisEuclideanDistance) # might fail on divide by 0
				'''
				vScaled = v
				#colName = 'tortuosity' + '_' + str(idx)
				dfDict = {colName: vScaled}
				df = pd.DataFrame(dfDict, index=range(len(v)))
				#print('tort dict')
				#display(df)
				dfList_tortuosity.append(df)
			
			if k == 'branchType':
				#print('branchlength shape:', v.shape)
				if len(v) > maxNumBranchType:
					maxNumBranchType = len(v)
				dfDict = {colName: v} # no scaling
				df = pd.DataFrame(dfDict, index=range(len(v)))
				dfList_branchType.append(df)
	
	#
	# make 2d array of branch length to save numpy .npy and reopen for analysis
	# todo: take each of these 'for idx, result' and put into one loop
	print('    === making numpy arrays to save .npy')
	nResult = len(resultsList)
	npBranchLength = np.empty((maxNumBranchLength,nResult))  
	npBranchLength[:] = np.nan
	npEuclideanDistance = np.empty((maxNumEuclideanDistance,nResult))  
	npEuclideanDistance[:] = np.nan
	npBranchType = np.empty((maxNumBranchType,nResult))  
	npBranchType[:] = np.nan
	npTortuosity = np.empty((maxNumTortuosity,nResult))  
	npTortuosity[:] = np.nan
	for idx, result in enumerate(resultsList): # resultsList can be across files
		for k,v in result['data'].items():
			if k == 'branchLength':
				#v = v * xVoxel # assuming x/y voxel is the same
				#vScaled = np.multiply(v, xVoxel)
				vScaled = v
				npBranchLength[:len(v),idx] = vScaled
			if k == 'euclideanDistance':
				#vScaled = np.multiply(v, xVoxel)
				vScaled = v
				npEuclideanDistance[:len(v),idx] = vScaled
			if k == 'tortuosity':
				vScaled = v
				npTortuosity[:len(v),idx] = vScaled
				'''
				# sloppy, assuming branchLength comes BEFORE euclideanDistance
				# will print: divide by zero encountered in true_divide
				npTortuosity[:len(v),idx] = np.divide(npBranchLength[:len(v),idx], npEuclideanDistance[:len(v),idx])
				'''
				'''
				print('npTortuosity idx:', idx)
				print(npTortuosity[:,idx])
				'''
			if k == 'branchType':
				npBranchType[:len(v),idx] = v # no scaling
			# never gets here
			'''
			if k == 'tortuosity':
				# will print: divide by zero encountered in true_divide
				npTortuosity[:len(v),idx] = np.divide(npBranchLength[:,idx], npEuclideanDistance[:,idx])
			'''
	print('   these should all have same shape ...')
	# branch length
	print('    npBranchLength.shape:', npBranchLength.shape)
	branchLengthFile = savePath + '_branchLength.npy'
	np.save(branchLengthFile, npBranchLength)

	# branch euclidean distance
	print('    npEuclideanDistance.shape:', npEuclideanDistance.shape)
	euclideanDistanceFile = savePath + '_euclideanDistance.npy'
	np.save(euclideanDistanceFile, npEuclideanDistance)
	# branch type

	print('    npBranchType.shape:', npBranchType.shape)
	branchTypeFile = savePath + '_branchType.npy'
	np.save(branchTypeFile, npBranchType)

	# tortuosity (on save .npy file has all nan?
	print('    npTortuosity.shape:', npTortuosity.shape)
	'''
	print('    npTortuosity:', np.nanmean(npTortuosity))
	print(npTortuosity[:,0])
	'''
	tortuosityFile = savePath + '_tortuosity.npy'
	np.save(tortuosityFile, npTortuosity)
			
	#
	# save .xlsx files
	
	# save just branch length
	branchLengthFile = 'branchLength_summary.xlsx'
	df = pd.concat(dfList_branchLength, axis=1) 
	sheet_name = sheetName
	if os.path.isfile(branchLengthFile):
		mode = 'a'
	else:
		mode = 'w'
	print('saving xlsx file sheet_name:', sheet_name, 'mode:', branchLengthFile)
	with pd.ExcelWriter(branchLengthFile, mode=mode) as writer:
		df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)
	
	# save just branch euclidean
	branchEuclideanFile = 'branchEuclidean_summary.xlsx'
	df = pd.concat(dfList_euclideanDistance, axis=1) 
	sheet_name = sheetName
	if os.path.isfile(branchEuclideanFile):
		mode = 'a'
	else:
		mode = 'w'
	print('saving xlsx file sheet_name:', sheet_name, 'mode:', branchEuclideanFile)
	with pd.ExcelWriter(branchEuclideanFile, mode=mode) as writer:
		df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)
	
	# save just branch type
	branchTypeFile = 'branchType_summary.xlsx'
	df = pd.concat(dfList_branchType, axis=1) 
	sheet_name = sheetName
	if os.path.isfile(branchTypeFile):
		mode = 'a'
	else:
		mode = 'w'
	print('saving xlsx file sheet_name:', sheet_name, 'mode:', branchTypeFile)
	with pd.ExcelWriter(branchTypeFile, mode=mode) as writer:
		df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)
	
	# save just tortuosity
	tortuosityFile = 'branchTortuosity_summary.xlsx'
	df = pd.concat(dfList_tortuosity, axis=1) 
	sheet_name = sheetName
	print('saving xlsx file sheet_name:', sheet_name, 'mode:', tortuosityFile)
	if os.path.isfile(tortuosityFile):
		mode = 'a'
	else:
		mode = 'w'
	with pd.ExcelWriter(tortuosityFile, mode=mode) as writer:
		df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)
