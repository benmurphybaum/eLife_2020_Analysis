#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


//Get ROI ------------------------------------
//Extracts the time-varying Calcium signal from pre-defined ROI coordinates.
//The scans and ROIs to be analyzed are selected from list boxes in a GUI called 'NT'.

//Used for FIGURES 2, 4, 6, 8
Function/WAVE NT_GetROI()
	STRUCT IMAGING img
	
	DFREF TP = root:Packages:twoP:examine
	DFREF NTF = root:Packages:NT
	DFREF NTD = root:Packages:NT:DataSets
	
	//Gets the parameters, scan data, and ROI data from the GUI (NT)
	initParam(img)
	
	//Make the ROI analysis folder if it doesn't already exist
	If(!DataFolderExists("root:ROI_analysis"))
		NewDataFolder root:ROI_analysis
	EndIf
	SetDataFolder root:ROI_analysis	
	
	//Make the output wave reference wave for passing the result onto another function
	Make/FREE/WAVE/N=0 outputWaveRefs
	
	Variable i,j,k,totalWaveCount = 0
	
	//Get the scan and background waves
	For(i=0;i<img.scan.num;i+=1)
		Variable ref = StartMSTimer
		
		switch(img.channel)
			case 1: //channel 1 only
				Wave theScan = img.scan.ch1[i] //signal fluorescence
				Wave theBgnd = img.scan.ch1[i] //background fluorescence
				break
			case 2: //channel 2 only
				Wave theScan = img.scan.ch2[i]
				Wave theBgnd = img.scan.ch2[i]
				break
			case 3: // ch1 / ch2
				Wave theScan = img.scan.ch1[i]
				Wave theBgnd = img.scan.ch2[i]
				break
			case 4: // ch2 / ch1
				Wave theScan = img.scan.ch2[i]
				Wave theBgnd = img.scan.ch1[i] 
				break
		endswitch
		
		//Get dendritic mask -- returns an image mask of the dendrites
		Wave mask = GetDendriticmask(theBgnd)
		Redimension/B/U mask
		
		//Get dark value from the area not within the dendritic mask
		ImageStats/R=mask theBgnd
		Variable darkVal = 0.9*V_avg
		
		For(j=0;j<img.roi.num;j+=1)
			String theROI = img.rois[j]
			
			
			String ROIFolder = "root:ROI_analysis:" + theROI
			
			If(!DataFolderExists(ROIFolder))
				NewDataFolder $ROIFolder
			EndIf
			
			//X and Y waves that define the ROI area
			Wave roiX = img.roi.x[j]
			Wave roiY  = img.roi.y[j]
		
			//Seed values for filling out the ROI mask
			Variable maskMax,maskMin,xSeed,ySeed
			WaveStats/Q theBgnd
			
			maskMin = WaveMin(roiX)
			maskMax = WaveMax(roiX)
			
			xSeed = maskMax + DimDelta(theBgnd,0)
			If(xSeed > IndexToScale(theBgnd,DimSize(theBgnd,0)-1,0))
				xSeed = IndexToScale(theBgnd,0,0)
			EndIf
			
			maskMin = WaveMin(roiY)
			maskMax = WaveMax(roiY)
			
			ySeed = maskMax + DimDelta(theBgnd,1)
			If(ySeed > IndexToScale(theBgnd,DimSize(theBgnd,1)-1,1))
				ySeed = IndexToScale(theBgnd,0,1)
			EndIf
			
			//Generate the ROI mask wave	
			SetDataFolder $ROIFolder			
			ImageBoundaryToMask ywave=roiY,xwave=roiX,width=(DimSize(theBgnd,0)),height=(DimSize(theBgnd,1)),scalingwave=theBgnd,seedx=xSeed,seedy=ySeed			
		
			Wave ROIMask = $(ROIFolder + ":M_ROIMask")	
			
			//Did the ROI mask actually get created?
			If(!WaveExists(ROIMask))
				DoAlert 0, "Couldn't find the ROI mask wave for: " + NameOfWave(theScan)
				continue
			EndIf
			
			//Make the raw ROI waves for signal and background
			Variable numFrames = DimSize(theScan,2)
			Make/O/FREE/N=(numFrames) ROI_Signal,ROI_Bgnd
			
			//Set all the scales of the ROI waves
			SetScale/P x,DimOffset(theScan,2),DimDelta(theScan,2),ROI_Signal,ROI_Bgnd
			
			//Average values over the ROI region
			For(k=0;k<numFrames;k+=1)
				ImageStats/M=1/P=(k)/R=ROImask theScan
				ROI_Signal[k] = V_avg
				
				ImageStats/M=1/P=(k)/R=ROImask theBgnd
				ROI_Bgnd[k] = V_avg
			EndFor		
			
			//Savitzky-Golay smoothing
			Smooth/S=2 (img.filter), ROI_Signal
			
			//Use median for the baseline, so it doesn't get pulled up or down from noisy values
			Variable	bsln = median(ROI_Bgnd,img.bsSt,img.bsEnd)
			
			//Absolute fluorescence or delta fluorescence?
			If(img.mode == 1)
			//∆F/F
				String outName = NameOfWave(theScan) + "_" + theROI + "_dF"
			ElseIf(img.mode == 2)
			//Abs
				outName = NameOfWave(theScan) +  theROI + "_" + "_abs"
			EndIf	
			
			//Make the dF wave
			Make/O/N=(numFrames) $outName
			Wave dF = $outName
			
			//Set all the scales of the dF waves
			SetScale/P x,DimOffset(theScan,2),DimDelta(theScan,2),dF
			
			//Calculate the ∆F/F or Absolute fluoresence ratios
			If(img.mode == 1) //dF
				dF = (ROI_Signal - bsln) / (bsln - darkVal)
			ElseIf(img.mode == 2) //abs
				dF = ROI_Signal
			EndIf
						
			//These are all the output ROI waves
			Redimension/N=(totalWaveCount + 1) outputWaveRefs
			outputWaveRefs[totalWaveCount] = dF
			totalWaveCount += 1
		EndFor
		
		print "Get ROI:",NameOfWave(theScan) + ",",StopMSTimer(ref) / (1e6),"s"
		
	EndFor
	
	//pass the output wave references on
	return outputWaveRefs
End


//Returns a mask wave for the input scan
//Used by GetROI()
Function/WAVE GetDendriticMask(theWave)
	Wave theWave //calcium imaging scan
	
	//Get max projection image
	MatrixOP/S/O maxProj = sumBeams(theWave)
	ImageStats maxProj
	
	//Min/max of the image
	Variable minVal = V_min
	Variable maxVal = WaveMax(maxProj)
	
	//Simple value thresholding based on the min/max values
	Variable threshold = minVal + (maxVal - minVal) * 0.05 //1.25 is a mask threshold and can be changed
	Multithread maxProj = (maxProj < threshold) ? 0 : maxProj
	
	//Eliminate isolated points to reduce noise
	Make/FREE/N=(5,5) block
	block = 0
	
	Variable rows,cols,i,j
	
	rows = DimSize(maxProj,0)
	cols = DimSize(maxProj,1)
	
	For(i=0;i<rows;i+=1)	
		For(j=0;j<cols;j+=1)
			//skip zeros
			If(maxProj[i][j] == 0)
				continue
			EndIf			
	
			//check for image edges
			If(i-2 < 0 || i+2 > rows-2 || j-2 < 0 || j+2 > cols-2)
				continue
			Else
				//Get data block surrounding point
				block = maxProj[i-2 + p][j-2 + q]
				
				//Check for isolated point and remove
				If(sum(block) < 3*maxProj[i][j])
					maxProj[i][j] = 0
				EndIf	
			EndIf
		
			block = 0
		EndFor
	EndFor
	
	//2D median filter 3x3
	MatrixFilter/N=3 median maxProj
	
	//Create mask wave
	String maskName = NameOfWave(theWave) + "_mask"
	If(strlen(maskName) > 31)
		maskName = "Scan_mask"
	EndIf
	
	Make/O/N=(rows,cols)/FREE theMask
	MultiThread theMask = (maxProj == 0) ? 0 : 1
	
	//Scaling
	SetScale/P x,DimOffset(theWave,0),DimDelta(theWave,0),theMask
	SetScale/P y,DimOffset(theWave,1),DimDelta(theWave,1),theMask
	
	return theMask
End


//Returns ROI parameters to the calling function (Get ROI)
Function initParam(img)
	STRUCT IMAGING &img
	
	DFREF RF = root:twoP_ROIS
	DFREF TP = root:Packages:twoP:examine
	DFREF NTI = root:Packages:NT:Imaging
	
	//These controls are SetVariable inputs in separate GUI panel named 'NT'.
	ControlInfo/W=NT channelSelect
	img.channel = V_Value
	
	ControlInfo/W=NT dFSelect
	img.mode = V_Value
	
	ControlInfo/W=NT baselineSt
	img.bsSt = V_Value
	
	ControlInfo/W=NT baselineEnd
	img.bsEnd = V_Value
	
	ControlInfo/W=NT peakSt
	img.pkSt = V_Value
	
	ControlInfo/W=NT peakEnd
	img.pkEnd = V_Value
	
	ControlInfo/W=NT filterSize
	img.filter = V_Value
	
	//ROI ListBox list and select waves
	Wave/T ROIListWave = NTI:ROIListWave
	Wave ROISelWave =  NTI:ROISelWave
	
	Wave/T ScanListWave = NTI:ScanListWave
	Wave ScanSelWave = NTI:ScanSelWave
	
	img.roi.num = sum(ROISelWave)
	img.scan.num = sum(ScanSelWave)
	
	//active ROIs used for the analsis and their position wave references
	Make/O/N=(img.roi.num)/T NTI:ROI_List_Analysis
	Make/O/N=(img.roi.num)/WAVE NTI:ROI_Coord_X
	Make/O/N=(img.roi.num)/WAVE NTI:ROI_Coord_Y
	
	Wave/T img.rois = NTI:ROI_List_Analysis
	Wave/WAVE img.roi.x = NTI:ROI_Coord_X
	Wave/WAVE img.roi.y = NTI:ROI_Coord_Y
	
	//active Scans channels used for the analysis
	Make/O/N=(img.scan.num)/WAVE NTI:Scan_List_Ch1
	Make/O/N=(img.scan.num)/WAVE NTI:Scan_List_Ch2
	Wave/WAVE/Z img.scan.ch1 = NTI:Scan_List_Ch1
	Wave/WAVE/Z img.scan.ch2 = NTI:Scan_List_Ch2
	
	//Check that there is a selection at all for scans and rois
	Variable i = 0
	If(DimSize(ROISelWave,0) == 0 || img.roi.num == 0 || img.scan.num == 0)
		Redimension/N=0 img.rois,img.scan.ch1,img.scan.ch2,img.roi.x,img.roi.y
		return 0
	EndIf
	
	//Fill out all the ROI name and get their position waves
	Variable count = 0
	Do
		If(ROISelWave[i] == 1)
			img.rois[count] = ROIListWave[i]
			img.roi.x[count] = RF:$(img.rois[count] + "_x")
			img.roi.y[count] = RF:$(img.rois[count] + "_y")
			count += 1
		EndIf
		i += 1
	While(i < DimSize(ROISelWave,0))
	
	//Fill out all the scan waves
	count = 0;i = 0
	Do
		If(ScanSelWave[i] == 1)
			img.scan.ch1[count] = $("root:twoP_Scans:" + ScanListWave[i] + ":" + ScanListWave[i] + "_ch1")
			img.scan.ch2[count] = $("root:twoP_Scans:" + ScanListWave[i] + ":" + ScanListWave[i] + "_ch2")
			count += 1
		EndIf
		i += 1
	While(i < DimSize(ScanSelWave,0))
	
End


//Holds parameters of the Scans and ROIs for call by functions
Structure IMAGING
	STRUCT ROI roi
	STRUCT SCAN scan
	uint16 channel
	uint16 mode
	uint16 bsSt
	uint16 bsEnd
	uint16 pkSt
	uint16 pkEnd
	uint16 filter
	Wave/T rois
EndStructure

Structure ROI
	Wave/WAVE x
	Wave/WAVE y
	uint16 num
EndStructure

Structure SCAN
	Wave/WAVE ch1
	Wave/WAVE ch2
	uint16 num
EndStructure

//VectorSum-----------
//Calculates a vector sum of the input wave, and returns the specified value (angle, radius, or DSI)
//Used for FIGURES 1, 3S2, 4, 5, 7, 8
Function VectorSum(inputWave,doPrint,returnItem,[scaled,angleWave,PN])
	Wave inputWave //tuning curve wave
	Variable doPrint //print the results
	String returnItem //'vAngle', 'DSI', or 'vRadius'
	Variable scaled //x scaling is the angle
	Wave angleWave //supply an angle wave
	Variable PN //preferred null Vector Sum
	
	SetDataFolder GetWavesDataFolder(inputWave,1)
	
	If(!DataFolderExists("root:var"))
		NewDataFolder root:var
	EndIf
	
	If(ParamIsDefault(angleWave))
	//angle not wave supplied
		Make/O/N=8 root:var:direction
		Wave angleWave = root:var:direction
		If(ParamIsDefault(scaled))
			//not scaled, assue 45° delta angle
			angleWave = 45*x
		Else
			Redimension/N=(DimSize(inputWave,0)) angleWave
			angleWave = DimOffset(inputWave,0) + DimDelta(inputWave,0) * x
		EndIf
	EndIf
	
	If(ParamIsDefault(PN))
		PN = 0
	Else
		PN = 1
	EndIf
	
	//PN vector sum, don't use full tuning curve
	If(PN == 1)
		Redimension/N=2 angleWave
		angleWave = 180 * x
	EndIf
	
	//User error check.
	If(cmpstr(returnItem,"vAngle") !=0 && cmpstr(returnItem,"vRadius") !=0 && cmpstr(returnItem,"DSI") !=0)
		DoAlert 0,"Must indicate return value of 'vAngle','vRadius', or 'DSI'."
		return -1
	EndIf
	
	Variable i,j,numCols
	Variable vSumX,vSumY,totalSignal
	Variable numAngles = DimSize(angleWave,0)
	
	numCols = DimSize(inputWave,1)
	Make/FREE/N=(numCols) angles,dsi_cols
	If(DimSize(angles,0) == 0)
		numCols += 1
		Redimension/N=(numCols) angles,dsi_cols
	EndIf
	angles = 0
	
	For(j=0;j<numCols;j+=1)
		//loop through each column of the input wave, in case there are multiple tuning curves, one per column
	
		//get data from each column of input wave
		Make/FREE/N=(DimSize(inputWave,0)) data
		data[][0] = inputWave[p][j]
		SetScale/P x,DimOffset(inputWave,0),DimDelta(inputWave,0),data
		
		vSumX = 0
		vSumY = 0
		totalSignal = 0
		
		Variable nullPt
		
		If(PN)
			//PN vector sum
			WaveStats/Q data
			nullPt = polarMath2(V_maxLoc,180,"deg","add")
			If(nullPt == 360)
				nullPt = 0
			EndIf
		
			nullPt = ScaleToIndex(data,nullPt,0)//180° off of the max value direction (preferred)
			
			vSumX += data[V_maxRowLoc] * cos(angleWave[0]*pi/180)
			vSumY += data[V_maxRowLoc] * sin(angleWave[0]*pi/180)
			totalSignal += data[V_maxRowLoc]
			
			vSumX += data[nullPt] * cos(angleWave[1]*pi/180)
			vSumY += data[nullPt] * sin(angleWave[1]*pi/180)
			totalSignal += data[nullPt]
			
		Else
			//full tuning curve vector sum
			For(i=0;i<numAngles;i+=1)
				If(numtype(data[i]) == 2) 
					continue
				EndIf
				vSumX += data[i]*cos(angleWave[i]*pi/180)
				vSumY += data[i]*sin(angleWave[i]*pi/180)
				totalSignal += data[i]
			EndFor
		EndIf
		
		Variable vRadius = sqrt(vSumX^2 + vSumY^2)
		Variable vAngle = -atan2(vSumY,vSumX)*180/pi
		Variable	DSI = vRadius/totalSignal
		
		If(vAngle < 0)
			vAngle +=360
		Endif
		
		vAngle = 360 - vAngle
		
		angles[j] = vAngle
		dsi_cols[j] = DSI
		
		If(doPrint)
			print "vAngle =",vAngle,"\r  vRadius =",vRadius,"\r  DSI =",DSI	
		EndIf
		
	EndFor
	
	If(cmpstr(returnItem,"vAngle") == 0)
		If(DimSize(inputWave,1) > 0)
			Make/O/N=(DimSize(inputWave,1)) vAng_columns
			Wave vAng_cols = vAng_columns
			vAng_cols = angles
		EndIf
		return vAngle
	ElseIf(cmpstr(returnItem,"vRadius") == 0)
		return vRadius
	ElseIf(cmpstr(returnItem,"DSI") == 0)
		If(DimSize(inputWave,1) > 0)
			Make/O/N=(DimSize(inputWave,1)) vDSI_columns
			Wave vDSI_cols = vDSI_columns
			vAng_cols = dsi_cols
		EndIf
		return DSI
	EndIf
	
End Function

//Mathematical operations using polar coordinates.
Function polarMath2(pnt1,pnt2,degrad,op)
	Variable pnt1,pnt2
	String degrad,op
	Variable angOut
	
	strswitch(op)
		case "add":
			angOut = pnt1 + pnt2
			break
		case "distance":
			Variable x1,y1,x2,y2,D,A
			 //linear distance between the points
			
			D = 2*pi*1*A/360
			
			If(!cmpstr(degrad,"deg"))
				x1 = cos(pnt1 * pi/180)
				y1 = sin(pnt1 * pi/180)
				x2 = cos(pnt2 * pi/180)
				y2 = sin(pnt2 * pi/180)
				D = sqrt( (x2 - x1)^2 + (y2 - y1)^2 )
				angOut = acos( (2 * (1^2) - D^2) / (2 * (1^2)) ) * 180/pi
			Else
				x1 = cos(pnt1)
				y1 = sin(pnt1)
				x2 = cos(pnt2)
				y2 = sin(pnt2)
				D = sqrt( (x2 - x1)^2 + (y2 - y1)^2 )
				angOut = acos( (2 * (1^2) - D^2) / (2 * (1^2)) )
			EndIf
			break
	endswitch
	
	strswitch(degrad)
		case "deg":
			angOut = (angOut > 360) ? (angOut - 360) : angOut
			angOut = (angOut < 0) ? (angOut + 360) : angOut
			break
		case "rad":
			angOut = (angOut > 2*pi) ? (angOut - 2*pi) : angOut
			angOut = (angOut < 0) ? (angOut + 2*pi) : angOut
			break
	endswitch
	
	return angOut
End

//Returns angular statistics on the waves
Function AT_angularStats(suffix)

	String suffix
	
	//Finds the wave paths for analysis - this is just a string list of wave paths.
	String theWaveList = getWaveNames()
	Variable i,numWaves = ItemsInList(theWaveList,";")
	
	//Set data folder to that of the first wave on the wavelist
	SetDataFolder GetWavesDataFolder($StringFromList(0,theWaveList,";"),1)
		
	Make/O/N=(numWaves) $"vAng_" + suffix + "_median"
	Wave medianWave = $"vAng_" + suffix + "_median"
	
	Make/O/N=(numWaves) $"vAng_" + suffix + "_angDev"
	Wave angDevWave = $"vAng_" + suffix + "_angDev"
	
	Make/O/N=(numWaves) $"vAng_" + suffix + "_mean"
	Wave meanWave = $"vAng_" + suffix + "_mean"
	
	For(i=0;i<numWaves;i+=1)
		Wave theWave = $StringFromList(i,theWaveList,";")
		SetDataFolder GetWavesDataFolder(theWave,1)
		
		StatsCircularMoments/MODE=2 theWave
		Wave angStats = W_CircularStats
		
		medianWave[i] = angStats[11] * 180 / pi
		angDevWave[i] = angStats[10] * 180 / pi
		meanWave[i] = angStats[8] * 180 / pi
	EndFor
End

//AT_Correlate-------------------
//Cross correlates corresponding waves in each wave set. Normalizes the output waveform to the average auto-correlation of each input wave, as per MATLAB's xcorr function
//Used in FIGURE 3
Function AT_Correlate(inputWaveSet_1,inputWaveSet_2,suffix)
	String inputWaveSet_1,inputWaveSet_2,suffix
	
	String waveList1,waveList2
	
	waveList1 = getWaveNames(ignoreWaveGrouping=1,dataSet=inputWaveSet_1) //wave to wave analysis, ignore wave groupings
	waveList2 = getWaveNames(ignoreWaveGrouping=1,dataSet=inputWaveSet_2) //wave to wave analysis, ignore wave groupings
	
	Variable i,numWaves = ItemsInList(waveList1,";")
	For(i=0;i<numWaves;i+=1) //loop each wave
		Wave theWave_1 = $StringFromList(i,waveList1,";")
		Wave theWave_2 = $StringFromList(i,waveList2,";")
		
		If(!WaveExists(theWave_1) || !WaveExists(theWave_2))
			print "Couldn't find the wave for index: " + num2str(i)
			continue
		EndIf
		
		//duplicate waves for autocorrelations
		Duplicate/FREE theWave_1,auto_1 
		Duplicate/FREE theWave_2,auto_2
		
		SetDataFolder GetWavesDataFolder(theWave_1,1)
		
		//get autocorrelations
		Correlate/NODC theWave_1,auto_1
		Correlate/NODC theWave_2,auto_2
		
		//average autocorrelation peak
		Variable auto_avg = (WaveMax(auto_1) + WaveMax(auto_2)) / 2
		
		If(!strlen(suffix))
			Abort "Must supply a suffix for the output correlation wave"
		EndIf
		
		//Duplicate the second wave to make an output correlation wave
		String outWaveName = NameOfWave(theWave_1) + "_" + suffix
		Duplicate/O theWave_2,$outWaveName
		Wave corr = $outWaveName
		
		//correlate the two waves
		Correlate/NODC theWave_1,corr
		
		//scale the correlation wave to the average auto-correlation
		corr /= auto_avg		
	EndFor
End

