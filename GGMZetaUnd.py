#!/usr/bin/env python
import argparse, configparser, os, time, sys, datetime
import numpy as np
import pandas as pd # sudo apt install python3-pandas
import inspect, ahrs # pip install ahrs
import multiprocessing as mp
from typing import Optional
def cmdLineParse():
	parser = argparse.ArgumentParser(description='Estimate height anomalies (ZETA)/geoid undulations (UND) from EGM harmonic coefficients.')	
	parser.add_argument('-ggm', '--ggm'   , type=str, required=True , dest='ggm' , help='Path to EGM/GGM file.')	
	parser.add_argument('-z2n', '--z2n'   , type=str, required=False, dest='z2n' , help='Path to Zeta-to-N file.', default=None)
	parser.add_argument('-i'  , '--i'     , type=str, required=True , dest='ifl' , help='Path to input file of latitude/longitude in degrees. Format: txt.')
	parser.add_argument('-o'  , '--o'     , type=str, required=True , dest='ofl' , help='Path to output file of latitude/longitude in degrees and zeta/und in meters. Format: xlsx.')
	parser.add_argument('-c'  , '--cfg'   , type=str, required=True , dest='cfl' , help='Path to config file. Format: .cfg')
	parser.add_argument('-n'  , '--nprocs', type=int, required=False, dest='npr' , help='Number of CPU cores/threats used in multiprocessing. Leave default if multiprocessing is not used.', default=1)
	inputs = parser.parse_args()
	return inputs

def prgrsTime(disstr: str, prgrs_prct: float, prgrs_runtime: float, barlen: float=50, end: bool=False): # Create a progress bar to track the program
	prgrs_esttime = prgrs_runtime / prgrs_prct * (100 - prgrs_prct)	# the remaining running time in seconds
	prgrs_len = int(prgrs_prct/100*barlen) # the length of the job done so far.
	eta_len   = barlen - prgrs_len # the length of the remaining job.
	prgrs_txt = disstr.expandtabs(8) + " |" + "="*prgrs_len + " "*eta_len + "| " + str(int(round(prgrs_runtime))) + " s (" + str(round(float(prgrs_runtime)/60, 1)) + " mins) / " + str(int(round(prgrs_esttime))) + " s (" + str(round(float(prgrs_esttime/60), 1)) + " mins)"	
	print (prgrs_txt) if end else print(prgrs_txt, end="\r")	

class MPSLICE: # Class to work with slicing data for multi processing
	def slice_data_block (self, DATA: np.ndarray, nbins: int): # Split data into bins used in multiprocesing by blocks of arrays. This has the same output as slice_data_idxlst but quicker
		""" Split DATA into bins used in multiprocesing
		Args:
			DATA (ndarray)		: Array (N, M) needs to be split into bins.			
			nbins (int)			: Number (n) of bins by which Array is split.
		Returns:
			sliceDATA (list)	: List of n numpy arrays
		"""
		aver, res = divmod(len(DATA), nbins) # quotient and remainder
		nums = [0] + [(aver+1) if bin<res else aver for bin in range(nbins)] # list: number of subdatasets (rows of the array DATA) in each bin
		nums = [sum(nums[:ii+1]) for ii in range(len(nums))] # list: sum of numbers of subdatasets from the beginning to each bin
		sliceDATA = [DATA[ii:jj, :] for ii, jj in zip(nums[:len(nums)-1], nums[1:])] # list of arrays
		sliceDATA = [[item for item in row] for row in sliceDATA] # list (n) of list of arrays		
		return sliceDATA # list of list of array
	
class DATA:
	def TXTRead(self, txtfile: str, split: str='row', dismisrow: list=[], encoding: Optional[str]=None): # read a .txt file returning a list with each item being a column or row
		"""
		Parameters:
			txtfile (str)	: Input .txt file.
			split   (str)	: Return a list of 'row' or a list of 'column'.
			dismisrow (list): List of indices (zero-base) of dimissed rows.
			encoding: (str)	: encoding characters as bytes in the input file: options: None (default will be used), "utf8", ...
		Returns:
			rlst or clst (list):
		"""
		if encoding == None: file = open(txtfile)
		else:                file = open(txtfile, encoding=encoding) # file = open(txtfile, encoding="utf8")
		rlst = [line.split() for ii, line in enumerate(file) if ii not in dismisrow] # list of item with each being a row in the .txt file		
		file.close()
		if   split == 'row': return rlst
		elif split == 'col':
			clst = [[] for _ in range(len(max(rlst, key=len)))] # list of col, each is empty
			for col in range(len(clst)): clst[col] = [item[col] if col < len(item) else '' for item in rlst]
			return clst
	def HarCoefRead(self, HarCoefFile: str, dismisrow: list=[], keycol: int=-1, encoding: Optional[str]=None, zonalcoefcorr: bool=True, aELL: float=np.nan, fELL: float=np.nan, omgELL: float=np.nan, GMELL: float=np.nan, output: str='dict', verbose: bool=False, verbosetab: str=''): # Read spherical harmonic coefficients from file
		""" Read fully normalized spherical harmonic coefficients from, e.g., EGM2008
		Parameters:
			HarCoefFile (str)	: A file storing spherical harmonic coefficients.
			dismisrow: (list)	: List of indices of rows in the input file HarCoefFile that will be dismissed.
			keycol: (int)		: The index of the key column. If no key column then default -1.
			encoding: (str)		: encoding characters as bytes in the input file: options: None (default will be used), "utf8", ...
			zonalcoefcorr (bool): Correction of zonal coefficients C20, C40, C60 if True.
			output (str)		: Output as a list of list or dictionary.
		Returns:
			HarCoef (list/dict): List of spherical harmonic coefficients: (n, m, Cnm, Snm, sigmaCnm, sigmaSnm) or dictionary: {(n, m): {'Cnm': Cnm, 'Snm': Snm, 'CCnm': CCnm, 'CSnm': CSnm}}
									sigmaCnm, sigmaSnm, CCnm, CSnm are optional depending on the number of columns in the input file.
		References:
			[1] Rapp (1982) Report No.334: A FORTRAN program for the computation of gravimetric quantities from high degree spherical harmonic expansions
			[2] NGA.STND.0036_1.0.0_WGS84 https://nsgreg.nga.mil/doc/view?i=4085 => https://kubic-nsg-standards-nsgreg-nsgreg-files-6lxvt.s3.us-east-1.amazonaws.com/doc/Document/NGA%20STND%200036_1%200%200_WGS84.pdf?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAVXR7TTKDXD6K4UH6%2F20241206%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20241206T132450Z&X-Amz-Expires=900&X-Amz-SignedHeaders=host&response-cache-control=900&response-content-disposition=inline&response-content-type=application%2Fpdf&X-Amz-Signature=958f710b6496485b08854d00e8703fd6d6a5f7a67c3790f7a2d4957a107719b4
			[3] NIMA WGS1984 http://geodesy.unr.edu/hanspeterplag/library/geodesy/wgs84fin.pdf
			[4] Moritz, H. (2000). Geodetic reference system 1980. Journal of geodesy, 74(1), 128-133.
		"""		
		outputs = {'LIST', 'DICT', 'DICTIONARY'}
		if output.upper() not in outputs: raise ValueError("Unsupported outputs. Use 'LIST', 'DICT' or 'DICTIONARY'.")
		HarCoef = self.TXTRead(txtfile=HarCoefFile, split='row', dismisrow=dismisrow, encoding=encoding) # list of list of string
		if   output.upper() == 'LIST'                : HarCoef = [[int(subitem) if (ii < 2) else float(subitem.replace('D', 'E')) for ii, subitem in enumerate(item[keycol+1:])] for item in HarCoef] # list in [n, m, Cnm, Snm, sigmaCnm, sigmaSnm] => sigmaCnm, sigmaSnm are optional
		elif output.upper() in ['DICT', 'DICTIONARY']: HarCoef = {(int(item[keycol+1]), int(item[keycol+2])): {'Cnm': float(item[keycol+3].replace('D', 'E')), 'Snm': float(item[keycol+4].replace('D', 'E'))} for item in HarCoef} # in {(n, m): 'Cnm': Cnm, 'Snm': Snm} => sigmaCnm/sigmaSnm have not coded so far.
		bELL = aELL*(1-fELL)
		if zonalcoefcorr: # Correction for C20/C40/C60 using J2/J4/J6 => see [1] Eq.(4), Eq.(28-31)		
			m = (omgELL**2 * aELL**2 * bELL)/GMELL; f = fELL
			J2 =  2/3  * (f*(1-f/2)-m/2*(1-2*f/7+11/49*f**2))
			J4 = -4/35 * f*(1-f/2)*(7*f*(1-f/2)-5*m*(1-2*f/7))
			J6 =  4/21 * f**2 * (6*f-5*m)
			if   output.upper() == 'LIST':
				idx = [ii for ii, item in enumerate(HarCoef) if item[0] == 2 and item[1] == 0][0]; HarCoef[idx][2] = HarCoef[idx][2] + J2/np.sqrt(5)	# C20
				idx = [ii for ii, item in enumerate(HarCoef) if item[0] == 4 and item[1] == 0][0]; HarCoef[idx][2] = HarCoef[idx][2] + J4/np.sqrt(9)	# C40				
				idx = [ii for ii, item in enumerate(HarCoef) if item[0] == 6 and item[1] == 0][0]; HarCoef[idx][2] = HarCoef[idx][2] + J6/np.sqrt(13)	# C60; From [1] Eq.(31) Pg.(7) J6/np.sqrt(2*6+1)
			elif output.upper() in ['DICT', 'DICTIONARY']:
				HarCoef[(2, 0)]['Cnm'] += J2/np.sqrt(5)		# C20
				HarCoef[(4, 0)]['Cnm'] += J4/np.sqrt(9)		# C40
				HarCoef[(6, 0)]['Cnm'] += J6/np.sqrt(13)	# C60; From [1] Eq.(31) Pg.(7) J6/np.sqrt(2*6+1)
		return HarCoef # list or dict
	def toXLSX(self, data, dataheader = None, dataindex = False, datasheetname = None, xlsxfile: Optional[str]=None): # method to save data to a .xlsx file
		"""
		Parameters:
			data (ndarray/list/dictionary/pandas dataframe): data to be saved in the list/numpy array/dictionary format, e.g., data = [['val11'], ['21', 'val22'], ['31', '32', '33']]
			dataheader (list of string): list of data column headers
			dataindex (boolean): whether or not printing row indices
			datasheetname (str): sheet name
			xlsxfile (str): xls/xlsx file with full path
			NOTE: This can be used to save data to a .xls format, but a warning will be shown when opening the file, just ignore the warning!
		"""
		if xlsxfile == None: print(f'toosBasic/EXCEL/toXLSX: xlsxfile is empty.'); sys.exit()
		elif os.path.splitext(xlsxfile)[1] != '.xlsx': print(f'DATA().toXLSX: xlsxfile must be in the .xlsx format.'); sys.exit()
		if dataheader == None: dataheader = False
		if datasheetname == None: datasheetname = 'Sheet1'
		if os.path.exists(xlsxfile): writer = pd.ExcelWriter(xlsxfile, engine='openpyxl', mode='a', if_sheet_exists='replace') # can be 'replace', 'overlay', 'new', 'error'
		else: writer = pd.ExcelWriter(xlsxfile, engine='xlsxwriter')		
		pd.DataFrame(data).to_excel(writer, sheet_name=datasheetname, header=dataheader, index=dataindex); writer.close()

class GEOID:
	def __init__(self, aELL: float=np.nan, fELL: float=np.nan, aGGM: float=np.nan, GMGGM: float=np.nan): # Constructor: Initialise.
		"""
		Parameters:			
		returns:
		References:
			[1] NGA.STND.0036_1.0.0_WGS84 https://nsgreg.nga.mil/doc/view?i=4085 => https://kubic-nsg-standards-nsgreg-nsgreg-files-6lxvt.s3.us-east-1.amazonaws.com/doc/Document/NGA%20STND%200036_1%200%200_WGS84.pdf?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAVXR7TTKDXD6K4UH6%2F20241206%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20241206T132450Z&X-Amz-Expires=900&X-Amz-SignedHeaders=host&response-cache-control=900&response-content-disposition=inline&response-content-type=application%2Fpdf&X-Amz-Signature=958f710b6496485b08854d00e8703fd6d6a5f7a67c3790f7a2d4957a107719b4
			[2] NIMA WGS1984 http://geodesy.unr.edu/hanspeterplag/library/geodesy/wgs84fin.pdf
			[3] Moritz, H. (2000). Geodetic reference system 1980. Journal of geodesy, 74(1), 128-133.
		"""
		self.aELL = aELL; self.fELL = fELL; self.aGGM = aGGM; self.GMGGM = GMGGM
	def ZetaUndHarCoef_starmap2Async      (self, geoLatLonEHgt=np.array([]), HarCoef: dict={}, HarCoefZeta2N: dict={}, maxdeg: Optional[int]=None, ZeroDegZeta: float=0.0, output: str='UND', verbose: bool=False, verbosetab: str='', verboseprct: int=1, nprocs=mp.cpu_count()): # Call ZetaUndHarCoef with multiprocessing starmap2Async
		""" inputs & output are the same as those in ZetaUndHarCoef """
		if verbose: print(verbosetab + 'Number of CPUs used: ' + str(nprocs))
		if nprocs < 2: raise ValueError('\r\n' + verbosetab + __name__ + '.' + inspect.currentframe().f_code.co_name + ": The number of workers 'nprocs' must be at least 2 \n")
		if nprocs > len(geoLatLonEHgt): raise ValueError('\r\n' + verbosetab + __name__ + '.' + inspect.currentframe().f_code.co_name + ": The number of workers 'nprocs' must be smaller than or equal to " + str(len(geoLatLonEHgt)) + " \n")
		if verbose: sttm = datetime.datetime.now(); print(verbosetab + 'Separate input data into subdatasets for multiprocessor started at: ' + str(sttm))
		inp_lists = MPSLICE().slice_data_block(DATA=geoLatLonEHgt, nbins=nprocs) # list of list of array
		if verbose: fntm = datetime.datetime.now(); print (verbosetab + 'Finished at: ' + str(fntm) + '; Running time: ' + str(fntm - sttm))
		del geoLatLonEHgt # to save memory
		if verbose: sttm = datetime.datetime.now(); print(verbosetab + 'Compute height anomaly/geoid undulation started at: ' + str(sttm))
		with mp.Pool(nprocs) as pool: Outputs = pool.starmap_async(self.ZetaUndHarCoef, [(np.array(item), HarCoef, HarCoefZeta2N, maxdeg, ZeroDegZeta, output, verbose, verbosetab, verboseprct, ii+1) for ii, item in enumerate(inp_lists)]).get()	# this calls pool.starmap_async, normally faster than pool.starmap
		ZetaUnd  = np.array([item for sublist in Outputs for item in sublist[0]])
		CountPnm =      sum([item[1] for item in Outputs])
		return ZetaUnd, CountPnm # ndarray (N, )
	def ZetaUndHarCoef                    (self, geoLatLonEHgt=np.array([]), HarCoef: dict={}, HarCoefZeta2N: dict={}, maxdeg: int=None, ZeroDegZeta: float=0.0, output: str='UND', verbose: bool=False, verbosetab: str='', verboseprct: int=1, worker: int=1):
		"""Compute height anomaly (a.k.a quasigeoid height) or geoid undulation (a.k.a geoidal height) from the fully-normalized spherical harmonic coefficients
		Args:
			geoLatLonEHgt (ndarray): Array (N, 3) of geodetic Latitude/Longitude/Ellipsoidal Height. Lat/Lon in degrees, Ellipsoidal Heights in meters.
			HarCoef (dict)		   : Dictionary of fully-normalized spherical harmonic coefficients of the Earth's gravitational potential: {(n, m): 'Cnm': Cnm, 'Snm': Snm}.
			HarCoefZeta2N (dict)   : Dictionary of fully-normalized spherical harmonic coefficients of Zeta to N: {(n, m): 'Cnm': Cnm, 'Snm': Snm}.
			maxdeg (int)		   : Compute upto maxdeg. Optional: if None then the maximum degree from the fully normalized spherical harmonic coefficients HarCoef will be applied.
			ZeroDegZeta (float)	   : Zero degree height anomaly: Default -0.41 metres for WGS84-EGM2008.
			output (str)		   : Output of Height Anomalies ('ZETA') or geoid undulation ('UND' or 'UNDULATION').
		Returns:
			ZetaUnd (ndarray)	   : Array (N, ) of height anomaly (a.k.a quasigeoid heights) or geoid undulation in meters.
		References:
			[1] NGA.STND.0036_1.0.0_WGS84 https://nsgreg.nga.mil/doc/view?i=4085 => https://kubic-nsg-standards-nsgreg-nsgreg-files-6lxvt.s3.us-east-1.amazonaws.com/doc/Document/NGA%20STND%200036_1%200%200_WGS84.pdf?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAVXR7TTKDXD6K4UH6%2F20241206%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20241206T132450Z&X-Amz-Expires=900&X-Amz-SignedHeaders=host&response-cache-control=900&response-content-disposition=inline&response-content-type=application%2Fpdf&X-Amz-Signature=958f710b6496485b08854d00e8703fd6d6a5f7a67c3790f7a2d4957a107719b4
			[2] NIMA WGS1984 http://geodesy.unr.edu/hanspeterplag/library/geodesy/wgs84fin.pdf
			[3] my paper: 2014-01-HUMG JMES-Computation of Height Anomaly in Vietnam territory ...-Bui et al.pdf
			[4] my paper: 2014-02-VJES-Refining Geoidal Heights Computed by Fully-normalised Harmonic Coefficients...-Bui et al
			[5] Hirt (2012): https://doi.org/10.1007/s00190-012-0550-y
			[6] Rapp (1997): https://doi.org/10.1007/s001900050096
			[7] "Description of Files Related to Using the EGM2008 Global Gravitational Model to Compute Geoid Undulations with Respect to WGS 84"
			NOTE: Pnm uses the spherical latitudes but normal gravity (gamma) uses the geodetic latitudes.
			NOTE: Pnm depends on geocentric spherical latitude but not longitude. So, geoLatLonEHgt is sorted in the ascending order of latitude before computation.
					This is advantageous in that Pnm of a point (i) will not be computed but use Pnm of the point (i-1) if the two points (i) & (i-1) are of the same latitude.
					This works normally in case of computing for gridded data in latitude/longitude. The results then need to be converted back to the original order.
		"""
		outputs = {'ZETA', 'UND', 'UNDULATION'} # 'UND' and 'UNDULATION' are the same.
		if output.upper() not in outputs: raise ValueError("modules.py/GEOID/ZetaUndHarCoef: Unsupported outputs. Use 'ZETA', 'UND' or 'UNDULATION'.")
		geoLatLonEHgt, sortidx = ARRAY().Sort(array=geoLatLonEHgt, sortcol=0) # sort the points in the ascending order of latitude
		ZetaUnd = np.zeros(len(geoLatLonEHgt)) # np.aray (N, ), Height Anomalies (Zeta) or Geoid Undulation (UND) depending on output		
		maxdegcoef = max([item[0] for item in list(HarCoef.keys())]); maxdeg = int(maxdegcoef) if maxdeg == None or maxdeg > maxdegcoef else maxdeg
		sLatLonRad = CoordConvTrans().geoLatLonEHgt2SphCoordConv(geoLatLonEHgt=geoLatLonEHgt, aELL=self.aELL, fELL=self.fELL, output='sLat') # np.array (N, 3): geocentric spherical latitude in degrees, longitude in degrees, geocentric radius
		tic, intvlprct = time.time(), max(1, round(float(len(geoLatLonEHgt))*verboseprct/100)); countPnm = 0
		for ii, item in enumerate(sLatLonRad):
			sLatDeg, LonRad, sRadius = item[0], np.radians(item[1]), item[2] # geocentric spherical latitude in degrees, longitude in radians, geocentric spherical radius in meters
			gLatDeg, LonDeg, EHgt    = geoLatLonEHgt[ii, 0], geoLatLonEHgt[ii, 1], geoLatLonEHgt[ii, 2]
			if ii == 0 or (sLatLonRad[ii, 0] != sLatLonRad[ii-1, 0]): Pnm = self.FullNormAssocLegendreFns(sLatDeg=sLatDeg, maxdeg=maxdeg); countPnm += 1 # Dictionary {(n, m): Pnm} computed from Geocentric Spherical Latitude. See [1] Eq.(6-4) Pg.(6-3)
			gamma_h = self.normal_gravity(geodlatdeg=gLatDeg, londeg=LonDeg, eHgt=EHgt, method='ahrs', higheHgt=False) # Normal gravity at point with ellipsoid h			
			for n in range(2, maxdeg+1): # degree
				temp = 0
				for m in range(n+1): temp += (HarCoef[(n, m)]['Cnm'] * np.cos(m*LonRad) + HarCoef[(n, m)]['Snm'] * np.sin(m*LonRad)) * Pnm[(n, m)] # order
				ZetaUnd[ii] += (self.aGGM/sRadius)**n * temp # ZetaUnd[ii] += (self.a/sRadius)**n * temp
			ZetaUnd[ii] *= self.GMGGM / (gamma_h * sRadius)
			ZetaUnd[ii] += ZeroDegZeta # Zero degree height anomaly			
			if output.upper()[:3] == 'UND': # compute zeta2N from Zeta-to-N_to2160_egm2008
				for n, m in HarCoefZeta2N: ZetaUnd[ii] += (HarCoefZeta2N[(n, m)]['Cnm'] * np.cos(m*LonRad) + HarCoefZeta2N[(n, m)]['Snm'] * np.sin(m*LonRad)) * Pnm[(n, m)] # see [7] Eq.(4) p.(1)
			if verbose and ((ii + 1) == 1 or (ii + 1) % intvlprct == 0 or (ii + 1) == len(sLatLonRad)):
				prgrs_prct, prgrs_runtime = float(ii + 1) / len(sLatLonRad) * 100, time.time() - tic
				prgrs_str = verbosetab + 'Worker: ' + str(worker) + '; Estimated point #: ' + str(ii + 1) + '/' + str(len(sLatLonRad)) + ' [' + str(int(round(prgrs_prct))) + '%]'
				prgrsTime(prgrs_str, prgrs_prct, prgrs_runtime, end=True) if (ii + 1) == len(sLatLonRad) else prgrsTime(prgrs_str, prgrs_prct, prgrs_runtime, end=False)
		ZetaUnd = ARRAY().Reverse(array=ZetaUnd, sortid=sortidx) # reverse the results ZetaUnd to the original order of the coordinate array geoLatLonEHgt
		return ZetaUnd, countPnm # ndarray (N, )'''	
	def FullNormAssocLegendreFns(self, sLatDeg: float, maxdeg: int): # Compute Fully Normalised Associated Legedre Functions providing Geocentric Spherical Latitude (Phi)
		"""Compute Fully Normalized Associated Legendre Functions up to maxdeg.
			Args:
				sLatDeg (float)			: Geocentric SPHERICAL latitude  in degrees.
				maxdeg (int)   			: Compute up to maximum degree maxdeg.
			Returns:
				PnmSinsLat (dictionary)	: Dictionary of Fully Normalized Associated Legendre Functions {(n, m): PnmSinsLat}.
			References:
				[1] my paper: 2014-01-HUMG JMES-Computation of Height Anomaly in Vietnam territory ...-Bui et al.pdf
				[2] my paper: 2014-02-VJES-Refining Geoidal Heights Computed by Fully-normalised Harmonic Coefficients...-Bui et al
				[3] Holmes & Featherstone (2002) https://doi.org/10.1007/s00190-002-0216-2
		"""
		sLatRad = np.radians(sLatDeg) # geocentric spherical latitude in radians
		PnmSinsLat = {(0, 0): 1, (1, 0): np.sqrt(3) * np.sin(sLatRad), (1, 1): np.sqrt(3) * np.cos(sLatRad)} # Dictionary; see [1] pg.(79), [2], [3] pg.(282)
		for n in range(2, maxdeg+1): # degree
			for m in range(n+1): # order
				if n == m : cnm = np.sqrt((2*m+1)/(2*m));								    PnmSinsLat[(n, m)] = np.cos(sLatRad) * cnm * PnmSinsLat[(n-1, m-1)] # see [1] Eq.(2), [2] Eq.(3), [3] Eq.(13)
				else      : anm = np.sqrt(((2*n-1)*(2*n+1))/((n-m)*(n+m)));				    PnmSinsLat[(n, m)] = anm * np.sin(sLatRad) * PnmSinsLat[(n-1,   m)] # see [1] Eq.(2), [2] Eq.(3), [3] Eq.(13)
				if n > m+1:	bnm = np.sqrt(((2*n+1)*(n+m-1)*(n-m-1))/((n-m)*(n+m)*(2*n-3))); PnmSinsLat[(n, m)] -= bnm * PnmSinsLat[(n-2, m)]					# see [1] Eq.(2), [2] Eq.(3), [3] Eq.(13)
		return PnmSinsLat # Dictionary {(n, m): PnmSinsLat}'''
	def normal_gravity(self, geodlatdeg: float, londeg: float, eHgt: float, method: str='ahrs', higheHgt: bool=True):
		"""Compute the normal gravity of a point on/above the WGS84 ellipsoid at a given latitude.
			Args:
				geodlatdeg (float): Geodetic latitude  in degrees.
				londeg (float)	  : Geodetic longitude in degrees.
				eHgt (float)	  : Ellipsoidal height in meters.
				method (str)	  : use the AHRS library or equations to compute.
				higheHgt (boolen) : use the equation working with high/low ellipsoidal height/altitude. Note: If higheHgt is selected then based on equations only, because AHRS has not applied these equations.
			Returns:
				gamma_h (float)	  : Normal gravity at the point P (latitude/longitude/eHeight) in m/s^2.
			References:
				[1] NGA.STND.0036_1.0.0_WGS84 https://nsgreg.nga.mil/doc/view?i=4085 => https://kubic-nsg-standards-nsgreg-nsgreg-files-6lxvt.s3.us-east-1.amazonaws.com/doc/Document/NGA%20STND%200036_1%200%200_WGS84.pdf?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAVXR7TTKDXD6K4UH6%2F20241206%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Date=20241206T132450Z&X-Amz-Expires=900&X-Amz-SignedHeaders=host&response-cache-control=900&response-content-disposition=inline&response-content-type=application%2Fpdf&X-Amz-Signature=958f710b6496485b08854d00e8703fd6d6a5f7a67c3790f7a2d4957a107719b4
				[2] NIMA WGS1984 http://geodesy.unr.edu/hanspeterplag/library/geodesy/wgs84fin.pdf
				[3] AHRS library: https://ahrs.readthedocs.io/en/latest/wgs84.html#
				[4] Moritz, H. (2000). Geodetic reference system 1980. Journal of geodesy, 74(1), 128-133.
		"""
		if (not higheHgt or eHgt == 0) and method.lower() == 'ahrs': gamma_h = ahrs.utils.WGS().normal_gravity(lat=geodlatdeg, h=eHgt) # use the AHRS library. see [3].
		else: # compute with equations instead of using the AHRS library. See [1], [2]
			latrad = np.radians(geodlatdeg) # to radians			
			if eHgt == 0 or (eHgt != 0 and not higheHgt): # on the Ellipsoid or Above/Below the ellipsoid with low ellipsoidal height: use the equation of low ellipsoidal height/altitude
				if self.ellipsoid.upper() == 'WGS84':
					gamma_e = 9.7803253359  # Gravity on the ellipsoid WGS84 at the equator (m/s^2), in [1], [2] Tab. 3.6, p. 3-10.
					gamma_p = 9.8321849379  # Gravity on the ellipsoid WGS84 at the pole    (m/s^2), in [1], [2] Tab. 3.6, p. 3-10.
				# Compute normal gravity on the ellipsoid (i.e., eHgt = 0) first
				gamma_h = (self.a*gamma_e*np.cos(latrad)**2 + self.b*gamma_p*np.sin(latrad)**2)/np.sqrt(self.a**2*np.cos(latrad)**2 + self.b**2*np.sin(latrad)**2) # (1): m/s2; see the Eq in [3].
				if eHgt != 0: m = self.omega**2 * self.a**2 * self.b / self.GM; gamma_h *= (1 - 2/self.a*(1+self.f+m-2*self.f*np.sin(latrad)**2)*eHgt + 3/self.a**2*eHgt**2) # see [1] Eq.(4-3), p.4-3; [2] Eq.(4-3), p.4-2; [3]
			elif eHgt != 0 and higheHgt: # Above/Below the ellipsoid with high ellipsoidal height: use the equation of high ellipsoidal height/altitude
				# Below use Eqs (4-4) to (4-15) in [1] (p. 4-3 to 4-5), [2]
				XYZ = CoordConvTrans().geoLatLonEHgt2XYZConv(np.array([[geodlatdeg, londeg, eHgt]])); x, y, z = XYZ[0,0], XYZ[0,1], XYZ[0,2]
				E = np.sqrt(self.a**2 - self.b**2)
				tmp = x**2 + y**2 + z**2 - E**2; u = np.sqrt(1/2*tmp*(1 + np.sqrt(1 + 4*E**2*z**2/tmp**2)))
				beta = np.arctan((z * np.sqrt(u**2 + E**2))/(u * np.sqrt(x**2 + y**2)))
				w = np.sqrt((u**2 + E**2*np.sin(beta)**2)/(u**2 + E**2))
				q = 1/2*((1 + 3*u**2/E**2)*np.arctan(E/u) - 3*u/E)
				q0 = 1/2*((1+3*self.b**2/E**2)*np.arctan(E/self.b)-3*self.b/E)				
				q_prime = 3*(1 + u**2/E**2)*(1 - u/E*np.arctan(E/u)) - 1
				gamma_beta = (self.omega**2*self.a**2/np.sqrt(u**2 + E**2)*q/q0 - self.omega**2*np.sqrt(u**2+E**2)) * (1/w * np.sin(beta) * np.cos(beta))
				gamma_u    = -1/w * (self.GM/(u**2+E**2) + self.omega**2*self.a**2*E/(u**2+E**2)*q_prime/q0*(1/2*np.sin(beta)**2-1/6) - self.omega**2*u*np.cos(beta)**2)
				gamma_h    = np.sqrt(gamma_beta**2 + gamma_u**2)
		return gamma_h # scalar in m/s2

class ARRAY: # Class to work with an ARRAY with methods, e.g., sort, reverse
	def Sort   (self, array: np.array, sortcol: int): # sort a 2D array by sortcol column returning the sorted array and sorted indices
		""" Sort rows of a 2D array in the ascending of a column.
		Parameters:
			array (ndarray)	   : a 2D array (N, M) of which the rows to be sorted.
			sortcol (int)	   : The index of the column based on which the array is sorted.
		returns:
			sortArray (ndarray): The sorted array (N, M).
			sortid (ndarray)   : An array of sorted indices (N, ).
		"""
		sortid = np.argsort(array[:, sortcol]) # indices in the ascending order of the sortcol column
		sortArray = array[sortid] # reorder the array by the sorted indices
		return sortArray, sortid
	def Reverse(self, array: np.array, sortid): # reverse an array to the original indices providing sortid
		""" Reverse an array back to the original indices by sortid. This is usually done after apply the above function Sort.
		Parameters:
			array (ndarray)	   : a 2D array (N, M).
			sortid (ndarray)   : An array of sorted indices (N, ).
		returns:
			revsArray (ndarray): The reversed array (N, M).
		"""
		revsArray = np.empty_like(array)
		revsArray[sortid] = array
		return revsArray

class CoordConvTrans: # Class: Coordinate conversion and transformation with methods:
	def geoLatLonEHgt2SphCoordConv(self, geoLatLonEHgt: np.ndarray, aELL: float, fELL: float, output: str='sLat'): # Conversion of geodetic coordinates to spherical coordinates
		"""
		Parameters:
			geoLatLonEHgt (ndarray): Array (N, 3) of geodetic Latitude/Longitude/Ellipsoidal Height values. Lat/Lon in degrees.
			output (str)		   : estimate spherical latitudes if 'slat', otherwise colatitudes will be computed
		Returns:
			LatLonRad (ndarray)	   : Array (N, 3) of Colatitude (theta) or spherical latitude (Phi), Longitude (lambda), Geocentric Radius (r) values.
										coLat is the angle from the Z axis to the line connecting the geocenter and the point. 
										sLat is the spherical latitude, which can be computed as: sLat = pi/2-coLat.
		"""
		XYZ = self.geoLatLonEHgt2XYZConv(LatLonEHgt=geoLatLonEHgt, aEll=aELL, fEll=fELL) # np.array (N, 3)
		r = np.sqrt(XYZ[:, 0]**2 + XYZ[:, 1]**2 + XYZ[:, 2]**2)  # np.array (N, ): Geocentric radius: r = sqrt(X^2+Y^2+Z^2)
		if output.upper() == 'SLAT': sLat = np.rad2deg(np.arcsin(XYZ[:, 2]/r))  # np.array (N, ); spherical latitudes in degrees: sLat = asin(Z/r) = atan(Z/sqrt(X^2+Y^2)) = 90 - coLat
		else: coLat = np.rad2deg(np.arccos(XYZ[:, 2]/r)) # np.array (N, ); colatitudes in degrees: coLat = acos(Z/r) = atan(sqrt(X^2+Y^2)/Z) = 90 - sLat
		return np.column_stack((sLat, geoLatLonEHgt[:, 1], r)) if output.upper() == 'SLAT' else np.column_stack((coLat, geoLatLonEHgt[:, 1], r)) # np.array (N, 3)
	def geoLatLonEHgt2XYZConv     (self, LatLonEHgt: np.ndarray, aEll: float, fEll: float): # Conversion of geodetic coordinates Lat/Lon/EllHgt to cartesian coordinates X/Y/Z
		"""
		Parameters:
			LatLonEHgt (ndarray): Array (N, 3) of Latitude/Longitude/Ellipsoidal Height values. Lat/Lon in degrees
		Returns:
			XYZ (ndarray)		: Array (N, 3) of cartesian coordinates X/Y/Z
		"""
		e2 = 2*fEll-fEll**2
		N = aEll/np.sqrt((1-e2*np.sin(np.radians(LatLonEHgt[:, 0]))**2)) # np.array (N, )
		X = (N+LatLonEHgt[:, 2])*np.cos(np.radians(LatLonEHgt[:, 0]))*np.cos(np.radians(LatLonEHgt[:, 1])) # np.array (N, )
		Y = (N+LatLonEHgt[:, 2])*np.cos(np.radians(LatLonEHgt[:, 0]))*np.sin(np.radians(LatLonEHgt[:, 1])) # np.array (N, )
		Z = ((1-e2)*N+LatLonEHgt[:, 2])*np.sin(np.radians(LatLonEHgt[:, 0])) # np.array (N, )'''
		return np.column_stack((X, Y, Z)) # np.array (N, 3)

if __name__ == '__main__':	
	print ()
	print ("##############################################################################################################################")
	print ("#                                                                                                                            #")
	print ("#    This is a tool used to estimate height anomalies (ZETA)/geoid undulations (UND) from EGM/GGM harmonic coefficients      #")
	print ("#                                                                                                                            #")
	print ("##############################################################################################################################")	

	sttm = datetime.datetime.now()
	# command line parser
	inputs  = cmdLineParse()
	Z2Nfile = inputs.z2n if inputs.z2n == None else os.path.abspath(inputs.z2n)
	inpfile, oupfile, cfgfile, GGMfile = os.path.abspath(inputs.ifl), os.path.abspath(inputs.ofl), os.path.abspath(inputs.cfl), os.path.abspath(inputs.ggm)
	nprocs = inputs.npr
	if nprocs > mp.cpu_count(): raise ValueError('\r\n' + __name__ + '.' + inspect.currentframe().f_code.co_name + ": The number of workers 'nprocs' must be smaller than or equal to " + str(mp.cpu_count()) + " \n")
	cwd     = os.getcwd()

	# Create a parser and read the config file
	config = configparser.ConfigParser(inline_comment_prefixes=('#', ';')); config.read(cfgfile)	
	skrwsINP = [ii for ii in range(config.getint('INPUT', 'skiprows', fallback=0))]
	aELL   = config.getfloat('ELLIPSOID', 'a', fallback=6378137.0)
	invf   = config.getfloat('ELLIPSOID', 'invf', fallback=298.257223563); fELL  = 1/invf
	GMELL  = float(config['ELLIPSOID']['GM'])
	omgELL = float(config['ELLIPSOID']['omega'])
	aGGM  = config.getfloat('GGM', 'a', fallback=6378137.0)
	GMGGM = float(config['GGM']['GM'])
	skrwsGGM    = [ii for ii in range(config.getint('GGM', 'skiprows', fallback=0))]
	keycolGGM   = config.getint('GGM', 'keycol', fallback=-1)
	encodingGGM = config.get('GGM', 'encoding', fallback=None)
	degGGM      = config.getint('CONSTANTS', 'deg', fallback=float('inf'))
	print(f'-- Equatorial radius of the reference ellipsoid [m]: {aELL}')
	print(f' - Flattening of the reference ellipsoid: 1/{invf}')
	print(f' - Geocentric Gravitational Constant of the reference ellipsoid [m3/s2]: {GMELL}')
	print(f' - Nominal Mean Angular Velocity of the Earth of the reference ellipsoid [rads/s]: {omgELL}')
	print(f'\n-- The scaling factor of the GGM model a [m]: {aGGM}')
	print(f' - The scaling factor of the GGM model GM [m3/s2]: {GMGGM}')

	print(f'\n-- Read coordinates Lat/Lon from the file            : {inpfile}')	
	LatLonEHgt = np.array([[float(item[0]), float(item[1]), 0] for item in DATA().TXTRead(txtfile=inpfile, split='row', dismisrow=skrwsINP)])
	print(f' - Read spherical harmonic coefficients from the file: {GGMfile}');	HarCoefDict = DATA().HarCoefRead(HarCoefFile=GGMfile, dismisrow=skrwsGGM, keycol=keycolGGM, encoding=encodingGGM, zonalcoefcorr=True, aELL=aELL, fELL=fELL, omgELL=omgELL, GMELL=GMELL, output='dict', verbose=True, verbosetab='\t-- ')
	if Z2Nfile != None:
		print(f' - Read spherical harmonic coefficients of Zeta2N from the file: {Z2Nfile}');
		HarCoefZeta2NDict, output = DATA().HarCoefRead(HarCoefFile=Z2Nfile, dismisrow=[], keycol=-1, zonalcoefcorr=False, output='dict'), 'UND'
	else:
		HarCoefZeta2NDict, output = {}, 'ZETA'
	degGGM = min(degGGM, max([item[0] for item in list(HarCoefDict.keys())]))	
	print(f'\n-- Zeta/Undulation will be computed up to degree: {degGGM}')
	print(f' - Compute Height Anomalies (Zeta)/Geoid Undulation (UND).')	
	if nprocs > 1: ZetaUnd, _ = GEOID(aELL=aELL, fELL=fELL, aGGM=aGGM, GMGGM=GMGGM).ZetaUndHarCoef_starmap2Async(geoLatLonEHgt=LatLonEHgt, HarCoef=HarCoefDict, HarCoefZeta2N=HarCoefZeta2NDict, maxdeg=degGGM, ZeroDegZeta=0, output=output, verbose=True, nprocs=nprocs)
	else: ZetaUnd, _ = GEOID(aELL=aELL, fELL=fELL, aGGM=aGGM, GMGGM=GMGGM).ZetaUndHarCoef              (geoLatLonEHgt=LatLonEHgt, HarCoef=HarCoefDict, HarCoefZeta2N=HarCoefZeta2NDict, maxdeg=degGGM, ZeroDegZeta=0, output=output, verbose=True)
	xlsxData = np.column_stack((LatLonEHgt[:, :2], ZetaUnd))
	dataheader = ['LAT (deg)', 'LON (deg)', 'ZETA (m)'] if output == 'ZETA' else ['LAT (deg)', 'LON (deg)', 'UND (m)']
	DATA().toXLSX(data=xlsxData, datasheetname='Sheet1', dataheader=dataheader, dataindex=True, xlsxfile=oupfile)
	print(f' - Results have been saved to the file: {oupfile}')
	print('\n-- Computation finished at: ', datetime.datetime.now())
	fntm = datetime.datetime.now()
	print (' - Total running time: ' + str(fntm - sttm))
	print ('-----------------------------------------------------------------------------------------')
	print ('=========================================================================================')