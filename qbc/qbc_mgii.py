from astropy.io import fits
from astropy.constants import c, G
from plots import plot_results_show, plot_overlap, plot_gc, plot_sum
from astropy.table import QTable
from astropy.table import Column
from astropy.cosmology import Planck15 as cosmo
import astropy.coordinates as coord
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import astropy.units as u
import numpy as np
import argparse
import pickle
import os

running_mode ='run'#default
# running_mode='spec_p'
#running_mode='test'
# running_mode='neil'

help_text = 'Looks for Qso Absorbers/Galaxy Clusters pairs and calculates dn/dz'
sign_off = 'Author: Hector Salas <hsalas@das.uchile.cl>'

parser = argparse.ArgumentParser(description=help_text, epilog=sign_off)

parser.add_argument('-minew', '--min_EW', type=float, default=0.6, dest='minEW', metavar='MINEW', help='Float. Minimum rest equivalent width of the MgII 2796 AA line to be considered. Default value 0.6 AA.',action='store')
parser.add_argument('-maxew', '--max_EW', type=float, default=2.0, dest='maxEW', metavar='MAXEW', help='Float. Maximum rest equivalent width of the MgII 2796 AA line to be considered. Default value 2.0 AA.',action='store')
parser.add_argument('-minmass', '--min_mass', type=float, default=13.6, dest='minmass', metavar='MINM', help='Float. Minimim cluster mass to be considered in units of log10(Msun). Default value 13.6',action='store')
parser.add_argument('-maxmass', '--max_mass', type=float, default=20.0, dest='maxmass', metavar='MAXM', help='Float. Maximum cluster mass to be considered in units of log10(Msun). Default value 20.0',action='store')
parser.add_argument('-minz', '--min_z', type=float, default=0, dest='minz', metavar='MINZ', help='Float. Minimum redshift to be considered. Default value 0.',action='store')
parser.add_argument('-max', '--max_z', type=float, default=10, dest='maxz', metavar='MAXZ', help='Float. Maximum redshift to be considered. Default value 10.',action='store')
parser.add_argument('-minip', '--min_IP', type=float, default=0.1, dest='minIP', metavar='MINIP', help='Float. Minimum impact parameter to be considered. Default value 0.1',action='store')
parser.add_argument('-maxip', '--max_IP', type=float, default=40, dest='maxIP', metavar='MAXIP', help='Float. Maximum impact parameter to be considered. Default value 40.',action='store')
parser.add_argument('-pr', '--plot_result', type=str, default='no', dest='plotresult', metavar='PR', help='''str. Whether to display a plot of the reuslts. Default 'no'.''', action='store')
parser.add_argument('-s', '--significance', type=float, default=3, dest='significance', metavar='S', help='Float. Significance of the min equivalent width',action='store')
parser.add_argument('-s2n', '--signal2noise', type=str, default='local', dest='signaltonoise', metavar='S2N', help="str. Signal to noise to be use, 'local' or 'global'. Default 'local'.",action='store')
parser.add_argument('-dvqso', '--min_dv_qso', type=float, default=12000.0, dest='mindvqso', metavar='DVQSO', help='Float. Minimum velocity difference  (in km/s) between the cluster and the Qso at the Qso redshift. Default value 12000 km/s',action='store')
parser.add_argument('-ls', '--lambda_survey', type=float, default=1215.67, dest='lsurvey', metavar='LS', help='Float. Lower wavelength  limit, standar choices are CIV  1548.0 AA and LyAlpha 1215.67 AA. Default value Lya 1215.67 AA',action='store')
parser.add_argument('-dls', '--delta_lambda_survey', type=float, default=6000.0, dest='dlsurvey', metavar='DLS', help='Float. Minimum distance blueward of ls (in km/s) to start looking for MGII absorbers. Default value 6000 km/s',action='store')
parser.add_argument('-zct', '--z_cl_type', type=str, default='all', dest='zcltype', metavar='ZCT', help='''str. Selected redshift type: 'spec', 'phot' or 'all'.''',action='store')
parser.add_argument('-d', '--distance', type=str, default='com', dest='d', metavar='D', help='str. Type of distance to to be used: comoving(com), proper(pro) or R_200(r200). Default value com',action='store')
parser.add_argument('-gr', '--grid_type', type=str, default='log', dest='gridtype', metavar='GR', help='''str. Grid type to be used, linear 'lin', logarithmic 'log' , same number hits per bin 'snh', or same number pairs per bin 'snp'.''', action='store')
parser.add_argument('-n', '--number_bins', type=float, default=1.0, dest='nbins', metavar='NB', help='float. Scale factor for the number of bins to use in the grid',action='store')
parser.add_argument('-m', '--mode', type=int, default=1, dest='mode', metavar='M', help='''int. Mode in which the program wil be run, with the following options. 
					0: Plot Redshift Distribution of MgII absorbres and Galaxy clusters
					1: Run dn/dz vs Impact Parameter (default) 
					2: Run dn/dz vs Equivalent Width for r< R200 & r>R200.
					3: Run dn/dz vs Redshift
					4: Run dn/(dzdw) vs Equivalent Width (in development).
					5, 6 : Run from saved files (in  development) 
''',action='store')
arguments = parser.parse_args()
min_dv_qso = arguments.mindvqso # in km/s
dls = arguments.dlsurvey
ls = arguments.lsurvey
min_mass = arguments.minmass # in log10(Msun)
max_mass = arguments.maxmass # in log10(Msun)
min_EW = arguments.minEW # in AA
max_EW = arguments.maxEW # in AA
min_z = arguments.minz 
max_z = arguments.maxz 
min_IP = arguments.minIP # in units of dd (Mpc comoving, Mpc proper or R200)
max_IP = arguments.maxIP # in units of dd (Mpc comoving, Mpc proper or R200)
z_cl_type = arguments.zcltype
s = arguments.significance
sn = arguments.signaltonoise
n = arguments.nbins
grt = arguments.gridtype
m = arguments.mode
dd = arguments.d # distance units
plot_result =arguments.plotresult

delta_v = 500.0*(u.km/u.s)
l_0 = 2796.3542699 # rest wavelength mgii in AA
lmin_sdss = 3800.0 
lmax_sdss =	9200.0
z_min = lmin_sdss/l_0-1.0
z_max = lmax_sdss/l_0-1.0

dataset = 'mge10e'+str(min_mass)
min_mass = (10**min_mass)*u.Msun
max_mass =  (10**max_mass)*u.Msun
min_EW = min_EW*u.AA
max_EW = max_EW*u.AA

alias_dict={
1548.0 : 'lim_civ',
1215.67 : 'lim_lya',
'ew' :  'rew_mgii_2796',
'ew_pares' : 'w_min_s1',
'mass' : 'cluster_mass',
'com' : 'sep_comoving',
'pro' : 'sep_proper',
'r200' : 'sep_200'
}

# ----- Functions -----
def create_dir(dir_name):
	'''Check if the directory dir_name exists, if it dosen't it creates it
	'''
	if os.path.isdir(dir_name):
		pass
	else:
		os.makedirs(dir_name)	

def names(alias):
	'''Gives the value for 'alias' 'according to 'alias_dict'  defined above'''
	try:
		name = alias_dict[alias]
		return name
	except KeyError:
		raise KeyError('''select a valid tag (-d) 'com' (comoving distance), 'pro' (proper distance) or '200' (r200 distance) or mass (cluster mass)''')

def dztodv(dz,z):
	'''Transform a redshift difference (dz) to a velocity difference at redshift z, the result is given in Km/s
	
	inputs: 
		dz:	Redshift difference (float).
		z:	Redshift (float).
	returns:
		dv:	Velocity difference in km/s (quantity).
	'''	
	dv = dz*c/(1.0+z)
	dv = dv.to(u.km/u.s)
	return(dv)

def dvtodz(dv,z):
	'''Takes a difference in velocity in km/s and transform it to a difference in redshift

	inputs: 
		dv:	Velocity difference in Km/s (quantity).
		z:	Redshift (float).
	returns:
		dz: 	Redshift difference (float).
	'''
	dv = dv.to(u.m/u.s)
	dz = dv*(1.0+z)/c
	dz = dz.to(u.dimensionless_unscaled)
	return(dz.value)

def w_min(sn, z, s):
	'''Gives the minimum resolution width as a function of redshift (z), signal to noise ratio (sn), Sloan resolution (sdss) and significance value (s)
	
	inputs: 
		sn: 	Signal to noise ratio (float).
		z: 	Redshift (float).
		s: 	Significance (float).
	returns:
		w: 	minimum resolution width (float).
	'''
	w = s*l_obs(l_0,z)/(sn*r_sdss(z))
	return(w)

def l_obs(l0,z):
	'''Gives the observed wavelength for l0 at redshift z.

	inputs: 
		l0: 	Rest wavelength (float).
		z: 	Redshift (float).
	returns:
		l: 	observed wavelength (float).
	'''
	l = l0*(1.0+z)
	return(l)

def r_sdss(z):
	'''Gives the Sloan resolution for the observed MgII wavelength (rest W. 2796 AA), based on a linear model with a resolution of 1500 at 3800 AA and 2500 at 9200 AA
	
	inputs: 
		z: 	Redshift (float).
	returns:
		y: 	Sloan resolution of MgII_2796 at redshift z (float). 
	'''
	x = l_obs(l_0,z)
	y = (2500.0-1500.0)/(lmax_sdss-lmin_sdss)*(x-lmin_sdss)+1500.0
	return(y)

def bins_edges(tabla, tipo, name, begin, end, n_bins, n):
	'''Calls to the corresponding functions to create the different bins edges.

	input: 
		tabla: 	Table with the data (QTable).
		tipo: 	Type of binning to be used (str).
		name:	Name of property to use in the binnig (str). Eg. sep comoving, sep proper, mass.
		begin: 	Starting point of the binning (float).
		end: 	Final point of the bining (float).
		n_bins: 	Number of bins (int).
		n: 	'scale' factor for n_bins, n_bins will be multiplied by this number (float).

	output:
		edges: Array with the bin's edges (numpy.ndarray).
	'''
	n = n*n_bins
	n = int(n)
	if tipo =='snh':
		edges = bin_same_number(hits_phot, name, begin, end, n)
	elif tipo == 'snp':
		edges = bin_same_number(pares_phot, name, begin, end, n)
	elif (tipo == 'lin' or tipo == 'log') or tipo =='neil':
		edges = bins(tipo, begin, end, n)
	elif tipo == 'lopez':
		edges = [0.05, 0.15,0.3,0.6, 1,2,3]
	elif tipo == 'one':
		edges = [begin, end]
	else:
		raise ValueError('''Choose a valid option for the input: 'lin'  or 'log' , 'snh'', 'snp'.''')
	return(edges)

def bins(tipo, begin, end, n):
	'''Creates an array with the edges of the bins to be used.

	input: 
		tipo: 	Type of bins to be used, linear 'lin', logarithmic 'log' (str).
		begin:	Initial point of the binning (float)
		end: 	Final point of the binning (float).
		n: 	Number of bins to use (int).
	returns:
		grid: 	Array with the bin's edges (numpy.ndarray).
	'''
	if tipo == 'lin':
		grid = np.linspace(begin, end, n+1)
	elif tipo == 'log':
		grid = np.logspace(np.log10(begin), np.log10(end), n+1)
	elif tipo == 'neil':
		grid = np.logspace(-1, np.log10(11.5), 13)		
	return(grid)

def bin_same_number(tabla, name, begin, end, n):
	''' Creates bins with the same number of elements

	input:
		table: 	Table with the data (QTable).
		name: 	Name of property to use to  eg. sep comoving, sep proper, mass (str)
		begin: 	Starting point of the binning (float).
		end: 	Final point of the bining (float).
		n: 	Number of bins (int).
	
	output:
		bins: 	Array with the bin's edges (numpy.ndarray).
	'''
	if name =='mass':
		if 'mass' in tabla.keys():
			tag =  'mass'
	else:
		tag = names(name)
	aux = tabla[tag]
	aux = np.asarray(aux)
	aux.sort()
	cond = aux < end
	aux = aux[cond]
	npt = len(aux)
	bins = np.interp(np.linspace(begin, npt, n + 1), np.arange(npt), aux)
	return(bins)

def err_plus(n):
	'''Gives an estimate of the upper confidence limit for n.

	input: 
		n:	Observed value (float).
	returns:
		ul:	Upper confidence limit (float).
	'''
	ucl = np.sqrt(n+3./4.)
	return(ucl)

def err_minus(n):
	'''Gives an estimate of the lower confidence limit for n.

	input: 
		n:	Observed value (float).
	returns:
		lcl:	 Lower confidence limit (float).
	'''
	if n == 0:
		lcl = 0
	else:
		lcl = np.sqrt(n-1./4.)
	return(lcl)

def max_angular_dist(d, clusters, qso):
	''' returns all cluster Qso pairs with an angular distance lower than d.

	inputs:
		d: 	Angular distance (quantity).
		cluster: 	Galaxy clusters catalog (Table).
		qso: 	Qso/MgII absorbers catalog (Table).
	returns:
		id_cl: 	Index of the pairs in the clusters table (int).
		id_mgii:	Index of the pairs in the qso table (int).
		sep2d: 	Angular separation of the pairs (Angle).
	'''
	mgii_coord = coord.SkyCoord(qso['ra'], qso['dec'], unit='deg')
	cluster_coord = coord.SkyCoord(clusters['ra'], clusters['dec'], unit='deg')
	id_mgii, id_cl, sep2d, _ = cluster_coord.search_around_sky(mgii_coord, d)
	return(id_cl, id_mgii, sep2d)

def cuts_pares(s, z_min, z_max, min_EW, min_mass, max_mass, min_z, max_z, min_IP, max_IP, dist,  tabla):      
	''' Takes a table with with cluster/qso pairs and returns only the pairs that obey the inputs limits.

	inputs:	
	 	s :	significance of MgII detection (float).
		z_min:	Minimum redshift redshift where detection is possible (float). Given by survey limits
		z_max:	Maximum redshift redshift where detection is possible (float). Given by survey limits
		min_EW:	Minimum equivalent width to consider in AA (Quantity).
		min_mass:	Minimum cluster mass to consider in solar mass units (Quantity).
		max_mass:	Maximum cluster mass to consider in solar mass units (Quantity).
		min_z:	Minimum redshift redshift to consider (float). Chosen by user.
		max_z:	Maximum redshift redshift to consider (float). Chosen by user.
		min_IP:	Minimum impact parameter to be considered in units of dist (float).
		max_IP:	Maximum Impact parameter to be considered in units of dist(float).
		dist:	Units of distance for impact parameter (str). 
				Acepted values:	sep_comoving, sep_proper or sep_200.
		tabla:	Table with all the clusters/qso pairs (QTable).
	returns:
		tabla: Table with the clusters/qso pairs (QTable).
		tabla_rejected: Table with the rejected clusters/qso pairs (QTable).
	'''
	print('Deleting  pairs that do not fit the criteria')
	tabla_rejected =  tabla['id_cluster', 'id_qso']
	col = Column( name='rejected', dtype = np.dtype((str, 25)), length=len(tabla))
	tabla_rejected.add_column(col)
	tabla_rejected['rejected'] = 'no'

	#removes pairs with impact parameter > max_IP in units of dist
	cond = tabla[dist] > max_IP
	tabla_rejected['rejected'][cond] = 'IP > max_IP'
	#removes pairs with impact parameter < min_IP in units of dist
	cond = tabla[dist] < min_IP
	tabla_rejected['rejected'][cond] = 'IP < min_IP'
	#removes pairs with z_cluster > than z_qso
	cond = tabla['z_cluster'] >= tabla['z_qso']
	tabla_rejected['rejected'][cond] = 'z_cluster >= Z_qso'
	#removes pairs with mass < min_mass
	cond = tabla['cluster_mass'] * u.Msun < min_mass
	tabla_rejected['rejected'][cond] = 'mass < min_mass '
	#removes pairs with mass > max_mass
	cond = tabla['cluster_mass'] * u.Msun >  max_mass
	tabla_rejected['rejected'][cond] = 'mass > max_mass '
	#deletes pairs  with redshift below z_min, related to overall redshift limits
	cond = tabla['z_cluster'] < z_min
	tabla_rejected['rejected'][cond] = 'z < z_min '
	#deletes pairs  with redshift above z_max, related to overall redshift limits
	cond = tabla['z_cluster'] > z_max
	tabla_rejected['rejected'][cond] = 'z > z_max '
	#deletes paris with redshift below min_z, related to cuts in redshift
	cond = tabla['z_cluster'] < min_z
	tabla_rejected['rejected'][cond] = 'z < min_z '
	#deletes pairs  with redshift above z_max, related to cuts in redshift
	cond = tabla['z_cluster'] > max_z
	tabla_rejected['rejected'][cond] = 'z > max_z '

	cond = tabla['z_cluster'] > tabla['z_max_qso']
	tabla_rejected['rejected'][cond] = 'z_cluster > z_max_qso'
	if ls == 1548.0:
		cond = tabla['z_cluster'] < tabla['z_min_qso_civ']
		tabla_rejected['rejected'][cond] = 'z_cluster < Z_min_qso'
	elif ls == 1215.67:
		cond = tabla['z_cluster'] < tabla['z_min_qso_lya']
		tabla_rejected['rejected'][cond] = 'z_cluster < Z_min_qso'
	else:
		cond = tabla['z_cluster'] < tabla['z_min_qso']
		tabla_rejected['rejected'][cond] = 'z_cluster < Z_min_qso'
	if sn == 'global':
		cond = tabla['w_min_s1']*s*u.AA >  min_EW  #NT: me da la impresion que esta condicion esta mala; favor revisar
	elif sn == 'local':
		cond = tabla['w_min_s1_local']*s*u.AA >  min_EW  #NT: me da la impresion que esta condicion esta mala; favor revisar
	tabla_rejected['rejected'][cond] = 'w_min > min_EW'	
	cond =	tabla_rejected['rejected'] == 'no'
	tabla = tabla[cond]			
	cond =	tabla_rejected['rejected'] != 'no'
	tabla_rejected = tabla_rejected[cond]

	# add corresponing units to columns
	tabla['sep_angular'].unit = u.deg
	tabla['ra_cluster'].unit = u.deg
	tabla['dec_cluster'].unit = u.deg
	tabla['ra_qso'].unit = u.deg
	tabla['dec_qso'].unit = u.deg
	tabla['cluster_mass'].unit = u.Msun
	tabla['dv_cluster'].unit = u.km/u.s	
	tabla['w_min_s1'].unit = u.AA
	tabla['sep_comoving'].unit = u.Mpc
	tabla['sep_proper'].unit = u.Mpc
	return(tabla, tabla_rejected)

def cuts_hits(s, z_min, z_max, min_EW, max_EW, min_mass, max_mass, min_z, max_z, min_IP, max_IP, dist,  tabla):
	'''Takes a table with  absorbers/cluster pairs and returns only the hits that obey the inputs limits.
	 
	 inputs:
	 	s :	significance of MgII detection (float).
		cluster: Galaxy clusters catalog (Table).
		z_min:	Minimum redshift redshift where detection is possible (float). Given by survey limits
		z_max:	Maximum redshift redshift where detection is possible (float). Given by survey limits
		min_EW:	Minimum equivalent width to consider in AA (Quantity).
		max_EW:	Maximum equivalent width to consider in AA (Quantity).
		min_mass:	Minimum cluster mass to consider in solar mass units (Quantity).
		max_mass:	Maximum cluster mass to consider in solar mass units (Quantity).
		min_z:	Minimum redshift redshift to consider (float). Chosen by user.
		max_z:	Maximum redshift redshift to consider (float). Chosen by user.
		min_IP:	Minimum impact parameter to be considered (float).
		max_IP:	Maximum Impact parameter to be considered (float).
		dist:	Units of distance for impact parameter (str). 
				Acepted values:	sep_comoving, sep_proper or sep_200.
		tabla:	Table with the clusters/qso pairs (QTable).

	returns:
		tabla: Table with the MgII absorbers/qso pairs (QTable).
		tabla_rejected: Table with the rejected MgII absorbers/qso pairs (QTable).
	'''
	print('Deleting  hits that do not fit the criteria')
	tabla_rejected =  tabla['id_cluster', 'id_mgii']
	col = Column( name='rejected', dtype = np.dtype((str, 25)), length=len(tabla))
	tabla_rejected.add_column(col)
	tabla_rejected['rejected'] = 'no'

	#removes pairs with impact parameter > max_IP in units of dist
	cond = tabla[dist] > max_IP
	tabla_rejected['rejected'][cond] = 'IP > max_IP'
	#removes pairs with impact parameter < min_IP in units of dist
	cond = tabla[dist] < min_IP
	tabla_rejected['rejected'][cond] = 'IP < min_IP'
	#removes pairs with z_cluster > than z_qso
	cond = tabla['z_cluster'] >= tabla['z_qso']
	tabla_rejected['rejected'][cond] = 'z_cluster >= Z_qso'
	#removes pairs with ew < min_EW
	cond =tabla['rew_mgii_2796']*u.AA < min_EW
	tabla_rejected['rejected'][cond] = 'rew_mgii_2796 < min_EW'
	#removes pairs with ew > max_EW
	cond = tabla['rew_mgii_2796'] *u.AA> max_EW
	tabla_rejected['rejected'][cond] = 'rew_mgii_2796 > max_EW'
	#removes pairs  with mass < min_mass
	cond = tabla['cluster_mass'] * u.Msun < min_mass
	tabla_rejected['rejected'][cond] = 'mass < min_mass '
	#removes pairs  with mass > max_mass
	cond = tabla['cluster_mass'] * u.Msun >  max_mass
	tabla_rejected['rejected'][cond] = 'mass > max_mass '	
	#deletes pairs   with redshift below z_min
	cond = tabla['z_cluster'] < z_min
	tabla_rejected['rejected'][cond] = 'z < z_min '
	#deletes pairs   with redshift above z_max
	cond = tabla['z_cluster'] > z_max
	tabla_rejected['rejected'][cond] = 'z > z_max '
	#deletes paris with redshift below min_z, related to cuts in redshift
	cond = tabla['z_cluster'] < min_z
	tabla_rejected['rejected'][cond] = 'z < min_z '
	#deletes pairs  with redshift above z_max, related to cuts in redshift
	cond = tabla['z_cluster'] > max_z
	tabla_rejected['rejected'][cond] = 'z > max_z '

	cond = tabla['z_cluster'] > tabla['z_max_qso']
	tabla_rejected['rejected'][cond] = 'z_cluster > z_max_qso'
	if ls == 1548.0:
		cond = tabla['z_cluster'] < tabla['z_min_qso_civ']
		tabla_rejected['rejected'][cond] = 'z_cluster < Z_min_qso'
	elif ls == 1215.67:
		cond = tabla['z_cluster'] < tabla['z_min_qso_lya']
		tabla_rejected['rejected'][cond] = 'z_cluster < Z_min_qso'
	else:
		cond = tabla['z_cluster'] < tabla['z_min_qso']
		tabla_rejected['rejected'][cond] = 'z_cluster < Z_min_qso'
	if sn == 'global':
		cond = tabla['w_min_s1']*s*u.AA >  min_EW  #NT: me da la impresion que esta condicion esta mala; favor revisar
	elif sn == 'local':
		cond = tabla['w_min_s1_local']*s*u.AA >  min_EW  #NT: me da la impresion que esta condicion esta mala; favor revisar
	tabla_rejected['rejected'][cond] = 'w_min > min_EW'			
	cond = abs(tabla['dv']) > tabla['dv_cluster']
	tabla_rejected['rejected'][cond] = 'dv > dv_cluster'
	cond =	tabla_rejected['rejected'] == 'no'
	tabla = tabla[cond]			
	cond =	tabla_rejected['rejected'] != 'no'
	tabla_rejected = tabla_rejected[cond]

	# add corresponing units to columns
	tabla['sep_angular'].unit = u.deg
	tabla['ra_cluster'].unit = u.deg
	tabla['dec_cluster'].unit = u.deg
	tabla['ra_qso'].unit = u.deg
	tabla['dec_qso'].unit = u.deg
	tabla['dv_cluster'].unit = u.km/u.s
	tabla['cluster_mass'].unit = u.Msun
	tabla['dv'].unit = u.km/u.s
	tabla['w_min_s1'].unit = u.AA
	tabla['rew_mgi_2853'].unit = u.AA
	tabla['rew_mgii_2803'].unit = u.AA
	tabla['rew_mgii_2796'].unit = u.AA
	tabla['rew_feii_2600'].unit = u.AA
	tabla['rew_feii_2374'].unit = u.AA
	tabla['rew_feii_2344'].unit = u.AA
	tabla['rew_feii_2586'].unit = u.AA
	tabla['rew_feii_2383'].unit = u.AA
	tabla['err_rew_mgi_2853'].unit = u.AA
	tabla['err_rew_mgii_2803'].unit = u.AA
	tabla['err_rew_mgii_2796'].unit = u.AA
	tabla['err_rew_feii_2600'].unit = u.AA
	tabla['err_rew_feii_2374'].unit = u.AA
	tabla['err_rew_feii_2344'].unit = u.AA
	tabla['err_rew_feii_2586'].unit = u.AA
	tabla['err_rew_feii_2383'].unit = u.AA
	tabla['sep_comoving'].unit = u.Mpc
	tabla['sep_proper'].unit = u.Mpc
	return(tabla, tabla_rejected)	

def grilla(tabla, grid, name):
	''' Separates the data contained in the table 'tabla' into the bins defined by 'grid', the results are sored in a dictionary where each key 
	correspond to the right edge of the corresponding bin.

	input:
		tabla:	Table with clusters/Qso pairs or with cluster/absorbers hits (QTable).
		grid: 	array with the edges of the bins to be used (numpy.ndarray).
		name: 	Name of property used in the grid (str). Either one of the three disctances: 'pro', 'com', 'r200' , or ew.
	returns:
		grid_result:	Dictionary with the table data separated by the bins defined by grid (dict).
	'''
	grid_result = {}
	tag = names(name) 
	if tag == 'sep_200':
		for i in range(len(grid)-1):
			aux = tabla.copy()
			cond = aux[tag] < grid[(-1*i)-1]
			aux = aux[cond]
			cond = aux[tag] >= grid[(-1*i)-2]
			aux = aux[cond]
			grid_result[grid[(-1*i) - 1]] = aux

	elif tag == 'w_min_s1':#revisar con mas atencion esta parte
		for i in range(len(grid)-1):
			aux = tabla.copy()
			cond = aux[tag]*s < grid[(-1*i)-2]*u.AA
			aux = aux[cond]
			grid_result[grid[(-1*i) - 1]*u.AA] = aux
	
	elif tag == 'rew_mgii_2796' :
		for i in range(len(grid)-1):
			aux = tabla.copy()
			cond = aux[tag] < grid[(-1*i)-1]*u.AA
			aux = aux[cond]
			cond = aux[tag] >= grid[(-1*i)-2]*u.AA
			aux = aux[cond]
			grid_result[grid[(-1*i) - 1]*u.AA] = aux
	
	else:
		for i in range(len(grid)-1):
			aux = tabla.copy()
			cond = aux[tag] < grid[(-1*i)-1]*u.Mpc
			aux = aux[cond]
			cond = aux[tag] >= grid[(-1*i)-2]*u.Mpc
			aux = aux[cond]
			grid_result[grid[(-1*i) - 1]*u.Mpc] = aux
	return(grid_result)
	
def g_c(tabla, zmin, zmax, delta):
	'''Calculates the 'cluster redshift path density' between zmin and zmax.

	inputs:	
		tabla:	Table with clusters/Qso pairs (QTable).
		zmin: 	Minimum redshift (float).
		zmax:	Maximum redshift(float).
		delta:	Step to be used in the redshift path calculation (float).
	returns:
		z:	Array of redshifts starting at zmin and ending at zmax with a step delta (numpy.ndarray).
		g:	Array containing the value of the redshift path density for each value of z (numpy.ndarray).
	'''
	z =  np.linspace(zmin, zmax, int((zmax-zmin)/delta), endpoint=True)
	g = np.zeros(len(z))

	if len(tabla) == 0:
		return(g)
	for row in tabla:
		aux =	np.zeros(len(z)) + 1.
		cond = z < row['z_cluster'] - row['dz_cluster'] 
		aux[cond] = 0
		cond = row['z_cluster'] + row['dz_cluster'] < z
		aux[cond] = 0
		g = g + aux
	return(z, g) 

def z_path(tabla, zmin, zmax):
	'''integration of the redshift path density in the redshift interval [z_min,z_max]
		
	inputs:	
		tabla:	Table with clusters/Qso pairs (QTable).
		zmin:  	Minimum redshift (float).
		zmax: 	Maximum redshift (float).
	returns:
		aux2: 	Redshift path (float).
		gc: 	Table with the data of the redshift path density(QTable).
	'''
	if len(tabla) == 0:
		gc = QTable(names=['z','gc'])
		return(0, gc)
	d = np.max(tabla['dz_cluster'])
	delta = 1.e-4
	aux1 = tabla.copy()
	aux1.sort('z_cluster')
	z, g = g_c(tabla, zmin, zmax, delta)
	gc = QTable([z,g], names=['z','gc'])
	aux2 =  np.sum(g)*delta
	return(aux2,gc)	

def redshift_path(target, zmin, zmax):
	'''Gives the redshift path for different values of the impact parameter.
	
	input:
		target:	Dictionary with tables for the different values of the impact parameter (dict).
		zmin:	Minimum redshift (float)
		zmax: 	Maximum redshift (float)
	returns:
		aux: 	Dictionary with the redshift path for each table in the target dictionary (dict).
		aux2: 	Dictionary with the redshift path density for each table in the target dictionary (dict).
	'''
	print('Obtaining Redshift path')	
	a = target.keys()
	a = list(a)
	a.sort()
	aux = {}
	aux2 = {}
	contador = 0
	for i in range(len(a)):	
		print('{}/{}'.format(contador, len(a)))		
		dz_c, gc = z_path(target[a[i]], zmin, zmax)
		contador += 1
		aux[a[i]] =  dz_c
		if type(a) == u.quantity.Quantity:
			aux2[a[i].value] = gc
		else:
			aux2[a[i]] = gc
	return(aux, aux2)

def results_table(x, target1,target2,target3):
	'''creates a QTable with the results
	
	inputs:	
		target1:	Dictionary with number of hits (dict).
		target2:	Dictionary with redshift path (dict).
		target3: 	Dictionary with number of pairs (dict).
	returns:
		results: 	Table summarizing the results stored in the dictionaries 'target1', 'target2' and 'target3' (QTable). 		
	'''
	a = target1.keys()
	a = list(a)
	a.sort()
	aux = {}
	contador = 0
	results = QTable()
	if (x == 'com' or x == 'pro') or x =='r200':	 
		bi = Column(np.zeros(len(a)), name='b_i')
		bf = Column(np.zeros(len(a)), name='b_f')
	elif x == 'ew':
		bi = Column(np.zeros(len(a)), name='ew_i')
		bf = Column(np.zeros(len(a)), name='ew_f')
	dn = Column(np.zeros(len(a)), name='dn')
	sigma_n_mas = Column(np.zeros(len(a)), name='sigma+_n')
	sigma_n_menos = Column(np.zeros(len(a)), name='sigma-_n')
	dz_col = Column(np.zeros(len(a)), name='dz')
	pares = Column(np.zeros(len(a)), name='#pares')
	dndz = Column(np.zeros(len(a)), name='dn/dz')
	dndz_err_mas = Column(np.zeros(len(a)), name='sigma+_dn/dz')
	dndz_err_menos = Column(np.zeros(len(a)), name='sigma-_dn/dz')	
	results.add_column(bi, index=0)
	results.add_column(bf, index=1)
	results.add_column(dn, index=2)
	results.add_column(sigma_n_mas, index=3)
	results.add_column(sigma_n_menos, index=4)
	results.add_column(dz_col)
	results.add_column(pares)
	results.add_column(dndz)
	results.add_column(dndz_err_mas)
	results.add_column(dndz_err_menos)
	for i in range(len(a)):	
		dn = len(target1[a[i]])
		dz = target2[a[i]]
		pares = len(target3[a[i]])
		if x =='r200':
			results['b_i'][i] = a[i - 1]
			results['b_f'][i] = a[i]
		elif x == 'com' or x == 'pro':	
			results['b_i'][i] = a[i - 1].value
			results['b_f'][i] = a[i].value
		elif x == 'ew':
			results['ew_i'][i] = a[i - 1].value
			results['ew_f'][i] = a[i].value
		results['dn'][i] = dn
		results['sigma+_n'][i] = err_plus(dn)
		results['sigma-_n'][i] = err_minus(dn)
		results['dz'][i] = dz
		results['#pares'][i] = pares
		if dz != 0:
			results['dn/dz'][i] = dn/dz
			results['sigma+_dn/dz'][i] = 1.0/dz*err_plus(dn)
			results['sigma-_dn/dz'][i] = 1.0/dz*err_minus(dn) 
		if results['dn'][i] > results['#pares'][i]:
			raise ValueError('This should not happen.')	
	if (x == 'com' or x == 'pro') or x =='r200':	 
		results['b_i'][0] = min_IP# pensar una mejor forma de hacer esto
	elif x=='ew':
		results['ew_i'][0] = min_EW.value
	return(results)

def cluster_table(cluster, z_min, z_max, min_mass, max_mass):
	'''Creates a QTable for the galaxy cluster data

	inputs:
		cluster: 
		z_min:	Lower redshift limit (float).
		z_max:	Upper redshift limit (float).
		min_mass:	Lower mass limit in solar mass units (Quantity).
		max_mass:	Upper mass limit in solar mass units (Quantity).
	returns:
		cluster:	Table with the selected cluster data (QTable).
		cluster_spec: 	Table with the selected cluster data, only spectroscopic redshift (QTable).
		cluster_phot: 	Table with the selected cluster data, only photometric redshift (QTable).
		cluster_rejected:	Table with the rejected cluster data (QTable).
		cluster_0: 	Table with all the cluster data (QTable).
	'''
	cluster = QTable(cluster, names = [name.lower() for name in cluster.names])
	cluster['z_lambda'].name = 'z_phot'
	try:#if 0:
		cluster['z_lambda_err'].name = 'z_phot_err'
		cluster['lambda'].name = 'richness'
		cluster['lambda_err'].name = 'richness_err'
	except KeyError:#else:
		cluster['z_lambda_e'].name = 'z_phot_err'
		cluster['lambda_chisq'].name = 'richness'
		cluster['lambda_chisq_e'].name = 'richness_err'
		cluster['bcg_spec_z'].name = 'z_spec'

	#''' add a column whit z=z_spec or z=z_phot when z_spec is not available'''
	cluster['z'] = cluster['z_spec']
	cond = cluster['z'] == -1
	cluster['z'][cond] = cluster['z_phot'][cond]
	
	#'''add mass estimate for redmapper'''
	mass = (10.**14. / (cosmo.H0.value/70.) ) * np.exp(1.72, dtype='f4') * np.power(cluster['richness']/60.,1.08, dtype='f4')
	cluster['mass'] = mass*u.Msun

	#'''Adds a column with r200 for each cluster'''
	rho_200 = 200. * cosmo.critical_density(cluster['z']).to('Msun/Mpc3') #in Msun/Mpc3	
	mass_200 = cluster['mass']# in Msun		
	r_200 = np.power(mass_200/rho_200/(4.*np.pi/3.), 1./3, dtype='f4')#.value #in Mpc
	cluster['r200_mpc'] = r_200
	sigma_v =np.sqrt(G*mass_200/r_200)

	# Adds a column with a velocity difference (dv) in the cluster table.
	# If the cluster have a spectroscopic redshift dv is obtained from rho_200 and mass_200 ussing the virial eq.
	# Otherwise the error in the photometric redshift is used
	cluster['delta_v'] = np.zeros(len(cluster), dtype='f4') + sigma_v
	cond = cluster['z_spec'] == -1
	cluster['delta_v'][cond] = dztodv(cluster['z_phot_err'][cond], cluster['z_phot'][cond])

	#add a column with dz
	cluster['dz'] = dvtodz(cluster['delta_v'],cluster['z'])

	col = Column( name='rejected', dtype = np.dtype((str,25)), length = len(cluster))
	cluster.add_column(col)
	cluster['rejected'] = 'no'

	#deletes clusters below z_min and above z_max
	cond = cluster['z'] < z_min
	cluster['rejected'][cond] = 'z < z_min '
	cond = cluster['z'] > z_max
	cluster['rejected'][cond] = 'z > z_max '

	#removes clusters with mass < min_mass
	cond = cluster['mass'] < min_mass
	cluster['rejected'][cond] = 'mass < min_mass '

	#removes clusters with mass > max_mass
	cond = cluster['mass'] >  max_mass
	cluster['rejected'][cond] = 'mass > max_mass '

	# import pdb; pdb.set_trace()
	cluster_0 = cluster.copy()
	cond = cluster['rejected'] != 'no'
	cluster_rejected = cluster[cond]
	cond = cluster['rejected'] == 'no'
	cluster = cluster[cond]
	#selects only the spectroscopic or photometric clusters
	cond = cluster['z_spec'] == -1
	cluster_phot = cluster[cond]
	cond = cluster['z_spec'] != -1
	cluster_spec = cluster[cond]	

	return(cluster, cluster_spec, cluster_phot, cluster_rejected, cluster_0)

def mgii_table(mgii, mgii_search, min_EW, max_EW, z_max_cluster, ls):
	'''Creates a QTable for the mgii absorbers and the data

		mgii:
		mgii_search:
		min_EW:	Lower rest frame equivalent width limit (float).
		max_EW: 	Upper rest frame equivalent width limit (float).
		z_max_cluster:  Maximum galaxy clusters redshit plus delta z (float). z_max_cluster= p.max(cumulo['z']+0.2).
	returns:
		mgii:	Table with the selected MgII absorbers data (QTable)
		mgii_search:	Table with the Qso searched for MgII absorbers data(QTable).
		mgii_rejected: 	Table with the rejected MgII absorbers data  (QTable).
		mgii_0: 	Table with all the MgII absorbers data (QTable)	.
	'''
	mgii = QTable(mgii,names=[name.lower() for name in mgii.names])

	mgii_search = QTable(mgii_search, names = [name.lower() for name in mgii_search.names])

	#'''Add the mean of the s/n ratio to the searched qso# 	
	sn_medio = np.mean(mgii['spec_snr_median'], dtype='f4')
	mgii['index_qso_hs'] = np.zeros(len(mgii), dtype='i4') - 1
	mgii_search['spec_snr_median'] = np.zeros(len(mgii_search), dtype='f4') + sn_medio
	mgii_search['index_qso'] = np.zeros(len(mgii_search), dtype='i4') - 1 
	mgii_search['index_qso_hs'] = np.arange(len(mgii_search), dtype='i4')
	mgii_search['id'] = np.arange(len(mgii_search), dtype='i4')
	match = max_angular_dist(3*u.arcsec, mgii_search, mgii)
	match = QTable(match, names=['id_search','id_qso','sep_ang'])
	index_search = mgii_search['id'][match['id_search']]
	index_mgii = match['id_qso']
	mgii_search['spec_snr_median'][index_search] = mgii['spec_snr_median'][index_mgii]
	mgii_search['index_qso'][index_search] = mgii['index_qso'][index_mgii]
	mgii['index_qso_hs'][index_mgii] = mgii_search['index_qso_hs'][index_search]
	mgii_search.remove_column('id')

	#adds a column with zmax
	mgii['zmax'] = np.minimum(l_0*(1.0 + mgii['zqso'] - dvtodz(min_dv_qso*u.km/u.s,mgii['zqso'])),lmax_sdss)/l_0 - 1.0
	mgii_search['zmax'] = np.minimum(l_0*(1.0 + mgii_search['zqso'] - dvtodz(min_dv_qso*u.km/u.s,mgii_search['zqso'])),lmax_sdss)/l_0 - 1.0

	#adds a column with zmin
	mgii['zmin'] = np.maximum(ls * (1.0 + mgii['zqso'] + dvtodz(dls*u.km/u.s, mgii['zqso'])),lmin_sdss)/l_0 - 1.0
	mgii_search['zmin'] = np.maximum(ls * (1.0 + mgii_search['zqso'] + dvtodz(dls*u.km/u.s,mgii_search['zqso'])),lmin_sdss)/l_0 - 1.0

	col = Column( name='rejected', dtype = np.dtype((str, 25)), length=len(mgii))
	mgii.add_column(col)
	mgii['rejected'] = 'no'

	cond = mgii['zabs'] > z_max_cluster+0.2
	mgii['rejected'][cond] = 'zabs > z_max_cluster + 0.2'

	#removes aborbers with ew < min_EW
	cond = mgii['rew_mgii_2796']*u.AA < min_EW
	mgii['rejected'][cond] = 'rew_mgii_2796 < min_EW'
	#removes aborbers with ew > max_EW
	cond = mgii['rew_mgii_2796']*u.AA > max_EW
	mgii['rejected'][cond] = 'rew_mgii_2796 > max_EW'

	mgii_0 = mgii.copy()
	cond = mgii['rejected'] != 'no'
	mgii_rejected = mgii[cond]
	cond = mgii['rejected'] == 'no'
	mgii = mgii[cond]

	return(mgii, mgii_search, mgii_rejected, mgii_0)

def load_catalogs(ls, mode='run'):
	'''Loads the Galaxy Clusters, MgII absorption and Qso catalogs

	add info about the used catalogs
	'''
	if mode == 'run':
		#catalogs to use for normal runing of the code
		cluster = fits.getdata('../data/clusters/dr8_run_redmapper_v5.10_lgt5_catalog.fit')
		if ls >= 1548.0:
			mgii = fits.getdata('../data/mgii/Trimmed_SDSS_DR7_107.fits')
		else:
			mgii = fits.getdata('../data/mgii/Expanded_SDSS_DR7_107.fits')
		mgii_search = fits.getdata('../data/mgii/Expanded_Searched_Quasars_SDSS_DR7_107.fits')

	elif mode =='test':
		#test catalogs with 100 random pairs
		cluster = fits.getdata('../../../../QbC/test_cluster.fits')
		mgii = fits.getdata('../../../../QbC/test_mgii.fits')
		mgii_search = fits.getdata('../../../../QbC/test_qso.fits')
		
	elif mode == 'test_table_0':
		#test catalog created with a kwon behavior
		cluster = fits.getdata('./test_tables/cumulo_0.fits')
		mgii = fits.getdata('./test_tables/mgii_0_pro.fits')
		mgii_search = fits.getdata('./test_tables/qso_0_pro.fits')

	elif mode == 'neil':
		#to compare with neil
		cluster = fits.getdata('../data/clusters/dr8_run_redmapper_v5.10_lgt5_catalog.fit')
		mgii = fits.getdata('../data/mgii/Expanded_SDSS_DR7_107.fits')
		mgii_search = fits.getdata('../data/mgii/QSObased_Expanded_SDSS_DR7_107.fits')
	return(cluster, mgii, mgii_search)

def load_master_tables(z_cl_type, mode):
	'''Loads the master tables

	input:
		z_cl_type:	Clustrer's redshift type, 'spec' or 'phot' (str)
		mode:	Type of tables to use, this will define whow the program is run. 
				mode='test' will load test table resulting  in a test running mode.
				mode='neil' will load tables with inputs similar to those used in neil's version of the code, resulting in a mode compare with those results
				mode='spec_p' will load a version of the spectroscopic clusters table in which their photometric values are used instead of the spectroscopic. photometric clusters will not be affected
				Anything else will result in loading the standard tables, resulting in a normalrun of the program
	output
		master_table_pares_xxxx master table containnig the pairs. xxxx=z_cl_type used 
		master_table_hits_xxxx master table containg the hits. xxxx=z_cl_type used
	'''
		
	if z_cl_type == 'phot':

		if  mode == 'test':
			#for testing
			master_table_hits_phot = fits.getdata( '../saved_files/tabla_maestra_hits_phot_test.fits')
			master_table_pares_phot = fits.getdata( '../saved_files/tabla_maestra_pares_phot_test.fits')
		elif  mode == 'neil':
			#to reproduce neil
			master_table_hits_phot = fits.getdata( '../saved_files/tabla_maestra_hits_phot_neil.fits')
			master_table_pares_phot = fits.getdata( '../saved_files/tabla_maestra_pares_phot_neil.fits')	
		else:
			#for normal running
			master_table_hits_phot = fits.getdata( '../saved_files/tabla_maestra_hits_phot.fits')
			master_table_pares_phot = fits.getdata( '../saved_files/tabla_maestra_pares_phot.fits')

		master_table_hits_phot = QTable(master_table_hits_phot)
		master_table_pares_phot = QTable(master_table_pares_phot)
		
		cond = master_table_hits_phot['id_qso_hs']!=-1
		master_table_hits_phot = master_table_hits_phot[cond]	

		if z_cl_type == 'phot':
			return(master_table_hits_phot, master_table_pares_phot)

	elif z_cl_type == 'spec':
		if  mode == 'test':
			#for testing
			master_table_hits_spec = fits.getdata( '../saved_files/tabla_maestra_hits_spec_test.fits')
			master_table_pares_spec = fits.getdata( '../saved_files/tabla_maestra_pares_spec_test.fits') 
		elif mode =='neil':
			#to reproduce neil
			master_table_hits_spec = fits.getdata( '../saved_files/tabla_maestra_hits_spec_neil.fits')
			master_table_pares_spec = fits.getdata( '../saved_files/tabla_maestra_pares_spec_neil.fits')
		elif mode =='spec_p':
			#to use the photometometric valus of the spectroscopic clusters for comparison
			master_table_hits_spec = fits.getdata( '../saved_files/tabla_maestra_hits_spec_p.fits')
			master_table_pares_spec = fits.getdata( '../saved_files/tabla_maestra_pares_spec_p.fits')
		else:
			#for normal running
			master_table_hits_spec = fits.getdata( '../saved_files/tabla_maestra_hits_spec.fits')
			master_table_pares_spec = fits.getdata( '../saved_files/tabla_maestra_pares_spec.fits')

		master_table_hits_spec = QTable(master_table_hits_spec)
		master_table_pares_spec = QTable(master_table_pares_spec)
		
		cond = master_table_hits_spec['id_qso_hs']!=-1
		master_table_hits_spec = master_table_hits_spec[cond]

		if z_cl_type =='spec':	
			return(master_table_hits_spec, master_table_pares_spec)
	else:
		raise ValueError('Select a valid option, "spec" or "phot". ')

def dndz_vs_x(x, dir_name, pares, hits, tipo, begin, end, n, z_min, z_max):
	'''Calls to the fucntions to calculate dn/dz = dn/dz(x)
 	
	Inputs 
		x:	Variable against wich calculate dn/dz (str). 
			Accepted values: 	'com', 'pro' or 'r200' for impact parameter. 
						'ew' for equvalent with
		dir_name:	Directory name fore save the results (str).
		master_table_pares:	Master tababe with all the lines of sigth pairs (QTable).
		master_table_hits:	Master tababe with all the possible hits (QTable).
		tipo: 	Type of grid to be used for x (str). 
			Accepted values:	log: Logarithmic grid.
						lin: Linear grid.
			    			snh: Same number of hits per bin.
			    			snp: Same number of pairs per bin.
			    			neil: Grid used in neil's code .
			    			one: One single bin.
		begin:	Starting point of the grid (float). 
		end:	Ending point of the grid (float). 
		n:	Number of point in the grid (int).
		z_min:		Minimum redshift redshift where detection is possible (float). Given by survey limits
		z_max:		Maximum redshift redshift where detection is possible (float). Given by survey limits

	Outputs:
		grid_pares:
		grid_hits:
		red_path:
		gc:
		results:
	'''
	#create bins
	n_bins = np.log2(len(hits)) #+ 1.
	if  tipo == 'snh':
		edges = bins_edges(hits, tipo, x, begin, end, n_bins, n)
	else:
		edges = bins_edges(pares, tipo, x, begin, end, n_bins, n)

	if x == 'ew':
		grid_pares = grilla(pares, edges, 'ew_pares')
	else:
		grid_pares = grilla(pares, edges, x)
	with open(dir_name + '/grid_pares'+'.pickle', 'wb') as f:
		pickle.dump(grid_pares, f, protocol=2)

	grid_hits = grilla(hits, edges, x)
	with open(dir_name + '/grid_hits'+'.pickle', 'wb') as f:
		pickle.dump(grid_hits, f, protocol=2)

	red_path, gc = redshift_path(grid_pares, z_min, z_max)			
	with open(dir_name + '/dz'+'.pickle', 'wb') as f:
		pickle.dump(red_path, f, protocol=2)
	with open(dir_name+ '/gc'+'.pickle', 'wb') as f:
		pickle.dump(gc, f, protocol=2)

	# import pdb; pdb.set_trace()
	results = results_table(x, grid_hits, red_path, grid_pares)
	with open(dir_name + '/results'+'.pickle', 'wb') as f:
		pickle.dump(results, f, protocol=2)

	return(grid_pares, grid_hits, red_path, gc, results)

# ------- Main -------
if __name__ == '__main__':

	limit_by = names(ls)
	dist_type = names(dd)

	grid_str =  '{}_n{:.1f}'.format(grt, n)
	mass_str = '-mass_10e'+str(np.log10(min_mass.value))+'_to_10e'+str(np.log10(max_mass.value))
	ew_str = '-rew_'+str(min_EW.value)+'_to_'+str(max_EW.value)
	z_str = '-z_{:.2f}_to_{:.2f}'.format(min_z, max_z)
	ip_str = '-ip_{}_{:1f}_to_{:1f}'.format(dd, min_IP, max_IP)

	if sn == 'global':
		sn_str = '-s_{:.1f}_global'.format(s)
	elif sn == 'local':
		sn_str = '-s_{:.1f}_local'.format(s)

	# read the catalogs
	cluster, mgii, mgii_search = load_catalogs(ls)
	# transform the catalogs from Fits to astropy QTables. Some cuts are made to the catalogs.
	cumulo, cumulo_spec, cumulo_phot,cumulo_rejected, cumulo_0 = cluster_table(cluster, z_min, z_max, min_mass, max_mass)
	mgii, qso, mgii_rejected, mgii_0 = mgii_table(mgii, mgii_search, min_EW, max_EW, 2.2827778, ls)

	#load master tables and make cuts corresponding to input values
	if m !=0 :
		
		# run for spec
		if z_cl_type == 'spec' or  z_cl_type =='all':
			#load spectroscopic master tables
			master_table_hits_spec, master_table_pares_spec = load_master_tables('spec', running_mode)
			#cuts corresponding to the input or default values  on the master tables
			pares_spec, pares_spec_rejected  = cuts_pares(s, z_min, z_max, min_EW, min_mass, max_mass, min_z, max_z, min_IP, max_IP, dist_type, master_table_pares_spec)
			hits_spec, hits_spec_rejected = cuts_hits(s, z_min, z_max, min_EW, max_EW, min_mass, max_mass, min_z, max_z, min_IP, max_IP, dist_type, master_table_hits_spec)

		#run for phot
		if z_cl_type =='phot' or z_cl_type == 'all':
			#load photmetric master tables 
			master_table_hits_phot, master_table_pares_phot = load_master_tables('phot', running_mode)
			#cuts corresponding to the input or default values  on the master tables
			pares_phot, pares_phot_rejected = cuts_pares(s, z_min, z_max, min_EW, min_mass, max_mass, min_z, max_z, min_IP, max_IP, dist_type, master_table_pares_phot)
			hits_phot, hits_phot_rejected = cuts_hits(s, z_min, z_max, min_EW, max_EW, min_mass, max_mass, min_z, max_z, min_IP, max_IP, dist_type, master_table_hits_phot)
	
	# run for dn/dz vs b, dn/dz vs ew and dn/dz vs z
	if m == 1 or m ==2 or m==3:

		# define directory name and assign  some variable values
		if m ==1:#dn/dz vs b
			dir_name = '../saved_files/dndz_v_b/'+grid_str+'-'+limit_by+mass_str+ew_str+z_str+ip_str+sn_str
			create_dir(dir_name)
			x_value = dd
			first = min_IP
			last = max_IP

		elif m == 2:#dn/dz vs ew
			dir_name = '../saved_files/dndz_v_ew/'+grid_str+'-'+limit_by+mass_str+ew_str+z_str+ip_str+sn_str
			create_dir(dir_name)
			x_value = 'ew'
			first = min_EW.value
			last = max_EW.value 

		elif m == 3:#dn/dz vs z
			dir_name = '../saved_files/dndz_v_z/'+grid_str+'-'+limit_by+mass_str+ew_str+z_str+ip_str+sn_str
			create_dir(dir_name)
			x_value = 'z' 
			first = min_z
			last = max_z

		#run for spectroscopic clusters
		if z_cl_type == 'spec' or z_cl_type == 'all':
			print('\nRuning for spectroscopic clusters\n')		
			#define sub directory
			dir_name_spec = dir_name+'/spec'
			if running_mode == 'spec_p':
				dir_name_spec = dir_name+'/spec_p'
			create_dir(dir_name_spec)
			#call to the function that makes the calculations
			grid_pares_spec, grid_hits_spec, red_path_spec, gc_spec, results_spec = dndz_vs_x(x_value, dir_name_spec, pares_spec, hits_spec, grt, first, last, n, z_min, z_max)

		#run for photometric clusters
		if z_cl_type == 'phot' or z_cl_type=='all':
			print('\nRuning for photometric clusters\n')
			#define sub directory
			dir_name_phot = dir_name+'/phot'	
			create_dir(dir_name_phot)
			#call to the functions that makes the calculations
			grid_pares_phot, grid_hits_phot, red_path_phot, gc_phot, results_phot = dndz_vs_x(x_value, dir_name_phot, pares_phot, hits_phot, grt, first, last, n, z_min, z_max)

		#Shows plot with the results
		if plot_result == 'yes':
			#create figure
			fig = plt.figure()
			legend = []

			#plot for sum of phot and spec results, NOT IN USE
			# if z_cl_type == 'all':
				# plot_sum(results_spec,results_phot,'r',-0.95)
			
			#plot for spectroscopic results
			if z_cl_type == 'spec' or z_cl_type == 'all':
				print('\nResults spec:\n')
				print(results_spec)
				plot_results_show(fig, x_value, results_spec, color='red')
				legend_spec = 'mass 10e{} Msun to 10e{} Msun, ew {} to {}, z {} to {}, s {}, {}, {}'.format(np.log10(min_mass.value), np.log10(max_mass.value), min_EW, max_EW, min_z, max_z, s, 'spec', limit_by)
				legend.append(legend_spec)
			
			#plot for photpmetric results
			if z_cl_type == 'phot' or z_cl_type=='all':
				print('\nResults phot:\n')
				print(results_phot)
				plot_results_show(fig, x_value, results_phot, color='blue')
				legend_phot = 'mass 10e{} Msun to 10e{} Msun, ew {} to {}, z {} to {}, s {}, {}, {}'.format(np.log10(min_mass.value), np.log10(max_mass.value), min_EW, max_EW, min_z, max_z, s, 'spec', limit_by)
				llegend.append(legend_phot)
			
			#add legend and display plot
			plt.legend(legend, loc='best', fontsize='medium' , numpoints=1)
			plt.show()
			fig.clf

	#Plot redshift distribution of MgII absorbers and galaxy clusters from the catalogs. 
	elif m == 0:
		# mass bins for galaxy clusters.
		cumulo_1 = cumulo_0[10**13.6 < cumulo_0['mass'].value]
		cumulo_1 = cumulo_1[cumulo_1['mass'].value < 10**14.0]
		cumulo_2 = cumulo_0[10**14.0 < cumulo_0['mass'].value]
		cumulo_2 = cumulo_2[cumulo_2['mass'].value < 10**14.2]
		cumulo_3 = cumulo_0[10**14.2 < cumulo_0['mass'].value]
		cl_mass = [cumulo_0, cumulo_1, cumulo_2, cumulo_3]
		cl_mass_label = ['Clusters', 'Clusters with 10^13.6'+'['+r'$M_{\odot}$'+'] < M < 10^14.0'+'['+r'$M_{\odot}$'+']','Clusters with 10^14.0'+'['+r'$M_{\odot}$'+'] < M < 10^14.2'+'['+r'$M_{\odot}$'+']', 'Clusters with 10^14.2'+'['+r'$M_{\odot}$'+'] < M']
		
		# Equivalent width bins for MgII absorbers.
		mgii_1 = mgii_0[0.6 < mgii_0['rew_mgii_2796']]
		mgii_1 = mgii_1[mgii_1['rew_mgii_2796'] < 1.0]
		mgii_2 = mgii_0[1.0 < mgii_0['rew_mgii_2796']]
		mgii_2 = mgii_2[mgii_2['rew_mgii_2796'] < 1.5]
		mgii_3 = mgii_0[1.5 < mgii_0['rew_mgii_2796']]
		mgii_3 = mgii_3[mgii_3['rew_mgii_2796'] < 2.0]
		mgii_rew = [mgii_0, mgii_1, mgii_2, mgii_3]
		mgii_rew_label = ['MgII 2796', 'MgII 2796 with 0.6 [AA] < REW < 1.0 [AA]', 'MgII 2796 with 1.0 < REW < 1.5 [AA]', 'MgII 2796 with 1.5 [AA] < REW < 2.0 [AA]']
		
		# Redshift type bins for galaxy clusters, photometric or spectroscopic
		cluster_list = [cumulo_0, cumulo_0[cumulo_0['z_spec'] != -1], cumulo_0[cumulo_0['z_spec'] == -1]]
		cluster_list_label = ['Clusters', 'Spectroscopic Clusters', 'Photometric Clusters']
		
		#call to function that makes the plot
		# plot_overlap(cl_mass, cl_mass_label, mgii_rew, mgii_rew_label)#Plot for mass bins and ew bins	
		plot_overlap(cluster_list, cluster_list_label, [mgii], ['MgII'])#Plot for redshift type bins and ew bins
		print('{} galaxy clusters, {} MgII absorbers and {} Qso'.format(len(cumulo_0), len(mgii_0), len(qso)))
	
	#Run fo dn/(dz/dw) vs ew (In development)
	elif  m==4:
		#Create directory
		dir_name = '../saved_files/dndzdw_v_ew/'+grid+'-'+limit_by+'-mass_10e'+str(np.log10(min_mass.value))+'_to_10e'+str(np.log10(max_mass.value))+'-rew_'+str(min_EW.value)+'_to_'+str(max_EW.value)+'-z_{:.2f}_to_{:.2f}'.format(min_z, max_z)+'-'+sn_str
		create_dir(dir_name)
		#create ew grid
		# ew_grid =

		#run for spectroscopic clusters
		if z_cl_type == 'spec' or z_cl_type == 'all':
			print('\nRuning for spectroscopic clusters\n')		
			#define sub directory
			dir_name_spec = dir_name+'/spec'
			if running_mode == 'spec_p':
				dir_name_spec = dir_name+'/spec_p'
			create_dir(dir_name_spec)

			# results_spec = QTable[]
			#iteration over ew grid
			for i in range(1, len(ew_grid)):
				left_lim = ew_grid(i - 1)
				right_lim = ew_grid(i)
				#call to dn/dz function for width
				# grid_pares_spec, grid_hits_spec, red_path_spec, gc_spec, results = dndz_vs_x(dd, dir_name_spec, pares_spec, hits_spec, grt, n, z_min, z_max)
				#add results to dn/(dzdw) vs ew table
				row = results[0]
				results_spec.add_row(row)
				pass	

		#run for photometric clusters
		if z_cl_type == 'phot' or z_cl_type=='all':
			print('\nRuning for photometric clusters\n')
			#define sub directory
			dir_name_phot = dir_name+'/phot'	
			create_dir(dir_name_phot)

			# results_phot = QTable[]
			for i in range(1, len(ew_grid)):
				left_lim = ew_grid(i - 1)
				right_lim = ew_grid(i)
				#call to dn/dz function for width
				# grid_pares_phot, grid_hits_phot, red_path_phot, gc_phot, results = dndz_vs_x(dd, dir_name_phot, pares_phot, hits_phot, grt, n, z_min, z_max)
				#add results to dn/(dzdw) vs ew table
				pass	
		#show results

	elif m == 5:
		# hay que modificar esta parte para que funcione con las modificaciones que se hcieron a las funciones
		if z_cl_type == 'spec' or z_cl_type == 'all':

			print('\nRuning for spectroscopic clusters\n')		
			
			dir_name_spec = dir_name+'/spec'	
			create_dir(dir_name_spec)

			n_bins = np.log2(len(hits_spec)) + 1.
			edges = bins_edges(n_bins, n)

			try:
				with open(dir_name_spec + '/grid_pares_'+str(dd)+'.pickle', 'rb') as f:
					grid_pares_spec = pickle.load(f)
			except UnicodeDecodeError:
				with open(dir_name_spec + '/grid_pares_'+str(dd)+'.pickle', 'rb') as f:
					grid_pares_spec = pickle.load(f, encoding='latin1')		
			except OSError:
				grid_pares_spec = grilla(pares_spec, edges, dd)
				with open(dir_name_spec + '/grid_pares_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(grid_pares_spec, f, protocol=2)
			except IOError:
				grid_pares_spec = grilla(pares_spec, edges, dd)
				with open(dir_name_spec + '/grid_pares_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(grid_pares_spec, f, protocol=2)

			try:
				with open(dir_name_spec + '/grid_hits_'+str(dd)+'.pickle', 'rb') as f:
					grid_hits_spec = pickle.load(f)
			except UnicodeDecodeError:
				with open(dir_name_spec + '/grid_hits_'+str(dd)+'.pickle', 'rb') as f:
					grid_hits_spec = pickle.load(f, encoding='latin1')		
			except OSError:
				grid_hits_spec = grilla(hits_spec, edges, dd)
				with open(dir_name_spec + '/grid_hits_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(grid_hits_spec, f, protocol=2)
			except IOError:
				grid_hits_spec = grilla(hits_spec, edges, dd)
				with open(dir_name_spec + '/grid_hits_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(grid_hits_spec, f, protocol=2)

			try:
				with open(dir_name_spec + '/dz_'+str(dd)+'.pickle', 'rb') as f:
					red_path_spec = pickle.load(f)
				with open(dir_name_spec + '/gc_'+str(dd)+'.pickle', 'rb') as f:
					gc_spec = pickle.load(f)
			except UnicodeDecodeError:
				with open(dir_name_spec + '/dz_'+str(dd)+'.pickle', 'rb') as f:
					red_path_spec = pickle.load(f, encoding='latin1')
				with open(dir_name_spec + '/gc_'+str(dd)+'.pickle', 'rb') as f:
					gc_spec = pickle.load(f, encoding='latin1')	
			except OSError:
				red_path_spec, gc_spec = redshift_path(grid_pares_spec, z_min, z_max)			
				with open(dir_name_spec + '/dz_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(red_path_spec, f, protocol=2)
				with open(dir_name_spec + '/gc_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(gc_spec, f, protocol=2)
				plot_gc(gc_spec, 'spec')
			except IOError:
				red_path_spec, gc_spec = redshift_path(grid_pares_spec, z_min, z_max)			
				with open(dir_name_spec + '/dz_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(red_path_spec, f, protocol=2)
				with open(dir_name_spec + '/gc_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(gc_spec, f, protocol=2)				
				plot_gc(gc_spec, 'spec')
			try:
				with open(dir_name_spec + '/results_'+str(dd)+'.pickle', 'rb') as f:
					results_spec = pickle.load(f)
			except UnicodeDecodeError:
				with open(dir_name_spec + '/results_'+str(dd)+'.pickle', 'rb') as f:
					results_spec = pickle.load(f, encoding='latin1')
			except OSError:
				results_spec = results_table(grid_hits_spec, red_path_spec, grid_pares_spec, dd)
				with open(dir_name_spec + '/results_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(results_spec, f, protocol=2)
				# plot_results(results_spec, dir_name_spec, dd)	
			except IOError:	
				results_spec = results_table(grid_hits_spec, red_path_spec, grid_pares_spec, dd)
				with open(dir_name_spec + '/results_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(results_spec, f, protocol=2)
				# plot_results(results_spec, dir_name_spec, dd)

		if z_cl_type == 'phot' or z_cl_type=='all':

			print('\nRuning for photometric clusters\n')

			dir_name_phot = dir_name+'/phot'
			create_dir(dir_name_phot)

			n_bins = np.log2(len(hits_phot)) + 1.
			edges = bins_edges(n_bins, n)

			try:
				with open(dir_name_phot + '/grid_pares_'+str(dd)+'.pickle', 'rb') as f:
					grid_pares_phot = pickle.load(f)
			except UnicodeDecodeError:
				with open(dir_name_phot + '/grid_pares_'+str(dd)+'.pickle', 'rb') as f:
					grid_pares_phot = pickle.load(f, encoding='latin1')		
			except OSError:
				grid_pares_phot = grilla(pares_phot, edges, dd)
				with open(dir_name_phot + '/grid_pares_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(grid_pares_phot, f, protocol=2)
			except IOError:
				grid_pares_phot = grilla(pares_phot, edges, dd)
				with open(dir_name_phot + '/grid_pares_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(grid_pares_phot, f, protocol=2)

			try:
				with open(dir_name_phot + '/grid_hits_'+str(dd)+'.pickle', 'rb') as f:
					grid_hits_phot = pickle.load(f)
			except UnicodeDecodeError:
				with open(dir_name_phot + '/grid_hits_'+str(dd)+'.pickle', 'rb') as f:
					grid_hits_phot = pickle.load(f, encoding='latin1')		
			except OSError:
				grid_hits_phot = grilla(hits_phot, edges, dd)
				with open(dir_name_phot + '/grid_hits_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(grid_hits_phot, f, protocol=2)
			except IOError:
				grid_hits_phot = grilla(hits_phot, edges, dd)
				with open(dir_name_phot + '/grid_hits_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(grid_hits_phot, f, protocol=2)

			try:
				with open(dir_name_phot + '/dz_'+str(dd)+'.pickle', 'rb') as f:
					red_path_phot = pickle.load(f)
				with open(dir_name_phot + '/gc_'+str(dd)+'.pickle', 'rb') as f:
					gc_phot = pickle.load(f)
			except UnicodeDecodeError:
				with open(dir_name_phot + '/dz_'+str(dd)+'.pickle', 'rb') as f:
					red_path_phot = pickle.load(f, encoding='latin1')
				with open(dir_name_phot + '/gc_'+str(dd)+'.pickle', 'rb') as f:
					gc_phot = pickle.load(f, encoding='latin1')	
			except OSError:
				red_path_phot, gc_phot = redshift_path(grid_pares_phot, z_min, z_max)			
				with open(dir_name_phot + '/dz_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(red_path_phot, f, protocol=2)
				with open(dir_name_phot + '/gc_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(gc_phot, f, protocol=2)
				plot_gc(gc_phot, 'phot')
			except IOError:
				red_path_phot, gc_phot = redshift_path(grid_pares_phot, z_min, z_max)			
				with open(dir_name_phot + '/dz_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(red_path_phot, f, protocol=2)
				with open(dir_name_phot + '/gc_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(gc_phot, f, protocol=2)	
				plot_gc(gc_phot, 'phot')
			try:
				with open(dir_name_phot + '/results_'+str(dd)+'.pickle', 'rb') as f:
					results_phot = pickle.load(f)
			except UnicodeDecodeError:
				with open(dir_name_phot + '/results_'+str(dd)+'.pickle', 'rb') as f:
					results_phot = pickle.load(f, encoding='latin1')
			except OSError:
				results_phot = results_table(grid_hits_phot, red_path_phot, grid_pares_phot, dd)
				with open(dir_name_phot + '/results_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(results_phot, f, protocol=2)
				# plot_results(results_phot, dir_name_phot, dd)	
			except IOError:	
				results_phot = results_table(grid_hits_phot, red_path_phot, grid_pares_phot, dd)
				with open(dir_name_phot + '/results_'+str(dd)+'.pickle', 'wb') as f:
					pickle.dump(results_phot, f, protocol=2)
				#plot_results(results_phot, dir_name_phot, dd)

	else:
		raise ValueError('''Invalid value for -m, choose '1', '2' or '3' 

					0: Plot Redshift Distribution of MgII absorbers and galaxy clusters
					1: Run dn/dz vs Impact Parameter (default) 
					2: Run dn/dz vs Equivalent Width for r< R200 & r>R200.
					3: Run dn/dz vs Redshift 
					4: Run dn/(dzdw) vs Equivalent Width (in development).
					5, 6 : Run from saved files (in  development) 
					''')