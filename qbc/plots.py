import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from xkcd_rgb import *
from astropy.io import fits
from astropy.table import Table,QTable, Column

#
plt.rcParams['figure.figsize'] = [4.0, 3.0]
plt.rcParams['xtick.labelsize'] = 17
plt.rcParams['ytick.labelsize'] = 17
plt.rcParams['xtick.major.size'] = 12.0
plt.rcParams['ytick.major.size'] = 12.0
plt.rcParams['xtick.minor.size'] = 6.0
plt.rcParams['ytick.minor.size'] = 6.0
plt.rcParams['axes.labelsize'] = 22

def plot_dNdz_vs_x(X, ax, tab, xlabel=None, ylabel=None, annotate=False, title=None, annotate_size=12, **kwargs):
	"""
	Parameters

		X : independent variable. impact parameter (b) or equivalent width (ew)
		ax : axis
			Matplotlib axis
		tab : QTable
			Table with the results, must contain certain columns
		outname : str, optional
			Name of the output file to save the figure if given
		xlabel : str, optional
			Label of x-axis if given 
		ylabel : str, optional
			Label of y-axis if given
		annotate : bool, optional
			Whether to annotate info (hits and dz)
		title: str, optional
			Title of the graph

	Return
	------
		Plots into ax.

	"""

	if X == 'b':
		x = (tab['b_i']+tab['b_f'])/2.0
		xerr = (tab['b_f']-tab['b_i'])/2.0
	elif X == 'ew':
		x = (tab['ew_i']+tab['ew_f'])/2.0
		xerr = (tab['ew_f']-tab['ew_i'])/2.0
	elif X == 'z':
		x = (tab['z_i']+tab['z_f'])/2.0
		xerr = (tab['z_f']-tab['z_i'])/2.0
	cond = tab['dn/dz'] > 0
	aux = tab[cond]
	y = np.zeros(len(tab))+0.01
	for i in range(len(y)):
		if i%2==0:
			y[i]+=19	
	# y = 31
	yerr = [[i['sigma-_dn/dz'], i['sigma+_dn/dz']] for i in tab]
	yerr = np.column_stack(yerr)  
	ax.errorbar(x, tab['dn/dz'], xerr=xerr, yerr=yerr, **kwargs)
	if xlabel is not None:
		ax.set_xlabel(xlabel)
	if ylabel is not None:
		ax.set_ylabel(ylabel)
	if title is not None:
		ax.set_title(title)
	if annotate:
		for ii, row in enumerate(tab):
			# import pdb; pdb.set_trace()
			s1 = '{:.0f}'.format(row['dn'])
			s2 = '{:.2f}'.format(row['dz'])
			s3 = '{:.0f}'.format(row['#pares'])
			if len(s1) > 5:
				s1 ='{:.1e}'.format(row['dn'])
			if len(s2) > 5:
				s2 ='{:.1e}'.format(row['dz'])
			if len(s3) > 5:
				s3 ='{:.1e}'.format(row['#pares'])		
			s = s1+'\n'+s2+'\n'+s3
			# s = '{:.0f}\n{:.2f}\n{:.0f}'.format(row['dn'], row['dz'], row['#pares'])
			x_annotate = x[ii]
			y_annotate = y[ii]
			if 'c' in  kwargs.keys():
				c_annotate = kwargs['c']
			else:
				c_annotate = 'grey'
			# ax.annotate(s, (x_annotate, y_annotate), rotation=0, ha='right', color=c_annotate, size=annotate_size, stretch='ultra-condensed')

def plot_neil():	
	from QbC_mgii_v6 import err_plus, err_minus, w_min
	files = {	
	'neil_spec':'../saved_files/neil-lim_civ-mass_10e13.6_to_10e20.0-rew_0.6_to_10.0-s_0.0/spec/results_pro.pickle',
	'neil_phot':'../saved_files/neil-lim_civ-mass_10e13.6_to_10e20.0-rew_0.6_to_10.0-s_0.0/phot/results_pro.pickle',
	}

	tables = {}
	for i in files:
		try:
			with open(files[i], 'rb') as f:
				tables[i]= pickle.load(f)        
		except UnicodeDecodeError:
			with open(files[i], 'rb') as f:
				tables[i] = pickle.load(f, encoding='latin1') 
	# import pdb;pdb.set_trace()
	aux = tables['neil_spec']['bi', 'bf', 'dn', 'dz', '#pares']
	aux['dn'] = aux['dn'] + tables['neil_phot']['dn']
	aux['dz'] = aux['dz'] + tables['neil_phot']['dz']	
	aux['#pares'] = aux['#pares'] + tables['neil_phot']['#pares']
	dndz = Column(np.zeros(len(aux))-1, name='dn/dz')
	sigma_n_mas = Column(np.zeros(len(aux))-1, name='sigma+_n')
	sigma_n_menos = Column(np.zeros(len(aux))-1, name='sigma-_n')
	dndz_err_mas = Column(np.zeros(len(aux))-1, name='sigma+_dn/dz')
	dndz_err_menos = Column(np.zeros(len(aux))-1, name='sigma-_dn/dz')
	aux.add_column(dndz)
	aux.add_column(dndz_err_mas)
	aux.add_column(dndz_err_menos)
	aux.add_column(sigma_n_mas, index=3)
	aux.add_column(sigma_n_menos, index=4)
	for i in range(len(aux)):	
		dn = aux['dn'][i]
		dz = aux['dz'][i]
		aux['sigma+_n'][i] = err_plus(dn)
		aux['sigma-_n'][i] = err_minus(dn)
		if dz != 0:
			aux['dn/dz'][i] = dn/dz
			aux['sigma+_dn/dz'][i] = 1.0/dz*err_plus(dn)
			aux['sigma-_dn/dz'][i] =1.0/dz*err_minus(dn)	
	
	fig = plt.figure(figsize=(10,10))
	# fig.suptitle(title, fontsize=12)
	ax = fig.add_subplot(111)
	s='[Mpc](proper)'
	plot_dNdz_vs_x('b', ax, aux , xlabel='Impact parameter'+s, ylabel='dn/dz', marker='o', annotate=True, c='r', capsize=0, lw=0.5, ms=18, ls='None')
	plt.hlines(0.3,0,10,linestyles = 'dashed' )
	ax.set_ymargin(0.2)
	ax.set_xmargin(0.2)
	ax.set_xlim(0,10)
	ax.set_ylim(0,10)
	ax.set_xscale('log')
	ax.set_yscale('log')
	# fig.savefig(table+'.png')
	plt.show()

def plots_wmin():
	# master_table_hits_spec = fits.getdata( '../../../../QbC/tabla_maestra_hits_spec_new.fits')
	master_table_pares_spec = fits.getdata( '../../../../QbC/tabla_maestra_pares_spec_new.fits')
	# master_table_hits_phot = fits.getdata( '../../../../QbC/tabla_maestra_hits_phot_new.fits')
	# master_table_pares_phot = fits.getdata( '../../../../QbC/tabla_maestra_pares_phot_new.fits')
	# master_table_hits_spec = QTable(master_table_hits_spec)
	master_table_pares_spec = QTable(master_table_pares_spec)	
	# master_table_hits_phot = QTable(master_table_hits_phot)
	# master_table_pares_phot = QTable(master_table_pares_phot)

	master_table_pares_spec['w_min_s1_local'] = w_min(master_table_pares_spec['s2n_local'], master_table_pares_spec['z_cluster'], 1)
	fig = plt.figure(figsize=(10,10))
	# fig.suptitle(title, fontsize=12)
	ax = fig.add_subplot(121)
	ax.set_xlabel('w_min_old')
	ax.set_ylabel('w_min_new')
	# plt.scatter(master_table_pares_spec['w_min_s1'],master_table_pares_spec['w_min_s1_local'])
	plt.hist2d(master_table_pares_spec['w_min_s1'],master_table_pares_spec['w_min_s1_local'], 500, vmax=6000)
	plt.colorbar()
	ax = fig.add_subplot(122)
	ax.set_xlabel('w_min_old')
	ax.set_ylabel('w_min_new - w_min_old')
	# plt.scatter(master_table_pares_spec['w_min_s1'],master_table_pares_spec['w_min_s1_local'] - master_table_pares_spec['w_min_s1'])
	plt.hist2d(master_table_pares_spec['w_min_s1'],master_table_pares_spec['w_min_s1_local'] - master_table_pares_spec['w_min_s1'], 500, vmax=6000)	
	plt.colorbar()
	plt.show()

def plot_results(table, name, dist, show=False):

		index_0 = name.index('mass')
		cond = table['dn/dz'] > 0
		aux = table[cond]
		if '/spec' in name:
			z = 'spectroscopic'
			o2 = '/spec'
			index_1 = name.index('/spec')
			color = 'r'
		elif '/phot' in name:
			z = 'photometric'
			o2 = '/phot'
			index_1 = name.index('/phot')
			color='blue'
		title = name[index_0:index_1]
		fig = plt.figure(figsize=(15,7.5))
		fig.suptitle(title, fontsize=20)
		ax = fig.add_subplot(111)
		if dist == 'pro' :
			s='[Mpc](proper)'
			o1 = '/pro'
		elif dist == 'com' :
			s='[Mpc](comoving)'
			o1 = '/com'
		elif dist == 'r200' :
			s='(R200)'
			o1 = '/r200'
		dir_name = '../figs/'+name[15:index_1]+o1+o2
		outname = dir_name+'/dndz.png'
		# import pdb;pdb.set_trace()
		if os.path.isdir(dir_name):
			pass
		else:
			os.makedirs(dir_name)
		plot_dNdz(ax, table , marker='o', annotate=True, c=color, capsize=0, lw=0.5, ms=18, ls='None')
		ax.set_ylabel('dn/dz', fontsize=18)
		ax.set_xlabel('Impact parameter'+s, fontsize=18)
		if 'log'  in outname:	
			ax.set_xlim(0.1, 40)
			ax.set_ylim(0.01, 30)
			ax.set_xscale('log')
			ax.set_yscale('log')
		else:
			ax.set_xlim(0, 40)
			ax.set_ylim(0.01, 30)
			ax.set_yscale('log')
		ax.legend([z], loc='best', fontsize='medium', scatterpoints=1)
		if show:
			pass
			#plt.show()
		else:
			fig.savefig(outname)
			print('Figure saved '+outname)
		fig.clf

def plot_results_show(fig, x, table, title='', color='cobalt'):

		# fig = plt.figure()
		color = get_xkcd_color(color)
		fig.suptitle(title, fontsize=20)
		ax = fig.add_subplot(111)
		ax.set_ylabel('dn/dz', fontsize=18)
		if (x == 'com' or x == 'pro') or x =='r200':
			ax.set_xlabel('Impact parameter', fontsize=18)
			x='b'
		elif x =='ew':
			ax.set_xlabel('Equivalent Width [AA]', fontsize=18)
		ax.set_xlim(0.1, 40)
		ax.set_ylim(0.01, 30)
		ax.set_xscale('log')
		ax.set_yscale('log')
		plot_dNdz_vs_x(x, ax, table , marker='o', annotate=True, c=color, capsize=0, lw=0.5, ms=18, ls='None')
		# plt.show()
		# fig.clf

def plot_sum(tabla1,tabla2,col,pos_t):
	'''Plots the number of hits divided by the redshift path as a function of the impact parameter for the sum of the results contained 
	in 'tabla1' and 'tabla2'.
	inputs:
		tabla1: Table with the results (QTable).
		tabla2: Table with the results (QTable).
		col: color (str).
		pos_t: posicion texto (float).
	'''
	angulo=-35
	if grt == 'log':
		angulo = -20
	for i in range(len(tabla1)):	
		aux = '{:.0f}/{:.2f}/{:.0f}'.format(tabla1['dn'][i]+tabla2['dn'][i],tabla1['dz'][i]+tabla2['dz'][i],tabla1['#pares'][i]+tabla2['#pares'][i])
		plt.annotate(aux,(tabla1['b'][i],pos_t), rotation= angulo,color=col)
	plt.scatter(tabla1['b'],(tabla1['dn']+tabla2['dn'])/(tabla1['dz']+tabla2['dz']), c=col, label='line3')

def plot_gc(gc, zct):
	'''Plots the redshift path density for different impact parameters in Mpc.

	inputs:
		gc:
		zct:
	'''
	dir_name = '../figs/'+grid+'-'+limit_by+'-mass_10e'+str(np.log10(min_mass.value))+'_to_10e'+str(np.log10(max_mass.value))+'-rew_'+str(min_EW.value)+'_to_'+str(max_EW.value)+'-'+'s_'+str(s)+'/'+str(dd)+ '/'+zct
	if os.path.isdir(dir_name):
		pass
	else:
		os.makedirs(dir_name)

	aux = gc.keys()
	aux = list(aux)
	aux.sort()
	
	for b in aux:
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)	
		ax.set_xlabel('Z')
		ax.set_ylabel(r'$g_{c}$')
		ax.plot(gc[b]['z'],gc[b]['gc'], ls='steps')
		ax.set_title('Redshift path density (b='+str(b)+' Mpc)')
		figname = dir_name +'/gc_'+'{:.3f}'.format(b)+'Mpc.png'
		fig.savefig(figname)
		#fig.show()
		fig.clf
		
	return()

def plot_overlap(clusters, clusters_label, mgii, mgii_label, z_bins=np.arange(0,2.6,0.005)):
	'''Plot an overlap of the histogram of the number of galaxy clusters per redshift bin and the histogram of the number of Qso per
	redshift bin.

	inputs:
		clusters: Llist of tables with the galaxy clusters data (QTable).
		clusters_label: List with labels for he clusters
		mgii: List of tables with the Qso data (QTable).
		mgii_label: List with mgii labels
		z_bins(optional): Array with the bins edges (numpy.ndarray).
	'''
	lw=[2,1,1,1]
	colors_cl = cm.winter(np.linspace(0, 1, len(mgii)+2))
	colors_mgii = cm.autumn(np.linspace(0, 1, len(mgii)+2))
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	for i in range(len(clusters)):
		ax.hist(clusters[i]['z'],z_bins,histtype='step', color=colors_cl[i], label=clusters_label[i],linewidth=lw[i])
	for i in range(len(mgii)):
		ax.hist(mgii[i]['zabs'],z_bins,histtype='step', color=colors_mgii[i], label=mgii_label[i],linewidth=lw[i])
	ax.set_ylabel('Frequency')
	ax.set_xlabel('Redshift')
	ax.legend()
	ax.set_xlim(0,1.25)
	ax.set_ylim(0.850)
	#plt.title(r'Frequecy vs Redshift')
	#fig.savefig('/home/hector/Dropbox/Tejos_Salas/figs/overlap.png')
	plt.show()


def plot_dist_cluster(cluster, spec, phot):
	'''Plot a histogram of the number of galaxy clusters per redshift bin separeted in photometric and spectroscopic.

	inputs:
		cluster: Table with the galaxy clusters data (QTable).
		spec: Table with the spectroscopic galaxy clusters data (QTable).
		phot: Table with the photometric galaxy clusters data (QTable).
		'''
	clusters = [spec, phot, cluster]
	cl_label = ['spec','phot', 'all']
	z_bins=np.arange(0,2.6,0.005)
	colors = cm.autumn(np.linspace(0, 1, len(clusters)+2))
	lw=[2,2,2.5]
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	for i in range(len(clusters)):
		ax.hist(clusters[i]['z'],z_bins,histtype='step', color=colors[i], label=cl_label[i],linewidth=lw[i])
	ax.set_ylabel('Frequency')
	ax.set_xlabel('Redshift')
	ax.legend()
	ax.set_xlim(0,1.25)
	ax.set_ylim(0.850)
	plt.title(r'Cluster distribution')
	#fig.savefig('/home/hector/Dropbox/Tejos_Salas/figs/cluster_distribution.png')
	plt.show()

def plot_ntr06(ax, N, W):
	"""
	Parameters

		ax : axis
			Matplotlib axis
	Return
	------
		Plots into ax.

	"""
	x = np.arange(0.01, 10, 0.01)
	y = N/W*np.exp(-x/W)
	ax.plot(x, y, 'k--', lw=2, label='_nolegend_')

if __name__ == '__main__':
	pass
	# plots_wmin()
	# plot_neil()