import os

"""This script will run QbC_mgii_v6.py for different cuts in 
different parameters. This will create results_tables and plots 
that will be stored in corresponding directories"""

# define the cuts
# cl_types = ['spec', 'phot']
cl_types = ['spec']
# sn_types = ['local', 'global']
sn_types = ['local']
#grid_types = ['log', 'lin', 'snp']
grid_types = [ 'log']
# dist_types = ['pro', 'com', 'r200']
dist_types = ['com']
#iplims =[(0.1, 12), (0.1, 40)]
iplims =[(0.1, 12)]
ewlims = [(0.6, 1.0), (1.0, 1.5), (1.5, 2.0), (0.6, 2.0)]
# ewlims = [(1.5, 2.0),(0.6, 2.0)]
masslims = [(13.6, 14.0), (14.0, 14.2), (14.20, 16.0), (13.6, 16.0)]
# zlims = [(0.36,0.44), (0.44,0.52), (0.52,0.6)]
zlims = [(0.36, 0.6)]
sig = [0.0, 1.0, 2.0, 3.0]
# sig = [0.0]
nbins = [0.5, 1., 2.]
# lims = [1215.67, 1548.0 ]#Lya, CIV
lims = [1215.67]#Lya

lims_dict={1215.67 : 'lya',1548.0 : 'civ'}

m = 1

prefix ='/data1/hsalas/projects/QbC_mgii/saved_files/'

f = open("list.txt", 'w')

counter = 0
for lim in lims:
	for n in nbins:
		for s in sig:
			for zlim in zlims:
				for masslim in masslims:
					for ewlim in ewlims:
						for iplim in iplims:
							for dist in dist_types:
								for grid in grid_types:
									for sn in sn_types:
										for cl_type in cl_types:
											s_sn_str = '-s_{:.1f}_{}'.format(s, sn)
											grid_n_str =  '{}_n{:.1f}'.format(grid, n)
											ip_str = '-ip_{}_{:.1f}_to_{:.1f}'.format(dist, iplim[0], iplim[1])
											ew_str = '-rew_{:.1f}_to_{:.1f}'.format(ewlim[0], ewlim[1])
											mass_str = '-mass_10e{}_to_10e{}'.format(masslim[0], masslim[1])
											z_str = '-z_{:.2f}_to_{:.2f}'.format(zlim[0], zlim[1])
											lim_str = 'lim_{}'.format(lims_dict[lim])
											ist_str = 'dist-{}'.format(dist)
											dir_name = grid_n_str+'-'+lim_str+mass_str+ew_str+z_str+ip_str+s_sn_str
											command_str = '-minew {} -maxew {} -minmass {} -maxmass {} -minz {} -maxz {} -minip {} -maxip {} -d {} -s {} -gr {} -n {} -zct {}  -m {} -s2n {} -ls {}'.format(ewlim[0], ewlim[1], masslim[0], masslim[1], zlim[0], zlim[1], iplim[0], iplim[1], dist, s, grid, n, cl_type, m, sn, lim)
											print(dir_name)
											counter += 1
											# run QbC
											f.write(dir_name) 
											print(' Running QbC...')
											os.system("python qbc_mgii.py {}".format(command_str))							
f.close()
print("We have written {} results".format(counter))