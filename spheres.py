'''
in terminal:

python spheres.py

'''
import random 
import math
import scipy
from scipy import stats

#return radius (in angstroms) of a sphere representing a protein
# of randomly generated chain length betwen 80-500 residues
def return_sphere(): 
	l = random.randint(80, 500)
	w = l*110
	r = 10*(0.066*(w**(1./3)))
	return l,r

#generate a set of cartesian coordinates for a random point on the
#surface of a sphere of radius r
pi = math.pi
def return_coords(r):	
	u = random.random()
	v = random.random()
	theta = math.acos(2*v-1)
	phi = 2*pi*u
	x = r*((math.sin(theta))*(math.cos(phi)))
	y = r*((math.sin(theta))*(math.sin(phi)))
	z = r*(math.cos(theta))
	coordinates = (x,y,z)
	return coordinates
	
#set of functions to get the distance between 2 points (norm of
#displacement vector)
def displace(p1,p2):
	x = p1[0] - p2[0]
	y = p1[1] - p2[1]
	z = p1[2] - p2[2]
	return (x,y,z)	
def norm(x): 
    return math.sqrt(sum(i**2 for i in x))		
def dist(p1,p2):
	v = displace(p1,p2)
	return norm(v)
	
#check distances between a given point and list of points to ensure
#they are all greater than 2.5 angstroms	
def check_dists(p, list):
	dists = []
	try:
		for i in list:
			d = dist(p,i)
			dists.append(d)
	except:
		pass
	clash = [x for x in dists if x < 2.50]	
	if len(clash)>0:
		return False
	else:
		return True						
		
#take a list of coordinates and assign each point as negative or 
#positive based on an input probability
#for selecting a negative point (R =1.212 meso, 1.130 thermo;
#P = 0.5479 meso, 0.5305 thermo)	
def assign_q(list, P_neg):
	neg = []
	pos = []
	for p in list:
		x = random.random()
		if x <= P_neg:
			neg.append(p)
		else:
			pos.append(p)
	return neg, pos			
	
#function to mimic analysis of charge distributions in proteins
#identifies 'attractive' and 'repulsive' interactions between points
#using input cutoff distance (5 and 7 A) and reports the 
#number of interaxns, branching, ratio of rep to atr,
#and the fraction of isolated points 
def properties(coords, neg, pos, d_cut):
	pairs = []
	anti_pairs = []
	for n in neg:
		for p in pos:
			if dist(n,p) <= d_cut:
				pairs.append((n,p))
		for i in neg:
			if 0 < dist(n,i) <= d_cut:
				anti_pairs.append((n,p))
	n_ip = len(pairs)
	n_ap = float(len(anti_pairs))
	try:
		ratio_rep_atr = n_ap/n_ip
	except:
		ratio_rep_atr = 0.0	
	n_points = float(len(coords))
	branched = []
	for i in coords:
		count = []
		for n,p in pairs:
			if i == n:
				count.append(n)
			elif i == p:
				count.append(p)
		n_interaxns = len(count)
		if n_interaxns > 1:
			branched.append(i)
	n_branched = float(len(branched))
	ip_points = []
	for a,b in pairs:
		if a not in ip_points:
			ip_points.append(a)
		if b not in ip_points:
			ip_points.append(b)			
	n_ip_points = len(ip_points)	
	try:
		frac_branched = n_branched/n_ip_points
	except:
		frac_branched = 0.0	
	n_iso = float(n_points - n_ip_points)
	frac_iso = n_iso/n_points
	return n_points,n_ip, ratio_rep_atr, frac_iso, frac_branched	
			
#######################################
meso_nq = []; thermo_nq = []
meso_nips = []; thermo_nips = []
meso_ratios = []; thermo_ratios = []
meso_isos = []; thermo_isos = []
meso_branches = []; thermo_branches = []
for x in range(10000):
	coords = []
	l,r = return_sphere()
	nqm_nres = random.gauss(0.234, 0.034)
	nqm = int(round(nqm_nres*l))
	while len(coords) < nqm:
			p = return_coords(r)
			if check_dists(p,coords) == True:
				coords.append(p)
			else:
				continue
	neg,pos = assign_q(coords, 0.5)
	np, n_ip, ratio_rep_atr, frac_iso, frac_branched = properties(coords, neg, pos, 5.0) 
	norm_nq = np/float(l)
	norm_nip = n_ip/float(l)
	meso_nq.append(norm_nq)
	meso_nips.append(norm_nip);meso_ratios.append(ratio_rep_atr)
	meso_isos.append(frac_iso);meso_branches.append(frac_branched)
	coords = [] ##########now repeat for therm dist on same sphere
	nqt_nres = random.gauss(0.279, 0.033)
	nqt = int(round(nqt_nres*l))
	while len(coords) < nqt:
			p = return_coords(r)
			if check_dists(p,coords) == True:
				coords.append(p)
			else:
				continue
	neg,pos = assign_q(coords, 0.5)
	np, n_ip, ratio_rep_atr, frac_iso, frac_branched = properties(coords, neg, pos, 5.0) 
	norm_nip = n_ip/float(l)
	norm_nq = np/float(l)
	thermo_nq.append(norm_nq)
	thermo_nips.append(norm_nip);thermo_ratios.append(ratio_rep_atr)
	thermo_isos.append(frac_iso);thermo_branches.append(frac_branched)	
	
####function to get stats on lists of values						
def stats(meso_list,thermo_list):
	meso_array = scipy.array(meso_list)
	thermo_array = scipy.array(thermo_list)
	meso_mean = scipy.mean(meso_array)
	meso_std = scipy.std(meso_array)
	thermo_mean = scipy.mean(thermo_array)
	thermo_std = scipy.std(thermo_array)
	p_val = scipy.stats.ttest_ind(thermo_array,meso_array)[1]
	s = '\nMm: ' + str(meso_mean) + '\nSTDm: ' + str(meso_std) + '\nMt: ' + str(thermo_mean) + '\nSTDt: ' + str(thermo_std) + '\nP: ' + str(p_val)
	return s
	
######get the stats and write to outfile
ofile = open('results5A-10k-P5050-thesisnorm-wnq.txt','w')
nip_stats = stats(meso_nips, thermo_nips)
s1 = 'Normalized N_ip: ' + nip_stats + '\n'
ofile.write(s1)
ratio_stats = stats(meso_ratios, thermo_ratios)
s2 = '\nRatio rep atr: ' + ratio_stats + '\n'
ofile.write(s2)
iso_stats = stats(meso_isos, thermo_isos)
s3 = '\nFraction isolated: ' + iso_stats + '\n'
ofile.write(s3)
branched_stats = stats(meso_branches, thermo_branches)
s4 = '\nFraction branched: ' + branched_stats + '\n'
ofile.write(s4)
meso_nip_data = str(meso_nips);thermo_nip_data = str(thermo_nips)
meso_ratio_data = str(meso_ratios);thermo_ratio_data = str(thermo_ratios)
meso_iso_data = str(meso_isos);thermo_iso_data = str(thermo_isos)
meso_branched_data = str(meso_branches);thermo_branched_data = str(thermo_branches)
s5 = '\nN_ip data: ' + '\nMESO: ' + meso_nip_data + '\nTHERMO: ' + thermo_nip_data
s6 = '\nRatio data: ' + '\nMESO: ' + meso_ratio_data + '\nTHERMO: ' + thermo_ratio_data
s7 = '\niso data: ' + '\nMESO: ' + meso_iso_data + '\nTHERMO: ' + thermo_iso_data
s8 = '\nbranched data: ' + '\nMESO: ' + meso_branched_data + '\nTHERMO: ' + thermo_branched_data
ofile.write(s5)
ofile.write(s6)
ofile.write(s7)
ofile.write(s8)
nq_stats = stats(meso_nq,thermo_nq)
s9 = '\nq over nres: ' + nq_stats + '\n'
ofile.write(s9)
ofile.close()


			

'''
to make a 3d scatter plot of the points:

from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
l,r = return_sphere()
coords = []
for i in range(1000):
	coords.append(return_coords(r))
xs = []
ys =[]
zs = []
for x,y,z in coords:
	xs.append(x)
	ys.append(y)
	zs.append(z)
fig = pyplot.figure()
ax = Axes3D(fig)	
ax.scatter(xs,ys,zs)
pyplot.show()



fig.savefig('plot.png')	
'''
	
	