
import random 
import math

###################
def return_sphere(): 
	l = random.randint(80, 500)
	w = l*110
	r = 10*(0.066*(w**(1./3)))
	return r, l

###########################
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
	
########################
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
						
		
#t#####################
def assign_q(p, P_neg):
	x = random.random()
	if x <= P_neg:
		return 'Negative'
	else:
		return 'Positive'
		
	
#f#############################
def properties(neg, pos, d_cut):
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
	return ratio_rep_atr	
			
#####################variation of R with increasing rho
outfile = 'Rvsrho_5kpts_2.txt'
f = open(outfile,'w')
for x in range(4):
	data_points = []
	r = return_sphere()
	A = 4*pi*r**2
	s1 = '\nRadius: ' + str(r) + '\nArea: ' + str(A) 
	f.write(s1)
	neg = []; pos = []
	for i in range(500):
		for k in range(10):
			p = return_coords(r)
			if assign_q(p, 0.5) == 'Negative':
				neg.append(p)
			else:
				pos.append(p)
		R = properties(neg, pos, 5.0)
		nq = float(len(neg)+len(pos))
		rho = nq/A
		data_points.append((R,rho))
	s2 = '\nData (R,rho): ' + str(data_points) + '\n'
	f.write(s2)
f.close()	

####################distribution of final R values
outfile = '1kpts.txt'
f = open(outfile,'w')
data_points = []
r = 18.
A = 4*pi*r**2
nq = 1000
rho = nq/A
s1 = '\nRadius: ' + str(r) + '\nArea: ' + str(A) + '\nRho: ' + str(rho)
f.write(s1)
for i in range(500):
	neg = []; pos = []
	for k in range(1000):
		p = return_coords(r)
		if assign_q(p, 0.5) == 'Negative':
			neg.append(p)
		else:
			pos.append(p)
	R = properties(neg, pos, 5.0)
	data_points.append(R)
s2 = '\ndat = ' + str(data_points) + '\n'
f.write(s2)
f.close()	



'''		
import matplotlib.pyplot as plt
x = []
y = []
for (R, rho) in data_points:
	x.append(rho)
	y.append(R)				
plt.plot(x, y, linewidth=2.0)
plt.ylabel('R')
plt.xlabel('Point Density')
plt.grid(True)
plt.show()
	
import scipy		
dat = []
datarray = scipy.array(dat)
hist = plt.hist(datarray, 20, alpha=0.5, color='b', label='2kpts')
plt.legend(loc=1)
plt.ylabel('Frequencies')
plt.grid(True)
plt.xlabel('R')
plt.show()
mean = scipy.mean(datarray)
std = scipy.std(datarray)
print mean, std
plt.savefig('.png')

f = open('5.txt','w')
for (a,b) in d5:
	s = '\n' + str(a) + '	' + str(b)
	f.write(s)
f.close()

thermarray = scipy.array(thermdat)
mesarray = scipy.array(mesdat)
thermophiles = plt.hist(thermarray, 20, alpha=0.5, color='b', label='56 pts')
mesophiles = plt.hist(mesarray, 20, alpha=0.5, color='g', label='46 pts')
plt.legend(loc=1)
plt.ylabel('Frequencies')
plt.grid(True)
plt.xlabel('R')
plt.show()
plt.set_title('R values for mean number of pts in mesophiles and thermophiles')
plt.savefig('.png')


'''
	

		


	

	
	