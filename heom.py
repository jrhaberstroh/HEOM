#!/usr/bin/env python
import csv
import numpy as np
import math
import parameters as p
import cPickle
import os

def HowMany(N):
	return math.factorial(N + 7)/(math.factorial(N) * math.factorial(7))

def Choose(a,b):
	return math.factorial(a)/(math.factorial(a-b) * math.factorial(b))

def Bucket(a,b):
	if (a < 0):
		return 0
	if (b < 1):
		return 0
	return math.factorial(a+b-1)/(math.factorial(a) * math.factorial(b-1))

def SOx(A,B):
	return A*B - B*A

def SOo(A,B):
	return A*B + B*A

def ThetaSO(Vj,B):
	return 2.*p.kt*p.reorg*SOx(Vj,B) - 1.0j*p.reorg*p.damp*SOo(Vj,B)

def Index(n_vec):
	output_index = 0
	i = 0
	n_sites = len(n_vec)
	remain = n_vec.sum()
	
	#The first subspace is for the total number
	# This is an especially useful subspace because it allows us to
	# back out the number of excitations from the index
	output_index += Bucket(remain - 1, n_sites + 1)
	while (remain > 0):
		#The next subspaces are for each index, done sequentially
		site_val = n_vec[i]
		while (site_val > 0):
			output_index += Bucket(remain, n_sites - i - 1) 
			#We do not use an extra bucket because we are in a number-restricted subspace
			site_val -= 1
			remain -= 1
		i += 1
	return output_index 

def REC_VectorLessThanN(N, n_sum, n_vec, i, nvec_list, n_list, nplu_list, nmin_list):
	if (i == len(n_vec)):
		#Base case
		ind = len(n_list)
		nvec_list.append(n_vec.copy())
		n_list.append(Index(n_vec))
		nplu_list.append(np.zeros(0,'complex'))
		nmin_list.append(np.zeros(0,'complex'))
		#print "n vector:",
		#print n_vec
		if (n_vec.sum() < N):
			for j in range(len(n_vec)):
				n_vec[j] += 1
				nplu_list[ind] = np.append(nplu_list[ind], Index(n_vec) + j*1j)
				#print n_vec,
				n_vec[j] -= 2
				if (n_vec[j] >= 0):
					nmin_list[ind] = np.append(nmin_list[ind], Index(n_vec) + j*1j)
					#print n_vec,
				#print 
				n_vec[j] += 1
			#print nplu_list[ind]

	if (i < len(n_vec)):
		#Recursion
		REC_VectorLessThanN(N,n_sum,n_vec,i+1, nvec_list, n_list, nplu_list, nmin_list)
		while (n_sum < N):
			n_sum += 1
			n_vec[i] += 1
			REC_VectorLessThanN(N,n_sum,n_vec,i+1, nvec_list, n_list, nplu_list, nmin_list)
		n_vec[i] = 0

def VectorLessThanN(N, n_vec, nvec_list, nind_list, nplu_list, nmin_list):
	REC_VectorLessThanN(N, sum(n_vec), n_vec, 0, nvec_list, nind_list, nplu_list, nmin_list)
	print "Largest Index: " 
	print max(nind_list)
	print "Number of Indexes: " 
	print len(nind_list)
	print "Indexes after culling duplicates: "
	print len(set(nind_list))


def HEOMRate(state, dt_st, nind_list, nvec_list, nplu_list, nmin_list, numMatx):
	for matx in range(numMatx):
		#Reset the dt matrix:
		dt_st[nind_list[matx]] = -1j*SOx(H_e, state[nind_list[matx]])
		if nvec_list[matx].sum() < N:
			for i in range(7):
				dt_st[nind_list[matx]] -= nvec_list[matx][i] * p.damp * state[nind_list[matx]]

		for plu_ind in nplu_list[matx]:
			# The imaginary component of the plus index stores the j value
			# The real component stores the index
			dt_st[nind_list[matx]] += 1j * SOx( V_j[int(plu_ind.imag)], state[int(plu_ind.real)] )
							
		for min_ind in nmin_list[matx]:
			dt_st[nind_list[matx]] += 1j * nvec_list[matx][int(min_ind.imag)] * ThetaSO( V_j[int(min_ind.imag)], state[int(min_ind.real)] )



############################################################

reader=csv.reader(open("H_electronic.csv","rb"),delimiter=',')
x=list(reader)
H_e=np.matrix(x).astype('float')
dt = 1./5309. #in spectroscopic units of cm for .5fs timestep 
n_time = 500   #number of timesteps
N = 10
numMatx= Bucket(N,7+1)  #7+1 buckets because there is the "null" bucket
cutoff = Bucket(N-1,7+1) #The cutoff index at and above which HEOM are truncated
print "Number of operators: ",
print numMatx
n_vec = np.array([0,0,0,0,0,0,0])
d_ind = Index(n_vec)
print "Index of the density matrix: ",
print d_ind


state = list()
statek1 = list()
statek2 = list()
statek3 = list()
k1 = list()
k2 = list()
k3 = list()
k4 = list()
for count in range(numMatx):
	state.append( np.matrix(np.zeros((7,7), 'complex')) )
	statek1.append( np.matrix(np.zeros((7,7), 'complex')) )
	statek2.append( np.matrix(np.zeros((7,7), 'complex')) )
	statek3.append( np.matrix(np.zeros((7,7), 'complex')) )
	k1.append( np.matrix(np.zeros((7,7), 'complex')) )
	k2.append( np.matrix(np.zeros((7,7), 'complex')) )
	k3.append( np.matrix(np.zeros((7,7), 'complex')) )
	k4.append( np.matrix(np.zeros((7,7), 'complex')) )

state[0][0,0] = 1 
print state[0]



# pickle data to skip steps
vector_pickle = "n_init.pkl"
index_pickle = "nind_init.pkl"
neighborplus_pickle = "nind_plus_init.pkl"
neighborminus_pickle = "nind_minus_init.pkl"
projection_pickle = "projection_init.pkl"
result_pickle = "result.pkl"



if os.path.isfile(vector_pickle) and \
   os.path.isfile(index_pickle) and \
   os.path.isfile(neighborplus_pickle) and \
   os.path.isfile(neighborminus_pickle) and \
   os.path.isfile(projection_pickle):
	print "Skipped initialization by loading pickles. Delete pickles to disable this!"
	f_vector = open(vector_pickle,'r')		
	f_index = open(index_pickle,'r')		
	f_plus = open(neighborplus_pickle,'r')
	f_minus= open(neighborminus_pickle,'r')
	f_proj = open(projection_pickle,'r')
	nvec_list = cPickle.load(f_vector)
	nind_list = cPickle.load(f_index)
	nplu_list = cPickle.load(f_plus)
	nmin_list = cPickle.load(f_minus)
	V_j       = cPickle.load(f_proj)



else:
	print "Computing the neighbor list..."
	nvec_list = list()
	nind_list = list()
	nplu_list = list()
	nmin_list = list()
	VectorLessThanN(N, n_vec, nvec_list, nind_list, nplu_list, nmin_list)
	print HowMany(10)
	print len(nind_list)
	print len(nplu_list)
	print len(nmin_list)
	
	print "Computing the projection operators..."
	
	V_j = list()
	for i in range(7):
		entry = np.matrix( np.zeros((7,7) , 'complex'))
		entry[i,i] = 1
		V_j.append(entry)
	
	
	print "Initialization complete. "

	f_vector = open(vector_pickle,'w')		
	f_index = open(index_pickle,'w')		
	f_plus = open(neighborplus_pickle,'w')
	f_minus= open(neighborminus_pickle,'w')
	f_proj = open(projection_pickle,'w')

	cPickle.dump(nvec_list, f_vector)
	cPickle.dump(nind_list, f_index)
	cPickle.dump(nplu_list, f_plus)
	cPickle.dump(nmin_list, f_minus)
	cPickle.dump(V_j, f_proj)

	print "Initialization Pickled, too!"
	




print "Beginning the main run with ",
print n_time,
print " steps of ",
print dt * 5309.,
print "fs."

traj = list()

	
for t_n in range(n_time):
	print "Step ",
	print t_n
	
	HEOMRate(state, k1, nind_list, nvec_list, nplu_list, nmin_list, numMatx)
	for matx in range(numMatx):
		statek1[nind_list[matx]] = state[nind_list[matx]] + (.5 * k1[nind_list[matx]] * dt)

	HEOMRate(statek1, k2, nind_list, nvec_list, nplu_list, nmin_list, numMatx)
	for matx in range(numMatx):
		statek2[nind_list[matx]] = state[nind_list[matx]] + (.5 * k2[nind_list[matx]] * dt)
	
	HEOMRate(statek2, k3, nind_list, nvec_list, nplu_list, nmin_list, numMatx)
	for matx in range(numMatx):
		statek3[nind_list[matx]] = state[nind_list[matx]] + (k3[nind_list[matx]] * dt)

	HEOMRate(statek3, k4, nind_list, nvec_list, nplu_list, nmin_list, numMatx)
	for matx in range(numMatx):
		state[nind_list[matx]] += 1./6. * dt * (k1[nind_list[matx]] + 2. * k2[nind_list[matx]] + 2. * k3[nind_list[matx]] + k4[nind_list[matx]])

	traj.append(state[0].copy())
	print state[0]

f_result = open(result_pickle,'w')
cPickle.dump(traj,f_result)
print state[0]







