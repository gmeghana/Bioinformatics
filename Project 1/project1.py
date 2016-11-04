from math import *
import sys

# safe_log function is needed to represent zero probability
def safe_log(x):
	if x==0:
		# python's interpretation of neginf is -1.7976931348623157e+308
		neginf=-1*sys.float_info.max
		return neginf
	else:
		return log10(x)

# LogSumP is useful when adding probabilities in log representation
def LogSumP(X, Y):
	if X >= Y: 
		return X + log10(1 + 10**(Y - X))
	else:
		return Y + log10(1 + 10**(X - Y))

# this function calculates the likelihood of the observations give different values of theta
def calculate_basic_score(prob_obs_unmut, prob_obs_mut, theta):
	return LogSumP(prob_obs_mut + safe_log(theta), prob_obs_unmut + safe_log(1 - theta))

# this function calculates nCr which is needed in the binomial expression
def nCr(n,r):
    f = factorial
    return f(n) / f(r) / f(n-r)

# this function calculates P(obs i|k,theta) from a specific person I
def calculate_prob_individual(prob_obs_unmut, prob_obs_mut, theta):
	sum = 0
	# loop through possible copy number of a specific person (k=0,1,2)
	for k in range(0,3):
		# use the equation for LogP(obs|theta) given in the basic variant scoring model, except replace theta with k/2
		# multiply binomial expression 
		sum = sum + nCr(2,k) * (theta**k) * (1-theta)**(2-k) * LogSumP(prob_obs_mut + safe_log(k/2.0), prob_obs_unmut + safe_log(1-k/2.0))
	return sum

# opens a file from stdin
file = open(sys.argv[1], 'r')

basic_sum_2 = 0.0
basic_sum_0 = 0.0

# matrix contains subset of observations i from a specific person I for all possible values of theta
# set the max range of the number of individuals to 1000 according to professor
matrix = [[0 for x in range(101)] for x in range(1000)]

# contains number of times each individual appears in a file input  
count_individual = {}

L_vector = [0] * 101
Q_vector = [0] * 101

# parse through each line in the file input and separate id, prob_obs_unmut, and prob_obs_mut
for line in file:
	words = line.split()
	id = int(words[0])
	prob_obs_unmut = float(words[1])
	prob_obs_mut = float(words[2])

	# calculate the sum of all likelihoods given theta=0.2
	basic_sum_2 = basic_sum_2 + calculate_basic_score(prob_obs_unmut, prob_obs_mut, 0.2)

	# calculate the sum of all likelihoods given theta=0.0
	basic_sum_0 = basic_sum_0 + calculate_basic_score(prob_obs_unmut, prob_obs_mut, 0.0)

	# for each person, loop through all possible values of theta and fill matrix with sum of all observations for a given person
	for j in range(0,101):
		 matrix[id][j] = matrix[id][j] + calculate_prob_individual(prob_obs_unmut, prob_obs_mut, j/100.0)

# compute the log probabilities and store in L_vector		
for k in range(0,101):
	for i in range (0,1000):
		# sum over all observations for a given theta and store in L_vector
		L_vector[k] += matrix[i][k]

# set the first value of the Q_vector as neginf
Q_vector[0] = safe_log(0)
for m in range(1,101):
	# need to calculate the value of Q_vector[100] so that it can be used to calculate the advanced score
	Q_vector[m] = LogSumP(Q_vector[m-1], LogSumP(L_vector[m-1],L_vector[m]) + log10((1/100.0)/2.0))

# calculate the basic score after calculating sum of all likelihoods
basic_score = basic_sum_2 - basic_sum_0
print "BASIC_SCORE ", basic_score

# calculate the advanced score
advanced_score = Q_vector[100] - L_vector[0]
print "ADVANCED_SCORE ", advanced_score

# calculate the log of the posterior probability density for theta = j/m (j ranges from 0 to 101 and m=100)
for j in range(0,101):
	print "LOG_P_PHI", j/100.0, L_vector[j] - Q_vector[100]
file.close()