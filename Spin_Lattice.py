#! /usr/bin/env python

__author__ = 'Sophia_Kronthaler_and_Tobias_Goeppel'

import random as ra
import time as ti
import numpy as np
import mytimer as mt


#timer = mt.create(10)


class SpinLattice:
	
	
	
	def __init__(self, n=20, j=1., h=0., ordered=True):
		
		
		'''
        member variable    |   description             
        ===================|===========================
        n                  |   number of spins         
        j                  |   coupling strength                 
        ordered            |   initial configuration  
        beta               |   1./Temperature
        hDirection         |  for h>= 0 all spins are up if ordered=True
        h                  |   absolute value of the external magnetic field
        config             |   container for the initial configuration (getConfig())
        E                  |   container for the total energy
		'''
		
		ra.seed(ti.time()) 
		
		self.n = n
		self.j = j  
		self.ordered = ordered
		self.beta = None
		self.hDirection = bool(h >= 0)
		self.h = abs(h)
		self.config = self.getConfig(ordered)
		self.E = self.getE()
		
		
		
		
	def changeT(self, newT):
	
		'''
		to adjust the member variable beta according to the newT: beta = 1./newT
		in addition computes the new (temperature dependent) weights
		returns None
		'''
		
		self.beta = 1./newT
		
		self.getWeight()
		
		return None
		
		
	def getIndex(self, i):
		
		'''
		simplifies life by returning the correct index with regard to periodic
		boundary conditions (python can deal with negative indices,
		so we do not have to implement the -1 case)
		'''
		if i < self.n: return i  
		
		else: return 0  
                        
		
	def getWeight(self):
		'''
		computes the Boltzmann weight factor needed in the Metropolis Algorithm
		changes the member variable weight
		returns None
		'''
		self.weight = np.exp(-np.dot(self.deltaE, self.beta))
		
		return None
		
	def getRandomSpin(self): 
		'''
		returns coordinate of a randomly chosen spin: int(si)
		'''
		return np.random.randint(self.N-1)
		
		
	def test(self, number_iterations = 1000, temperature = 1., ):
		'''
		test routine
		number of iterations: number of spin flips
		prints out several control parameters
		'''
		
		print "Initial configuration (ordered = True, disordered = False): ", self.ordered
		print "Spinlattice size: ", self.n
		print "initial energy", self.E
		print "initial magnetisation: ", spinLattice.getMag()
		self.printSpinConfig()
		print "Temperature: ", temperature
		self.changeT(temperature)
		
		print "weights: ", self.weight
		print "energy difference", self.deltaE
		
		for i in range(number_iterations):
		
			self.getNewConfig(self.getRandomSpin())
		
		print "Number of itearions: ", number_iterations
		print "final energy", self.E
		print "final magnetisation: ", self.getMag()
		self.printSpinConfig()
		
		

class SpinLattice_1d(SpinLattice):
	'''
	for the one dimensional case(h=0)
	'''
	
	
	def __init__(self, n=20, j=1., h=0, ordered=True):
		
		'''
		member variable |    description
		=================|============================================================
		hStrong          |  field is strong if h > j (to distinguish 
		                 |  between "strong" and "weak" fields when
		                 |  checking in the Metropolis Algorithm weather
		                 |  the new state is accepted directly or not
		weight           |  container for Boltzmann weight factor, saved in the beginnig 
		                 |  (getWeight())
		deltaE           |  container for energy difference
		
		'''
		
		SpinLattice.__init__(self, n, j, h, ordered)
		
		self.N = n
		self.hStrong = bool(self.h > 2*self.j)
		self.weight = None
		self.deltaE = self.getDeltaE()
				
		
	def getDeltaE(self):
		
		return 4*self.j
		
	def getConfig(self, ordered):
		'''
		returns the initial configuration: array(bool)
		if ordered=True: all spins up
		if ordered=False: random spin config
		'''
		
		if ordered is True: return [self.hDirection for i in range(self.n)]
		
		else: return [ra.choice([True, False]) for i in range(self.n)]		
		
		
	def printSpinConfig(self):
		'''
			prints the curren spin configuration
			up-spins are represented by "*", down-spins by "|"
		'''
		
		for i in range(self.n):
			
			print "*" if self.config[i] == True else "|", 
			
		print
		
		return None
		
		
		
	def getMag(self):
		'''
		returns the magnetization of the current spin configuration:
		int(Mag)
		'''
		
		Mag = 0
		
		for si in range(self.n):
			
			if self.config[si-1] is True: Mag += 1
			
			else:
				
				Mag -= 1
				
		return Mag
	
	
	def flipSpin(self, si):
		'''
		performs the spin flip of the spin si
		returns None
		'''
		
		self.config[si] = not self.config[si]
		
		return None
		
		
	def getE(self):
		'''
		returns the total energy of the initial configuration
		parallel: -j
		anti-parallel: +j
		'''
		E = 0.
		
		for si in range(self.n):
			
			if self.config[si-1] == self.config[si]:
				
				E -= self.j
				
			else:
				
				E += self.j
				
		return E
        
        
	def getNewConfig(self, si):
		'''
		this function computes the energy difference due to the prospective flip of
		the spin si and flips the randomly chosen spin with if the metropolis algorithm
		condition is full filled
		'''
		
		if self.config[si-1] != self.config[self.getIndex(si+1)]:
			
			self.flipSpin(si)
			
		elif self.config[si-1] == self.config[self.getIndex(si+1)] and self.config[si] == self.config[si-1]:
				
			if ra.random() <= self.weight:
				
				self.flipSpin(si)
				
				self.E += self.deltaE
				
		else:
			
			self.flipSpin(si)
			
			self.E -= self.deltaE
			
		return None


	
class SpinLattice_1d_h(SpinLattice_1d):
	
	'''
	for the one dimensional case with an external magnetic field
	'''
		
	def getE(self):
		
		'''
			calculates the total energy for a given temperature of the
			initial spin configuration
			is called in changeT()
			only compute the initial energy once in the beginning, during
			the simulation we only add or distract small amounts of energy (deltaE)
		'''

		E = 0.
		
		for si in range(self.n):
			
			if self.config[si-1] == self.config[si]:
				
				E -= self.j
				
			else: E += self.j
			
		
		for si in range(self.n):
			
			if self.config[si] == self.hDirection:
				
				E -= self.h
				
			else: E += self.h
			
		return E
		
		
	def getDeltaE(self): 
		'''
		returns the energy difference of flipping one spin
		deltaE[0]:   4*j + 2*h
		deltaE[1]:   4*j - 2*h
		deltaE[2]: - 4*j - 2*h
		deltaE[3]: - 4*j + 2*h
		deltaE[4]:         2*h
		deltaE[5]:       - 2*h
		'''
		
		deltaE = np.array([4*self.j + 2*self.h, 4*self.j - 2*self.h, -4*self.j - 2*self.h, -4*self.j + 2*self.h, 2*self.h, - 2*self.h ]) 
			
		return deltaE
		
	def getNewConfig(self, si):
		'''
		returns the new configuration due to the Metropolis algorithm
		'''
		
		if self.config[si-1] == self.config[self.getIndex(si+1)]:
			
			if self.config[si] == self.config[si-1]:
				
				if self.config[si] == self.hDirection:
					
					if ra.random() <= self.weight[0]:
						
						self.flipSpin(si)
						
						self.E += self.deltaE[0]
						
				else:
					
					if self.hStrong is True:
						
						self.flipSpin(si)
						
						self.E += self.deltaE[1]
						
					elif ra.random() <= self.weight[1]:
						
						self.flipSpin(si)
						
						self.E += self.deltaE[1]
						
			else:
				
				if self.config[si] != self.hDirection:
					
					self.flipSpin(si)
					
					self.E += self.deltaE[2]
					
				else:
					
					if self.hStrong is False:
						
						self.flipSpin(si)
						
						self.E += self.deltaE[3]
						
					elif ra.random() <= self.weight[3]:
						
						self.flipSpin(si)
						
						self.E += self.deltaE[3]
						
		else:
			
			if self.config[si] == self.hDirection:
				
				if ra.random() <= self.weight[4]:
					
					self.flipSpin(si)
					
					self.E += self.deltaE[4]
					
			else:
				
				self.flipSpin(si)
				
				self.E += self.deltaE[5]




class SpinLattice_2d(SpinLattice):
	
	
	def __init__(self, n=20, j=1., h=0, ordered=True):
		

	
		SpinLattice.__init__(self, n, j, h, ordered)

		'''
		member variable   |   description
		==================|==========================================
		N                 |   total amount of spins: n*n
		weigth            |   container for Boltzmann weight factor
		                  |   saved in the beginnig, array of length(5)
		deltaE            |   container for energy difference,
		                  |   array of length(5), calculated in getDeltaE()
		'''



		if h != 0:
			
			print "We only implement the h field free model in 2d, try h=0."
			
			exit()
		
		self.N = self.n * self.n
		self.weight = np.zeros(5)
		self.deltaE = self.getDeltaE()
		
	def getConfig(self, ordered):
		'''
		returns the initial configuration: array(bool)
		if ordered=True: all spins up
		if ordered=False: random spin config
		
		'''
		                           
		config = np.ones((self.n, self.n), dtype=bool)
		
		if ordered is True:
			
			return config
			
		else:
			
			for i in range(self.n):
				
				for j in range(self.n):
					
					
					config[i, j] = ra.choice([True, False])
					
			return config
			
			
	def getE(self):
		'''
		returns the total energy of the initial configuration
		parallel: -j
		anti-parallel: +j
		'''
		
		E = 0   # for parallel spins += -j, for anti-parallel spins += j
		
		for si in range(self.n):
			
			for sj in range(self.n):
				
				if self.config[si, sj-1] == self.config[si, sj]:
					
					E -= self.j
					
				if self.config[si, sj-1] != self.config[si, sj]:
					
					E += self.j
					
				if self.config[sj-1, si] == self.config[sj, si]:
					
					E -= self.j
					
				if self.config[sj-1, si] != self.config[sj, si]:
					
					E += self.j
					
		return E
		
	
	def getMag(self):
		'''
		returns the magnetization of the current spin configuration:
		int(Mag)
		'''
		
		Mag = 0
		
		for si in range(self.n):
			
			for sj in range(self.n):
				
				if self.config[si, sj] == True: Mag += 1
				
				else: Mag -= 1
				
		return Mag





	def flipSpin(self, i, j):
		'''
		performs the spin flip of the spin si
		returns None
		'''

		self.config[i, j] = not self.config[i, j]
		
		return None



	def getDeltaE(self): 
		'''
		returns the energy difference of flipping one spin
		DeltaE[0]=   8*j
		DeltaE[1]=   4*j
		DeltaE[2]=   0
		DeltaE[3]= - 4*j
		DeltaE[4]= - 8*j
		'''

		DeltaE = np.array([8*self.j, 4*self.j, 0, -4*self.j, -8*self.j])  
		
		return DeltaE



	def getNewConfig(self, sN):
		
		'''
		returns the new configuration due to the Metropolis algorithm
		'''
			
		i, j = divmod(sN, self.n)  
		
		case = 4 
		if self.config[i, j] == self.config[i, j-1]: case -= 1
		if self.config[i, j] == self.config[i-1, j]: case -= 1
		if self.config[i, j] == self.config[i, self.getIndex(j+1)]: case -= 1
		if self.config[i, j] == self.config[self.getIndex(i+1), j]: case -= 1
		
		
		if case >= 2:
			
			self.flipSpin(i, j)
			
			self.E += self.deltaE[case]
			
		else:
			
			if np.random.rand() <= self.weight[case]:
				self.flipSpin(i, j)
				
				self.E += self.deltaE[case]
				
			else: pass
			
		return None



	def printSpinConfig(self):
		'''
		prints the curren spin configuration
		up-spins are represented by "*", down-spins by "o"
		'''

		for i in range(self.n):
			
			for j in range(self.n):
				
				print "* " if self.config[i][j] == True else "o ",

			print

		return None

		
		


if __name__ == "__main__":

	
	spinLattice = SpinLattice_2d()
	
	spinLattice.test(temperature=3.)
	
