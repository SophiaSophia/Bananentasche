#! /usr/bin/env python


__author__ = 'Sophia_Kronthaler_and_Tobias_Goeppel'

import random as ra
import time as ti
import numpy as np
import mytimer as mt


timer = mt.create(10)


class SpinLattice:

	
	def __init__(self, n=20, j=1., h=0., ordered=True):
		
		'''
		
        variable    |   description             
        ============|===========================
        n           |   number of spins         
        j           |   coupling strength       
        h           |   magnetic field          
        ordered     |   initial configuration   
		
		'''
		
		ra.seed(ti.time()) 
		
		self.n = n  
		self.j = j  
		self.ordered = ordered
		self.beta = None
		
		
		
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
		this functions simplifies life by returning the correct index
		# index is not out of range case
		# index out if range -> periodic boundary conditions
        # (python can deal with negative indices, so we do not have to implement the -1 case)
		'''
		if i < self.n: return i  
		
		else: return 0  
                        
		
	def getWeight(self):  # this function computes the Boltzmann weight factor needed in the Metropolis Algorithm
		
		self.weight = np.exp(-np.dot(self.deltaE, self.beta))  # this is the Boltzmann weight
		
		return None  # this is a void-like function
		
		
		

class SpinLattice_1d(SpinLattice):
	
	

	def __init__(self, n=20, j=1., h=0, ordered=True):
		
		SpinLattice.__init__(self, n=20, j=1., h=0, ordered=True)
		
		self.weight = None
		
		self.hDirection = bool(h >= 0)  	# direction of the h field (the >= ensures that all spins are up, if h==0
                                            # and if ordered==True. (we want the ordered state for h==0 to be
                                            # up, up, ..., up)
                                            
		self.config = self.getConfig(ordered)
		
		self.deltaE = 4*self.j  # h==0 implementation, the energy difference per bond is 2j
		
		

	def getConfig(self, ordered):  # True means up, False means down
		
		if ordered is True: return [self.hDirection for i in range(self.n)]  # all spins up
		
		else: return [ra.choice([True, False]) for i in range(self.n)]  # random spin config
	

	def getRandomSpin(self): 
		'''
			returns coordinate of a randomly chosen spin
		'''
		return np.random.randint(self.N-1) 
		
		
		
	def printSpinConfig(self):
		'''
			print all spins
			up-spins are represented by an *, down-spins by an |
		'''
		
		for i in range(self.n):
			
			print "*" if self.config[i] == True else "|", 
			
		print
		
		return None
		
		
		
	def getMag(self):
		
		Mag = 0
		
		for si in range(self.n):
			
			if self.config[si-1] is True: Mag += 1
			
			else:
				
				Mag -= 1
				
		return Mag
	
	
	def flipSpin_(self, si):  # this function performs the spin flip of the spin si
		
		self.config[si] = not self.config[si]  # True (up) becomes False (down) and False (down) becomes True (up)
		
		return None  # this is a void-like function		
		
		
	def getE(self):  # this function computes the total energy of the initial configuration
		
		E = 0
		
		for si in range(self.n):  # iterate over all bonds
			
			if self.config[si-1] == self.config[si]:  # up-up and down-down bonds lower the energy
				
				E -= self.j
				
			else:  # mixed bonds raise the energy
				
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
	
	def __init(self, n=20, j=1., h=0, ordered=True):
		
		SpinLattice_1d.__init__(self, n=20, j=1., h=0, ordered=True)
		
		self.h = abs(h)
		
		self.hStrong = bool(self.h > 2*self.j)
		'''
		we say that the field is strong if h > j. We have to distinguish  between "strong" and "weak" fields when
		checking in the Metropolis Algorithm weather the new state is accepted directly or not
		'''
	
		self.deltaE = self.getDeltaE_h()
	
		self.E = self.getE()
		'''
		compute the initial energy. We only compute the initial energy once in the beginning.
		during the simulation we only add or distract small amounts of energy (deltaE)
		'''
		
	def getE(self):
		
		'''
			calculates the total energy for a given temperature of the initial spin configuration, is called in changeT()
		'''

		E_h = self.getE_()
		
		for si in range(self.n):
			
			if self.config[si] == self.hDirection:
				
				E_h -= self.h
				
			else: E_h += self.h
			
		return E_h
		
	def getDeltaE(self): 
		
		'''TO SPEED UP THE SIMULATION WE TRY TO PERFORM AS LEES CALCULATION DURING THE LOOPS AS
           POSSIBLE. WE THUS CALCULATE THE DELTA ENERGIES BEFORE GOING INTO A LOOP
		'''
		
		deltaE_h = []
		
		deltaE_h.append(4*self.j + 2*self.h)
		deltaE_h.append(4*self.j - 2*self.h)
		deltaE_h.append(-4*self.j - 2*self.h)
		deltaE_h.append(-4*self.j + 2*self.h)
		deltaE_h.append(2*self.h)
		deltaE_h.append(- 2*self.h)
		return deltaE_h
		
	def getNewConfig(self, si):   # returns the new configuration due to the Metropolis algorithm
		
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
		
		SpinLattice.__init__(self, n=20, j=1., h=0, ordered=True)


		if h != 0:  # we only implement the h field free model in 2d
			
			print "We only implement the h field free model in 2d, try h=0."
			
			exit()  # shut down the program
		
		self.N = self.n * self.n  # total number of spins
		
		self.config = self.getConfig(ordered)
		
		self.weight = np.zeros(5)  # this are the Boltzmann weight factors. We initialize them with 0 since we do
		                           # not know the temperature yet.
		                           # TO SPEED UP THE SIMULATION WE MAKE USE OF NUMPY ARRAYS
		                           
		self.deltaE = self.getDeltaE()
		
	def getConfig(self, ordered): # returns initial spinlattice configuration
	
		# ordered = True: all spins up (True)
		# ordered = False: disordered initial configuration, spin direction is randomly
		# assigned
		                           
		config = np.ones((self.n, self.n), dtype=bool)  # TO SPEED UP THE SIMULATION WE MAKE USE OF NUMPY ARRAYS
		
		if ordered is True:
			
			return config
			
		else:
			
			for i in range(self.n):
				
				for j in range(self.n):
					
					
					config[i, j] = ra.choice([True, False])
					
			return config
			
			
	def getE(self): # calculates the total energy for a given temperature of the initial spin configuration is called in changeT()
		
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
		
		Mag = 0
		
		for si in range(self.n):
			
			for sj in range(self.n):
				
				if self.config[si, sj] == True: Mag += 1
				
				else: Mag -= 1
				
		return Mag





	def flipSpin(self, i, j):  # flips a specific spin with coordinates [i, j]

		self.config[i, j] = not self.config[i, j]
		
		return None



	def getDeltaE(self):  # TO SPEED UP THE SIMULATION WE TRY TO PERFORM AS LEES CALCULATION DURING THE LOOPS AS
							# POSSIBLE. WE THUS CALCULATE THE DELTA ENERGIES BEFORE GOING INTO A LOOP

		DeltaE_2 = np.array([8*self.j, 4*self.j, 0, -4*self.j, -8*self.j])  # to speed up the simulation we make use of
                                                                            # numpy arrays instead of python standard
                                                                            # arrays

		return DeltaE_2



	def getNewConfig(self, sN):  # choosing the random spin: 10% (function call getRandomSpin() )

		'''
		# IMPLEMENTATION OF THE METROPOLIS ALGORITHM: if there are more than one anti-parallel neighbours the energy
        # remains constant or is lowered by a spin flip and we can accept the flip directly
        # otherwise we have to roll the dice to decide whether we flip the spin or not
        '''

		#timer[0].start('getNewConfig_2d')
		
		#timer[1].start('divmod') # (7%)
		i, j = divmod(sN, self.n)  # to speed up the simulation we generate only one random number in the range [0, n*n]
                                   # and recover the spin coordinates by a simple modulo operation
        #timer[1].stop()

        #timer[2].start('case') # (4%)
		case = 4  # there are four possible nearest neighbour (NN) configurations, case is the number of anti-parallel
                  # neighbouring spins
        #timer[2].stop()

        #timer[3].start('which case?') # (80%)
		if self.config[i, j] == self.config[i, j-1]: case -= 1
		if self.config[i, j] == self.config[i-1, j]: case -= 1
		if self.config[i, j] == self.config[i, self.getIndex(j+1)]: case -= 1  # getIndex is not very costly
		if self.config[i, j] == self.config[self.getIndex(i+1), j]: case -= 1
		#timer[3].stop()
		
		if case >= 2:
			
			self.flipSpin(i, j)
			
			self.E += self.deltaE[case]
			
		else:
			
			if np.random.rand() <= self.weight[case]:
				self.flipSpin(i, j)
				
				self.E += self.deltaE[case]
				
			else: pass
			
		return None



	def printSpinConfig(self):  # prints the spin configuration (up spins = *, down spins = o)

		for i in range(self.n):
			
			for j in range(self.n):
				
				print "* " if self.config[i][j] == True else "o ",

			print

		return None



if __name__ == "__main__":



	spinLattice = SpinLattice_1d_h(h=-1.)
	spinLattice.printSpinConfig()


	'''
    spinLattice = SpinLattice( dimension=2)


    number_iterations = 100000
    temperature = 0.5

    print "Dimension: ", spinLattice.dimension
    print "Initial configuration (ordered = True, disordered = False): ", spinLattice.ordered
    print "Spinlattice size: ", spinLattice.n
    print "Temperature: ", temperature
    print "initial energy", spinLattice.E
    spinLattice.printSpinConfig()


    spinLattice.changeT(temperature)

    timer[9].start('number_interations')

    print "magnetisation: ", spinLattice.getMag()

    for i in range(number_iterations):

        spinLattice.getNewConfig(spinLattice.getRandomSpin())

    timer[9].stop()

    print "Number of itearions: ", number_iterations
    print "final energy", spinLattice.E
    spinLattice.printSpinConfig()

    #print "magnetisation: ", spinLattice.getMag()

    mt.table()
    
	'''    
