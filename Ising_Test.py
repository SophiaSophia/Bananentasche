import numpy as np
import matplotlib.pyplot as plt
import Spin_Lattice as spl



class Ising:
	
	def __init__(self, dimension = 1, n = 20, j=1., h=0., ordered = True):
		
	
		self.perSpin = 1./n**dimension
		
		self.T_Liste = None
		
		self.Mag_MeanList = np.zeros(0)
		self.Mag_VarList = np.zeros(0)
		self.E_MeanList = np.zeros(0)
		self.susz = np.zeros(0)
		self.capacity = np.zeros(0)
		
		if dimension == 1 and h== 0.:
			
			self.spinlattice = spl.SpinLattice_1d(n, j, h, ordered)
		
		if dimension == 1 and h!= 0.:
			
			self.spinlattice = spl.SpinLattice_1d_h(n, j, h, ordered)
		
		if dimension == 2:
			
			self.spinlattice = spl.SpinLattice_2d(n, j, h, ordered)
		
		

	def setTemperatureRange(self, a, b, steps ):
	
		self.T_List = np.arange(0.1, 4, 0.05)
		
		return None
		

	def simulate(self,sweep = 1000, num_sweeps = 50):
	
		for t in range(self.T_List.size):

			self.spinlattice.changeT(self.T_List[t])

			print "temperature step", t

			MagList = np.zeros(num_sweeps)
			EList = np.zeros(num_sweeps)

			for i in range(num_sweeps):

				for j in range(sweep):

					self.spinlattice.getNewConfig(self.spinlattice.getRandomSpin())

				MagList[i] = self.spinlattice.getMag()
				EList[i] = self.spinlattice.E

			self.Mag_MeanList = np.	append(self.Mag_MeanList, np.mean(MagList))
			self.Mag_VarList = np.append(self.Mag_VarList, np.var(MagList))
			self.E_MeanList = np.append(self.E_MeanList, np.mean(EList))
		
		print self.Mag_MeanList
		print self.Mag_VarList
		print self.E_MeanList
		
		return None
	
	

	def getObservables(self):
		
		
		self.E_MeanList *= self.perSpin
		
		self.Mag_MeanList *= self.perSpin

		self.susz = self.Mag_VarList / self.T_List * self.perSpin
		
		for i in range(self.T_List.size):
			
			self.capacity = np. append(abs(self.E_MeanList[i] - self.E_MeanList[i+1])/ abs(self.T_List[i] - self.T_List[i+1])
			
		return None

		
	#def saveObservables(self):
		
	#	file = open ("E_MeanList", "w")
	#	for item in thelist:
	#		print>>file, item


if __name__ == "__main__":
	
	ising = Ising()
	ising.setTemperatureRange(1,2,0.5)
	ising.simulate(100,10)
	
	ising.getObservables()



'''

#print Mag_MeanList, T_List

plt.plot(T_List, Mag_MeanList, 'ro')
plt.show()


plt.plot(T_List, E_MeanList, 'ro')
plt.show()


#print Mag_VarList

plt.plot(T_List, susz,'ro')
plt.show()

plt.plot(T_List, capacity, 'ro')
plt.show()
'''
