import numpy as np
import matplotlib.pyplot as plt
import Spin_Lattice as spl
import os as os
import time as ti



class Ising:
	
	def __init__(self, dimension = 1, n = 200, j=1., h=0., ordered = True):
		
	
		self.perSpin = 1./n**dimension
		
		self.T_Liste = None
		self.dimension = dimension
		self.n = n
		self.j = j
		self.h = h
		self.ordered = ordered
		self.sweep = None
		self.num_sweeps = None
		
		self.Mag_MeanList = np.zeros(0)
		self.Mag_VarList = np.zeros(0)
		self.E_MeanList = np.zeros(0)
		self.susz = np.zeros(0)
		self.capacity = np.zeros(0)
		self.filename = None
		
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
	
		self.sweep = sweep
		self.num_sweeps = num_sweeps
		
	
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
		
		return None
	
	

	def getObservables(self):
		
		
		self.E_MeanList *= self.perSpin
		
		self.Mag_MeanList *= self.perSpin

		self.susz = self.Mag_VarList / self.T_List * self.perSpin
		
		for i in range(self.T_List.size-1):
			
			self.capacity = np. append(self.capacity, abs(self.E_MeanList[i] - self.E_MeanList[i+1])/ abs(self.T_List[i] - self.T_List[i+1]))
			
		return None

		
	def saveObservables(self):
	
		self.filename = "./Results_"+ti.strftime("%d.%m.%Y_%H.%M")
	
		os.mkdir(self.filename)
		
		
		file = open (self.filename + "/E_MeanList", "w")
		np.savetxt(file, self.E_MeanList, fmt='%.10f')
		file.close()
	
		file = open (self.filename + "/Mag_MeanList", "w")
		np.savetxt(file, self.Mag_MeanList, fmt='%.10f')
		file.close()

		file = open (self.filename + "/Suszeptibility", "w")
		np.savetxt(file, self.susz, fmt='%.10f')
		file.close()
		
		file = open (self.filename + "/Heat_Capacity", "w")
		np.savetxt(file, self.capacity, fmt='%.10f')
		file.close()
		
		file = open (self.filename + "/Temperature", "w")
		np.savetxt(file, self.T_List, fmt='%.10f')
		file.close()
		
		file = open (self.filename + "/SpinLattice_Settings", "a")
		file.write("Dimension: " + str(self.dimension) + "\n")
		file.write("Lattice size: " + str(self.n) + "\n")
		file.write("Coupling constant: " + str(self.j) + "\n")
		file.write("External magnetic field: " + str(self.h) + "\n")
		file.write("Starting with an ordered configuration: " + str(self.ordered) + "\n")
		file.write("Sweep size: " + str(self.sweep) + "\n")
		file.write("Number of sweeps: " + str(self.num_sweeps) + "\n")
		
	def savePlots(self):
	

		plt.plot(self.T_List, self.Mag_MeanList, 'ro')
		plt.xlabel('$Temperature$ $t$')
		plt.ylabel('$Magnetisation$ $per$ $spin$ $m$')
		plt.savefig(self.filename + "/Mag_Mean")
		plt.close()
		
		plt.plot(self.T_List, self.E_MeanList, 'ro')
		plt.xlabel('$Temperature$ $t$')
		plt.ylabel('$Energy$ $per$ $spin$ $e$')
		plt.savefig(self.filename + "/E_Mean")
		plt.close()
		
		plt.plot(self.T_List, self.susz,'ro')
		plt.xlabel('$Temperature$ $t$')
		plt.ylabel('$Suszeptibility$ $per$ $spin$ $\chi $')
		plt.savefig(self.filename + "/Suszeptibility")
		plt.close()

		plt.plot(self.T_List[:-1], self.capacity, 'ro')
		plt.xlabel('$Temperature$ $t$')
		plt.ylabel('$Heat$ $capacity$ $per$ $spin$ $c$')
		plt.savefig(self.filename + "/Heat_Capacitiy")
		plt.close()
		
		
if __name__ == "__main__":
	
	ising = Ising(dimension=1, n=200)
	ising.setTemperatureRange(1,2,0.5)
	ising.simulate(1000, 100)
	
	ising.getObservables()
	
	ising.saveObservables()
	ising.savePlots()


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
