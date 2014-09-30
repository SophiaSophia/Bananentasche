import Spin_Lattice as sl
import numpy as np
import matplotlib.pyplot as plt


n = 100
num_sweeps = 50
sweep = 100000
perSpin = 1./n

T_List = np.arange(0.1, 4, 0.05)
Mag_MeanList = np.zeros(T_List.size)
Mag_VarList = np.zeros(T_List.size)
E_MeanList = np.zeros(T_List.size)


#spinlattice.printSpinConfig()
spinlattice = sl.SpinLattice(n=n, ordered=True, dimension=2)

for t in range(T_List.size):

    spinlattice.changeT(T_List[t])

    print t

    MagList = np.zeros(num_sweeps)
    EList = np.zeros(num_sweeps)

    for i in range(num_sweeps):

        for j in range(sweep):

            spinlattice.getNewConfig(spinlattice.getRandomSpin())

        MagList[i] = spinlattice.getMag()
        EList[i] = spinlattice.E

    Mag_MeanList[t] = np.mean(MagList)
    Mag_VarList[t] = np.var(MagList)
    E_MeanList[t] = np.mean(EList)

susz = Mag_VarList / T_List * perSpin

Mag_MeanList *= perSpin

E_MeanList *= perSpin

capacity = np.zeros(T_List.size)

for i in range(T_List.size-1):

    capacity[i] = abs(E_MeanList[i] - E_MeanList[i+1])/ abs(T_List[i] - T_List[i+1])


capacity[T_List.size-1] = 1.

print
print capacity.size
print T_List.size
print

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
