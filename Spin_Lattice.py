__author__ = 'Sophia_Kronthaler_and_Tobias_Goeppel'

import random as ra
import time as ti
import numpy as np
import mytimer as mt





timer = mt.create(10)


class SpinLattice:



    def __init__(self, n=10, j=1., h=0, ordered=True, dimension=1):

        ra.seed(ti.time())  # seed the random number generator

        self.n = n  # number of spins

        self.j = j  # coupling strength

        self.dimension = dimension

        self.ordered = ordered

        self.beta = None  # the temperature will be set by the Metropolis Algorithm


        if dimension == 1:  # one dimensional case

            self.N = self.n * 1  # total number of spins. This becomes important in the 2d model

            self.weight = None  # this is the Boltzmann weight factor. We initialize it with None since we do not know
                                # the temperature yet

            self.config = self.getConfig(ordered)  # get the initial configuration

            self.printSpinConfig = self.printSpinConfig

            self.flipSpin = self.flipSpin_

            self.hDirection = bool(h >= 0)  # direction of the h field (the >= ensures that all spins are up, if h==0
                                            # and if ordered==True. (we want the ordered state for h==0 to be
                                            # up, up, ..., up)

            if h == 0:  # use the function implementations for the h==0 case. (Function names end with an _)

               self.getE = self.getE_  # use h==0 implementation

               self.deltaE = 4*self.j  # h==0 implementation, the energy difference per bond is 2j

               self.getNewConfig = self.getNewConfig_  # use h==0 implementation


            else:  # use the function implementations for the h!=0 case. (Function names end with an _h)

                self.h = abs(h)  # the absolute value of h

                self.hStrong = bool(self.h > 2*self.j)  # we say that the field is strong if h > j. We have to
                                                        # distinguish  between "strong" and "weak" fields when
                                                        # checking in the Metropolis Algorithm weather the new state is
                                                        # accepted directly or not

                self.getE = self.getE_h  # use h!=0 implementation

                self.deltaE = self.getDeltaE_h()  # use h!=0 implementation

                self.getNewConfig = self.getNewConfig_h  # use h!=0 implementation


        if dimension == 2:

            if h != 0:  # we only implement the h field free model in 2d

                print "We only implement the h field free model in 2d, try h=0."

                exit()  # shut down the program


            self.N = self.n * self.n  # total number of spins

            self.config = self.getConfig_2d(ordered)

            self.getE = self.getE_2d

            self.weight = np.zeros(5)  # this are the Boltzmann weight factors. We initialize them with 0 since we do
                                       # not know the temperature yet.
                                       # TO SPEED UP THE SIMULATION WE MAKE USE OF NUMPY ARRAYS

            self.deltaE = self.getDeltaE_2d()

            self.getNewConfig = self.getNewConfig_2d

            self.printSpinConfig = self.printSpinConfig_2d

            self.flipSpin = self.flipSpin_2d


        if dimension > 2:  # we only implement the one d and the two d model

            print "d > 2 not implemented, try d=1 or d=2"

            exit()  # shut down the program


        self.E = self.getE()  # compute the initial energy. We only compute the initial energy once in the beginning.
                              # during the simulation we only add or distract small amounts of energy (deltaE)



    def getIndex(self, i):  # this functions simplifies life by returning the correct index

        if i < self.n: return i  # index is not out of range case

        else: return 0  # index out if range -> periodic boundary conditions
                        # (python can deal with negative indices, so we do not have to implement the -1 case)



    def getConfig(self, ordered):  # True means up, False means down

        if ordered is True: return [self.hDirection for i in range(self.n)]  # all spins up

        else: return [ra.choice([True, False]) for i in range(self.n)]  # random spin config



    def changeT(self, newT):  # get a new Temperature

        self.beta = 1./newT  # change the member variable T to the new T

        self.getWeight()  # compute the new (temperature dependent) weights
                          # TO SPEED UP THE SIMULATION WE TRY TO PERFORM AS LEES CALCULATION DURING THE LOOPS AS
                          # POSSIBLE. WE THUS CALCULATE THE WEIGHTS BEFORE GOING INTO A LOOP

        return None  # this is a void-like function



    def getWeight(self):  # this function computes the Boltzmann weight factor needed in the Metropolis Algorithm

        self.weight = np.exp(-np.dot(self.deltaE, self.beta))  # this is the Boltzmann weight

        return None  # this is a void-like function



    def getRandomSpin(self): return np.random.randint(spinLattice.N-1)  # get coordinate of a randomly chosen spin
                                                                        # NUMPY RANDOM NUMBER GENERATOR IS TWICE AS
                                                                        # FAST AS THE MATH RANDOM NUMBER GENERATOR!


    def printSpinConfig(self):  # this function is only used for debugging

        for i in range(self.n):  # print all spins

            print "*" if self.config[i] == True else "|",  # up-spins are represented by an *, down-spins by an |

        print  # go th a new line

        return None  # this is a void-like function




###################################################################################
########################### member functions for d == 1 ###########################
###################################################################################


    def flipSpin_(self, si):  # this function performs the spin flip of the spin si

        self.config[si] = not self.config[si]  # True (up) becomes False (down) and False (down) becomes True (up)

        return None  # this is a void-like function




#########################################
###### member functions for h == 0 ######
#########################################



    def getE_(self):  # this function computes the total energy of the initial configuration

        E = 0

        for si in range(self.n):  # iterate over all bonds

            if self.config[si-1] == self.config[si]:  # up-up and down-down bonds lower the energy

                E -= self.j

            else:  # mixed bonds raise the energy
                E += self.j

        return E



    def getNewConfig_(self, si):  # this function computes the energy difference due to the prospective flip of
                                  # the spin si and flips the randomly chosen spin with if the metropolis algorithm
                                  # condition is full filled

        if self.config[si-1] != self.config[self.getIndex(si+1)]:  # if the left and the right neighbouring spin
                                                                   # have different directions, we directly accept the
                                                                   # flip of th spin in the middle

            self.flipSpin(si)  # perform the spin flip

        elif self.config[si-1] == self.config[self.getIndex(si+1)] and self.config[si] == self.config[si-1]:

            if ra.random() <= self.weight:  # if the considered spin and its two neighbours have the same direction
                                            # a spin flip will raise the energy. According to the Metropolis Algorithm
                                            # we only accept the spin flip with a certain probability given by the
                                            # Boltzmann weight. We thus draw a random number and compare it to the
                                            # corresponding Boltzmann weight

                self.flipSpin(si)  # perform the spin flip

                self.E += self.deltaE  # adjust the energy

        else:

            self.flipSpin(si)  # if the neighbours to the left and to right do have the same direction and if the
                               # considered spin in the middle does have the opposite direction, we directly accept
                               # the spin flip since the energy will be lowered.

            self.E -= self.deltaE  # adjust the energy

        return None  # this is a void-like function




#########################################
###### member functions for h != 0 ######
#########################################



    def getE_h(self):  # calculates the total energy for a given temperature of the initial spin configuration
                       # is called in changeT()

        E_h = self.getE_()

        for si in range(self.n):

            if self.config[si] == self.hDirection:
                E_h -= self.h

            else: E_h += self.h

        return E_h



    def getDeltaE_h(self):  # TO SPEED UP THE SIMULATION WE TRY TO PERFORM AS LEES CALCULATION DURING THE LOOPS AS
                            # POSSIBLE. WE THUS CALCULATE THE DELTA ENERGIES BEFORE GOING INTO A LOOP

        deltaE_h = []

        deltaE_h.append(4*self.j + 2*self.h)
        deltaE_h.append(4*self.j - 2*self.h)
        deltaE_h.append(-4*self.j - 2*self.h)
        deltaE_h.append(-4*self.j + 2*self.h)
        deltaE_h.append(2*self.h)
        deltaE_h.append(- 2*self.h)

        return deltaE_h



    def getNewConfig_h(self, si):   # returns the new configuration due to the Metropolis algorithm

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




###################################################################################
########################### member functions for d == 2 ###########################
###################################################################################



    def getConfig_2d(self, ordered):  # returns initial spinlattice configuration
                                      # ordered = True: all spins up (True)
                                      # ordered = False: disordered initial configuration, spin direction is randomly
                                      #                  assigned

        config = np.ones((self.n, self.n), dtype=bool)  # TO SPEED UP THE SIMULATION WE MAKE USE OF NUMPY ARRAYS

        if ordered is True: return config

        else:

            for i in range(self.n):

                for j in range(self.n):

                    config[i, j] = ra.choice([True, False])

            return config



    def getE_2d(self):  # calculates the total energy for a given temperature of the initial spin configuration
                        # is called in changeT()

        E = 0   # for parallel spins += -j
                # for anti-parallel spins += j

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



    def flipSpin_2d(self, i, j):  # flips a specific spin with coordinates [i, j]

        self.config[i, j] = not self.config[i, j]

        return None



    def getDeltaE_2d(self):  # TO SPEED UP THE SIMULATION WE TRY TO PERFORM AS LEES CALCULATION DURING THE LOOPS AS
                             # POSSIBLE. WE THUS CALCULATE THE DELTA ENERGIES BEFORE GOING INTO A LOOP

        DeltaE_2 = np.array([8*self.j, 4*self.j, 0, -4*self.j, -8*self.j])  # to speed up the simulation we make use of
                                                                            # numpy arrays instead of python standard
                                                                            # arrays

        return DeltaE_2



    def getNewConfig_2d(self, sN):  # choosing the random spin: 10% (function call getRandomSpin() )

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


        # IMPLEMENTATION OF THE METROPOLIS ALGORITHM: if there are more than one anti-parallel neighbours the energy
        # remains constant or is lowered by a spin flip and we can accept the flip directly
        # otherwise we have to roll the dice to decide whether we flip the spin or not

        if case >= 2:

            #timer[4].start('flipSpin') # (2%)
            self.flipSpin(i, j)  # flipping the spin, in both cases
            # #timer[4].stop()

            #timer[7].start('self.E += self.deltaE[case]') # (2%)
            self.E += self.deltaE[case]  # getting the change in energy, in both cases
            #timer[7].stop()

        else:

            if np.random.rand() <= self.weight[case]:

                #timer[4].start('flipSpin') # (1%)
                self.flipSpin(i, j)
                #timer[4].stop()

                #timer[7].start('self.E += self.deltaE[case]') # (2%)
                self.E += self.deltaE[case]
                #timer[7].stop()

            else: pass


        #timer[0].stop()

        return None



    def printSpinConfig_2d(self):  # prints the spin configuration (up spins = *, down spins = o)

        for i in range(self.n):

            for j in range(self.n):

                print "* " if self.config[i][j] == True else "o ",

            print

        return None



if __name__ == "__main__":



    spinLattice = SpinLattice(j= 1., n=20, ordered=False, dimension=2)


    number_interations = 100000
    temperature = 0.5

    print "Dimension: ", spinLattice.dimension
    print "Initial configuration (ordered = True, disordered = False): ", spinLattice.ordered
    print "Spinlattice size: ", spinLattice.n
    print "Temperature: ", temperature
    print "initial energy", spinLattice.E
    spinLattice.printSpinConfig()


    spinLattice.changeT(temperature)

    for i in range(number_interations):

        spinLattice.getNewConfig(spinLattice.getRandomSpin())

    print "Number of itearions: ", number_interations
    print "final energy", spinLattice.E
    spinLattice.printSpinConfig()

    mt.table()