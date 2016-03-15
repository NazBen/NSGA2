import numpy as np
import matplotlib.pyplot as plt
import random
from testFunctions import testFunc1, testFunc2

def computeBeta(u, eta=1):
    if u < 0.5:
        return (2.*u)**(1./(eta + 1.))
    else:
        return 1. / ((2.*(1. - u) ) ** (1./(eta + 1.)))


def computeDelta(r, eta=1):
    """
    TODO : check for the value of eta
    """
    if r < 0.5:
        return (2 * r)**(1 / (eta + 1)) - 1
    else:
        return 1 - (2 * (1 - r))**(1 / (eta + 1))


class Individual(object):
    def __init__(self, objective_functions, input_range):
        """
        An individual conserve the information about a gen. His input values, output values through the objective functions 
        """
        self.setObjectiveFunctions(objective_functions)
        self.setInputRange(input_range)

        self.rank = np.inf
        self.distance = 0.
        self._func_eval_count = 0

    def setInputValues(self, values):
        """
        """
        if type(values) in [float, int]:
            assert self._input_dimension == 1, "Input dimension problem"
        else:
            dim = len(values)
            assert dim == self._input_dimension, "Input dimension problem"

        self._input_values = values

        if any(self._input_values < self._lower_bounds) or any(self._input_values > self._upper_bounds):
            print "Out of Band"

    def __str__(self):
        string = "Input value: %s\n" % (str(self.input_value))
        string += "Output values: %s\n" % (str(self.output_values))
        return string

    def setInputRange(self, input_range):
        """
        Set the dimension of the input variables [min, max] * dim
        """
        if type(input_range) is not np.array:
            input_range = np.array(input_range)

        errorMessage = "The input range has no lower or upper bounds"
        try:
            dim, p = input_range.shape
            assert p == 2, errorMessage
        except:
            dim = 1
            assert len(input_range) == 2, errorMessage

        self._input_dimension = dim
        self._lower_bounds = input_range[:, 0]
        self._upper_bounds = input_range[:, 1]
        self._input_range = input_range

    def setObjectiveFunctions(self, objective_functions):
        """
        Set the objective functions and save the number.
        """
        if type(objective_functions) in [list, tuple]:
            num_objectives = len(objective_functions)
        else:
            raise TypeError("Objective functions are not in a list or a tuple")

        self._objective_functions = objective_functions
        self._num_objectives = num_objectives

    def evaluate(self):
        """
        Evaluate the current input_values with the objective functions
        """
        self._output_values = np.empty((self._num_objectives, ))

        for i, func in enumerate(self._objective_functions):
            self._output_values[i] = func(self._input_values)

        self._func_eval_count += self._num_objectives

    def init_sort(self):
        """
        Init of the number of diminated indiv
        """
        self.dominated_indivs = [] # Init le list of dominant individuals
        self.domination_count = 0

    def __rshift__(self, other):
        """
        self indiv dominates
        """               
        # Arrays of boolean to see which indivs are greater, lower or equals
        equal = self._output_values == other._output_values
        greater = self._output_values > other._output_values

        # If at least one obj functions of the current indiv is greater than the
        # obj of the other indiv, then the current indiv does not 
        # dominated the other.
        if any(greater):
            return False
        # If there is no obj functions of the current indiv that are 
        # greater than the obj functions of the other indiv
        # and the objectives are not all equals, then the current indiv
        # dominate the other
        elif not any(greater) and ( sum(equal) <= 1 ):
            return True
        else:
            return False
    def __lshift__(self, other):
        """
        """
        return other >> self


class Population(list):
    def __init__(self, pop_size=0):
        """
        A population is a list of croomosones, but with specific methods to sort them, and replace them easily
        """

        self.setPopSize(pop_size)
        self._func_eval_count = 0

    def initialize(self):
        """
        """
        # Clear the list
        del self[:]

        # Create the population
        self._create_population()

        # Initialize the population
        self._initialize_population()

    def _create_population(self):
        """
        Create the population of individual
        """
        for _ in range(self._pop_size):
            # We create the element and add it in the population object
            indiv = Individual(self._objective_functions, self._input_range)  # Construction of the object
            self.append(indiv)
    
    def _initialize_population(self):
        """
        """
        for indiv in self:
            rand = np.random.random(self._input_dimension)
            value = self._lower_bounds + (self._upper_bounds - self._lower_bounds)*rand
            indiv.setInputValues(value)

    def evaluate(self):
        """
        """
        for indiv in self:  # For each individual
            indiv.evaluate()  # Evaluate it
            # Increment the evaluation cost
            self._func_eval_count += indiv._func_eval_count
            
    def nondominated_sort(self):
        """
        """
        # Initializing the front
        self.front = [[]]

        iFront = 0 # Index of the front
        for indiv_i in self:  # For each individual
            # We need to check if this individial is dominant
            # For that we compare this one with every other individual.
            pop_except_i = self[:]  # A copy of the population
            pop_except_i.remove(indiv_i)  # We remove the current indiv
            indiv_i.init_sort()  # Initialise the sorting

            # For each indiv except the i-th one
            for indiv_j in pop_except_i:
                # If individual i dominates j
                if indiv_i >> indiv_j:
                    # Then i dominates j
                    # We append it in the list
                    indiv_i.dominated_indivs.append(indiv_j)

                # i is dominant to j if
                elif indiv_i << indiv_j :
                    # If at least one func leads to j < i and 
                    indiv_i.domination_count += 1

            # indiv i is not dominated by any other
            # i is ont the front
            if indiv_i.domination_count == 0:
                indiv_i.rank = iFront
                self.front[iFront].append(indiv_i)
                indiv_i.iFront = iFront

        # While the front is not empty
        while self.front[iFront]:
            Q = [] # Subset of individuals
            # For each indiv in the front
            for indiv_i in self.front[iFront]:
                # For each indiv j dominated by i
                for indiv_j in indiv_i.dominated_indivs:
                    indiv_j.domination_count -= 1

                    if indiv_j.domination_count == 0:
                        indiv_j.rank = iFront + 1
                        Q.append(indiv_j)

            iFront += 1
            self.front.append(Q)
            for indiv in Q:
                indiv.iFront = iFront

    def computeCrowdingDistance(self):
        """
        """       
        for k, iFront in enumerate(self.front):  # For each front
            if iFront:
                numIndivInFront = len(iFront)
                front_output_values = np.empty((numIndivInFront, self._num_objectives))
                
                # For each indiv in the front, get the objective values
                for i, indiv in enumerate(iFront):
                    indiv.distance = 0.
                    front_output_values[i, :] = indiv._output_values

                # The indiv in the front are sorted
                front_sortedID = front_output_values.argsort(axis=0)
                for m in range(self._num_objectives):
                    # Border distances setted to inf
                    iFront[front_sortedID[0, m]].d = np.inf
                    iFront[front_sortedID[-1, m]].d = np.inf

                    # Get the min and max of objfuncs
                    fmin = iFront[front_sortedID[0, m]]._output_values[m]
                    fmax = iFront[front_sortedID[-1, m]]._output_values[m]

                    sorted = front_sortedID[1:-1, m]
                    # For each indiv in the front, compute the distance
                    # Ugly but it works...
                    # TODO: change this
                    for k, ID in enumerate(sorted):
                        iFront[ID].distance += (iFront[front_sortedID[k+2, m]]._output_values[m] - iFront[front_sortedID[k, m]]._output_values[m]) / (fmax - fmin)
                        iFront[ID].distanceComputed = True

    def setObjectiveFunctions(self, objective_functions):
        """
        Set the objective functions and save the number.
        """
        if type(objective_functions) in [list, tuple]:
            num_objectives = len(objective_functions)
        else:
            raise TypeError("Objective functions are not in a list or a tuple")

        self._objective_functions = objective_functions
        self._num_objectives = num_objectives

    def setInputRange(self, input_range):
        """
        Set the dimension of the input variables [min, max] * dim
        """
        if type(input_range) is not np.array:
            input_range = np.array(input_range)

        errorMessage = "The input range has no lower or upper bounds"
        try:
            dim, p = input_range.shape
            assert p == 2, errorMessage
        except:
            dim = 1
            assert len(input_range) == 2, errorMessage

        self._input_dimension = dim
        self._lower_bounds = input_range[:, 0]
        self._upper_bounds = input_range[:, 1]
        self._input_range = input_range

    def setPopSize(self, pop_size):
        """
        Set the population size
        """
        assert type(pop_size) is int, "pop_size must be integer. Given type is %s" % (type(pop_size))
        self._pop_size = pop_size

    def getOutputValues(self):
        """
        """
        outputValues = np.empty((self._pop_size, self._num_objectives))

        for i in range(self._pop_size):
            outputValues[i, :] = self[i]._output_values
        return outputValues

    def getParetoFront(self):
        """
        """
        front = self.front[0]
        frontSize = len(front)

        pareto = np.empty((frontSize, self._num_objectives))
        for i, indinv in enumerate(front):
            pareto[i, :] = indinv._output_values
        return pareto

    def clear(self):
        """
        """
        for indiv in self:
            self.remove(indiv)
    

class NSGA2:
    def __init__(self, 
                 objective_functions, 
                 input_range,
                 pop_size=100, 
                 num_generation=10, 
                 prob_cross=0.9,
                 pool_size=50,
                 tour_size=10,
                 eta_c=1.,
                 eta_m=1.,
                 history=True):
        """

        """
        self._pop_size = pop_size
        self._pool_size = pool_size
        self._tour_size = tour_size
        self._prob_cross = prob_cross
        self._input_range = input_range
        self._num_generation = num_generation
        self._history = history
        self._eta_c = eta_c
        self._eta_m = eta_m

        self._num_objectives = len(objective_functions)
        self._objective_functions = objective_functions

        self.initialize()

        if self._history:
            self._population._history = history
            self._previouspops = Population()

    def initialize(self):
        """

        """
        # Create the population
        self._population = Population(self._pop_size)

        # Set the input range
        self._population.setInputRange(self._input_range)

        # Set the objective functions
        self._population.setObjectiveFunctions(self._objective_functions)

        # Initialize the population
        self._population.initialize()

    def run(self, withPlot=False):
        """

        """
        # Evaluate the initiale population
        self._population.evaluate()

        # TODO : add more conditions
        iLoop = 0
        # While the stop condition are not respected
        while not (iLoop == self._num_generation) :
            print "loop : ", iLoop

            # The population is sorted
            self._population.nondominated_sort()

            # The crowding distances are computed for each individual
            self._population.computeCrowdingDistance()

            # Select the parents
            self._tournament_selection()

            # Cross and mutate the parents
            self._genetic_operator()
            
            # The new population is sorted again
            self._population.nondominated_sort()

            # Same for crowding distances
            self._population.computeCrowdingDistance()

            # We actualize of the population
            self._actualizePopulation()
                       
            if withPlot:
                self.plotFunctionalSpace()

            iLoop += 1

    def _actualizePopulation(self):
        """
        """
        newPopulation = self._population
        del newPopulation[:]

        placeLeft = self._pop_size
        for iFront in self._population.front:
            numInFront = len(iFront)

            if numInFront < placeLeft:
                newPopulation.extend(iFront)
                placeLeft -= numInFront
            else:                
                # Sort considering the distance
                distances = [indiv.distance for indiv in iFront]
                sorted = np.argsort(distances, axis=0)[::-1]

                # Append only the individual in the remain space
                for i in sorted[:placeLeft]:
                    newPopulation.append(iFront[i])
                
                # There is no more space, we stop
                break
                               
        if self._history:
            for indiv in self._population:
                if indiv not in self._previouspops:
                    self._previouspops.append(indiv)

        self._population = newPopulation

    def _cross_individuals(self, indiv_i, indiv_j):
        """
        """
        u = np.random.rand()
        # Compute beta
        beta = computeBeta(u, self._eta_c)

        # Value of the parents
        p1 = indiv_i._input_values
        p2 = indiv_j._input_values

        children = []
        # For each child (there is 2 children)
        for i in range(2):
            # Create the child
            child = Individual(self._objective_functions, self._input_range)
            # We set the input value
            child._input_values = 0.5*((1 - beta)*p1 + (1 + beta)*p2) 
            child.evaluate()  # We evaluate
            children.append(child)  # We append to the list
            # The second got -beta
            beta *= -1

        return children
        
    def _mutate_individual(self, indiv_i):
        """
        """
        u = np.random.rand()
        # Compute Delta
        delta = computeDelta(u, self._eta_m)

        # Compute the new value
        value = indiv_i._input_values + (indiv_i._upper_bounds - indiv_i._lower_bounds)*delta

        # Create the child
        child = Individual(self._objective_functions, self._input_range)

        # We set the input value
        child._input_values = value 
        child.evaluate() # We evaluate

        return child

    def _genetic_operator(self):
        """
        Here we want to improve our population. The selection of best individuals should already have
        been done. Now we want to mutate or cross these individuals.
        A fixed propability is fixed for selecting the the individuals are crossed or mutated.
        """
        # Randomly select parents
        first_parent = np.random.choice(self._parents, self._pool_size)
                
        for indiv_i in first_parent:
            # We check with a certain probability if we cross them or
            if np.random.rand() <= self._prob_cross:  # We cross them                
                # Tmp list of individual
                tmp = self._parents[:]
                # We delete the 1st parent
                tmp.remove(indiv_i)
                # Randomly choiced 2nd parent
                indiv_j = np.random.choice(tmp)
                # We add them in the population
                self._population.extend(self._cross_individuals(indiv_i, indiv_j))
            else:
                # We add them in the population
                self._population.append(self._mutate_individual(indiv_i))

    def _tournament_selection(self):
        """
        When the individuals are sorted and their ranks and distances are computed, we can select them.
        The selection is done by a "tournament". The pool_size represent the number of individuals selected
        for the next step (mutation or cross). Each selected elements is done by randomly selecting tour_size
        individuals. By comparing these individuals we select the best one. Which is the one with the lower rank
        or the higher distance (when the ranks are equals).
        """
        # List of selected elements
        selected = []

        # The range of candidate (used for the random selecting only)
        listCandidate = range(self._pop_size) # List of the candidate

        # Start of the loop
        for i in range(self._pool_size):  # For each space in the pool
            # We select randomly a fixed number of candidates for the pre-selection
            candidates = random.sample(self._population, self._tour_size)

            # Initialisation of the ranks and distances vectors
            ranks = np.zeros(self._tour_size, dtype=int)
            distances = np.zeros(self._tour_size, dtype=np.float)

            # We get the rank and distance informations
            for j, candidate in enumerate(candidates):  # For each candidate
                ranks[j] = candidate.rank
                distances[j] = candidate.distance

            # We look for the index of the one (or multiple) minimum ranks.
            mins = np.argmin(ranks)

            if mins.size > 1:
                # There is multiple minimum, we selecte using the distances
                # We compare the distance and get the element id with max distance
                mins = np.where(distances == np.min(distances[mins]))[0]
                # The max distance win

            # He won
            selected.append(candidates[mins]) 
                
        self._parents = selected

    def plotFunctionalSpace(self):
        """

        """
        fig, ax = plt.subplots()
        paretoFront = self._population.getParetoFront()
        out = self._population.getOutputValues()
        #deleted = self._population.deletedPop.getOutputValues()
        ax.plot(out[:, 0], out[:, 1], 'b.')
        #ax.plot(deleted[:, 0], deleted[:, 1], 'k.')
        ax.plot(paretoFront[:, 0], paretoFront[:, 1], 'o')
        ax.set_xlabel("f1(x)") 
        ax.set_ylabel("f2(x)")
        ax.grid()
        if self._objective_functions == testFunc2:
            ax.set_xlim(-20., -14.)
            ax.set_ylim(-12., 2.)
        
        fig.tight_layout()
        plt.show()


np.random.seed(0)
test = 2
if test == 1:
    list_func = testFunc1  # objective functions
    dim = 1 # problem dimension
    xmin, xmax = -10.**3, 10.**3
elif test == 2:
    list_func = testFunc2  # objective functions
    dim = 3 # problem dimension
    xmin, xmax = -5., 5.

# range of the inputs
input_range = np.array([[xmin, xmax]] * dim)

# parameters
pop_size = 100 # population size
num_generation = 50
nsga = NSGA2(list_func, input_range, pop_size, num_generation, pool_size=50, tour_size=10, prob_cross=0.9, eta_c=0.5)
nsga.run(withPlot=False)
nsga.plotFunctionalSpace()

#u = np.linspace(0., 1., 100)
#beta = computeBeta(u, 10)
#fig, ax = plt.subplots()
#ax.plot(u, beta)