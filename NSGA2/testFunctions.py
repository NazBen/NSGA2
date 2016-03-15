import numpy as np

def testFunc1_1(x):
    return x**2

def testFunc1_2(x):
    return (x - 2.)**2

def testFunc2_1(x):
    """
    """
    return (-10*np.exp(-0.2*np.sqrt(x[:-1]**2 + x[1:]**2 ))).sum()

def testFunc2_2(x):
    return (abs(x)**0.8 + 5*np.sin(x**3)).sum()

testFunc1 = [testFunc1_1, testFunc1_2]
testFunc2 = [testFunc2_1, testFunc2_2]
