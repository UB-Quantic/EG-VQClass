import numpy as np
import matplotlib.pyplot as plt

####################################
# Loss function and Target Function
####################################

def Cp(point, label, angles):
    """Computes the cost value of a single input.
    
    Args.
        point (dim=2 float): coordinates of input point.
        label (int): what class does point belong to.
        angles (array float): values of free parameters.

    Ret.
        cp (float): cost value for a single input.
    """
    output = circ(point, angles)
    cp = np.linalg.norm(label/3 - f(output))
    cp = 0.5*cp*cp
    return cp

def C(data, angles):
    """Computes the cost function over the whole data set.

    Args.
        data (array float): set of input points.
        angles (array float): values of free parameters.
    Ret.
        c (float): cost value for the entire set.
    """
    c = 0
    n = len(data[0])
    for i in range(n):
        c+=Cp(data[0][i], data[1][i], angles)
    c/=n
    return c

def cost(angles):
    """Computes the cost over the training data set.

    Args.
        angles (array float): values of free parameters.
        
    Ret.
        ret (float): cost value for the entire set.
    """
    ret = C(training_data, angles)
    return ret

def f(state):
    """Computes the target function for a given output.

    Args.
        state (dim=2 complex): final state of quantum system.

    Ret.
        p0 (float): probability of the state being |0>.
    """
    p0 = state[0] * np.conj(state[0])
    return p0
