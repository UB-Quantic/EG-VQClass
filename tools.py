import numpy as np
import matplotlib.pyplot as plt

def accuracy(data, angles):
    """Computes the accuracy of the classification for a set of free parameters.
    
    Args.
        data (array float): set of points as pairs coordinates-label.
        angles (array float): values of free parameters.
        
    Ret.
        accuracy value.
    """
    results = [[np.argmax(f(circ(x, angles))), y]
           for (x, y) in zip(data[0], data[1])]
    return sum(int(x == y) for (x, y) in results)

def test (data, angles):
    """Organise the data in arrays that can be plotted.
    
    Args.
        data (array float): set of points as pairs coordinates-label.
        angles (array float): values of free parameters.
        
    Ret.
        circles (array of arrays): points organised according to which circle they belong to.
        right (array float): points which were correctly taged.
        wrong (array float): points which were incorrectly taged.
    """
    
    circles = [[],[],[],[]]
    right, wrong = [], []
    
    results = [[np.argmax(f(circ(x, angles))), y]
               for (x, y) in zip(data[0], data[1])]
    
    for (x, y) in zip(data[0], results):
        circles[y[0]].append(x)
        if(y[0] == y[1]):
            right.append(x)
        else:
            wrong.append(x)
            
    return circles, right, wrong 

def plot (data, angles):
    """Creates two scatter plots, one depicting the classification itself, the other with the accuracy.
    
    Args.
        data (array float): set of points as pairs coordinates-label.
        angles (array float): values of free parameters.
    """
    circles, right, wrong = test(data, angles)
    cx = [
        0.,
        -1.,
        0.4,
        0.4
    ]

    cy = [
        0.,
        0.,
        0.4,
        -0.6
    ]

    r = [
        2.,
        0.8,
        0.4,
        0.6
    ]
    
    plt.figure(figsize = (9, 4))

    x = [[],[],[],[]]
    y = [[],[],[],[]]
    z = [[],[]]
    t = [[],[]]

    for i in range(len(circles)):
        for point in circles[i]:
            x[i].append(point[0])
            y[i].append(point[1])


    ax = plt.subplot(121)
    circle = []
    for i in range(4):
        circle.append(plt.Circle((cx[i], cy[i]), r[i], color='k', fill=False))
        ax.add_artist(circle[-1])
    ax.plot(x[0], y[0], 'bo', x[1], y[1], 'ko', x[2], y[2], 'mo', x[3], y[3], 'yo', markersize = 3)
    plt.xlabel('x')
    plt.ylabel('y')


    for rights in right:
        z[0].append(rights[0])
        z[1].append(rights[1])

    for wrongs in wrong:
        t[0].append(wrongs[0])
        t[1].append(wrongs[1])

    bx = plt.subplot(122)
    for i in range(4):
        circle.append(plt.Circle((cx[i], cy[i]), r[i], color='k', fill=False))
        bx.add_artist(circle[-1])
    bx.plot(z[0], z[1], 'go', t[0], t[1], 'ro', markersize = 3) 
    plt.xlabel('x')

    plt.suptitle('sucess rate {:2.2f}%'.format(len(right) / 10))

############################################
# Loss function & Target Function
############################################

def Cp(point, label, angles):
    """Computes the cost value of a single input.
    
    Args.
        point (dim=2 float): coordinates of input point.
        label (int): label of input point.
        angles (array float): values of free parameters.
        
    Ret.
        cost value for a single input.
    """
    output = circ(point, angles) # circ implements the circuit and returns the final state
    return 0.5*np.linalg.norm(vectorized_result(label)-f(output))**2 # f composes the target function from the final state

def C(data, angles):
    """Computes the cost function over the whole data set.
    
    Args.
        data (array float): set of points as pairs coordinates-label.
        angles (array float): values of free parameters.
    Ret.
        cost value for the entire set.
    """
    c = 0
    n = len(data[0])
    for i in range(n):
        c += Cp(data[0][i], data[1][i], a) # Cp defined right above
    return c / n

def cost(angles):
    """Computes the cost over the training data set. Useful for minimizing with SciPy.
    
    Args.
        angles (array float): values of free parameters.
        
    Ret.
        ret (float): cost value for the entire set.
    """
    ret = C(training_data, angles)
    print(ret)
    return ret

def vectorized_result(j):
    """Converts a digit 0,1,2,3 into a corresponding desired output [1,0,0,0],...,[0,0,0,1].
    
    Args.
        j (int): position where a 1 should be found.
    
    Ret.
        e (dim=4 int): vectorized form of j.
    """
    e = np.zeros((4, 1))
    e[j] = 1.0
    return e


def f(state):
    """Computes the target function for a given output (for 4 qubits).
    
    Args.
        state (dim=16 complex): final state of quantum system.
        
    Ret.
        target (dim=4 float): target function used for learning.
    """
    state = state * np.conj(state)
    target = np.zeros(4)
    target = [
        sum(state[i] for i in range(8,15)),
        sum(state[i]+state[i+8] for i in range(4,7)) + state[7],
        sum(state[i]+state[i+1] for i in range(2,14,4)) + state[14],
        sum(state[i] for i in range(1,15,2))
    ]
    return target
