import numpy as np

def write(file, n):
    """Creates a random data set.

    Args.
        file(str): file name where data is written.
        n (int): number of generated points.
    """

    circle_x = [
        0,
        -0.5,
        -1,
        1
    ]

    circle_y = [
        0,
        -0.5,
        1,
        0
    ]

    circle_r = [
        2,
        0.5,
        1,
        1
    ]

    x, y, a = [],[],[]
    for i in range(n):
        x.append(np.random.uniform(-1,1))
        y.append(np.random.uniform(-1,1))
        for (cx, cy, cr, i) in zip(circle_x, circle_y, circle_r,
                                   range(len(circle_x))):
            if((x[-1]-cx)*(x[-1]-cx) + (y[-1]-cy)*(y[-1]-cy)<=(cr*cr)):
                rest = i
        a.append(rest)

    with open(file,'w') as f:
        for i in range(n):
            f.write(
                "%+1.8f"%x[i]+"\t"+"%+1.8f"%y[i]+"\t"+"%d"%a[i]+"\n")

def read(file, n, m):
    """Takes and formats the data from a given file.

    Args.
        file (str): name of the file the data should be read from.
        n (int): size of the training set.
        m (int): size of the test set.

    Ret.
        training_data (array float): data set used for training.
        test_data (array float): data set used for testing.
    """
    x = np.ndarray(shape = (n + m, 2), dtype = float)
    a = np.ndarray(n + m, dtype = int)

    f=open(file,"r")
    lines=f.readlines()
    i = 0
    for line in lines:
        aux = line.split('\t')
        x[i] = [float(aux[0]), float(aux[1])]
        a[i] = aux[2]
        i += 1
    f.close()

    training_set = x[:n]
    training_output = a[:n]
    test_set = x[n : n+m]
    test_output = a[n : n+m]

    training_data = (training_set, training_output)
    test_data = (test_set, test_output)

    return (training_data, test_data)
