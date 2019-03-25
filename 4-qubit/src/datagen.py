import numpy as np

def writecircles (file, n):
    """Creates a random data set, suitable for our task.

    Args:
            file (str): file name where we write.
            n (int): number of points to be generated.
    """

    circle_x = [
        0.,
        -1.,
        0.4,
        0.4
    ]

    circle_y = [
        0.,
        0.,
        0.4,
        -0.6
    ]

    circle_r = [
        2.,
        0.8,
        0.4,
        0.6
    ]

    x, y, a = [], [], []
    for i in range(n):
        x.append(np.random.uniform(-1, 1))
        y.append(np.random.uniform(-1, 1))
        for (cx, cy, cr, i) in zip(circle_x, circle_y, circle_r, range(len(circle_x))):
            if((x[-1] - cx) ** 2 + (y[-1] - cy) ** 2 <= cr ** 2):
                rest = i
        a.append(rest)

    with open(file,'w') as f:
        for i in range(n):
            f.write("%+1.8f"% x[i] + "\t" + "%+1.8f"%y[i] + "\t" + "%d"%a[i] + "\n")

def read(file, n, m):
    """Takes and formats the data from a given file.

    Args.
        file (str): name of the file the data should be read from.
        n (int): size of the training set.
        m (int): size of the test set.

    Ret.
        training_data (array float): set of points as pairs coordinates-label used for training.
        test_data (array float): set of points as pairs coordinates-label used for testing.
    """
    x = np.ndarray(shape = (n + m, 2), dtype = float)
    a = np.ndarray(n + m, dtype = int)


    f=open(file,"r")
    lines=f.readlines()
    i = 0
    for line in lines:
        b = float(line.split('\t')[0])
        c = float(line.split('\t')[1])
        x[i] = [b, c]
        a[i] = line.split('\t')[2]
        i += 1
    f.close()

    training_set = x[:n]
    training_output = a[:n]
    test_set = x[n : n+m]
    test_output = a[n : n+m]

    training_data = (training_set, training_output)
    test_data = (test_set, test_output)

    return (training_data, test_data)
