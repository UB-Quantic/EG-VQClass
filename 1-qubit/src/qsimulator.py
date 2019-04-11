import numpy as np
import math
import cmath
import random
Pi = math.pi

"-------------------------------------------"
"GO TO THE END, ACCURACY2 UNDER CONSTRUCTION"
"-------------------------------------------"

class QC(object):

    def __init__(self, qubits, depth):
        self.size = qubits
        self.depth = depth
        """
        The quantum state is initialized with all qubits at 0.
        """
        self.state = [0]*2**self.size
        self.state[0] = 1
        """
        The angles are initialized at random with normal distribution
        """
        self.angles = [np.random.randn(3) for i in range(self.depth)]

    def initialize(self):
        """Restores the quantum state to the initial one.
        """
        self.state = [0]*2**self.size
        self.state[0] = 1

    ###############################
    # 1-Qubit Gates
    ###############################
            
        
    def h(self, m):
        """Apply the Hadamard Gate on the m'th qubit.
        Args.
            m (int): the qubit we apply our gate on.
        """
        s = 1/np.sqrt(2)
        if m>=self.size: raise ValueError('Qubit does not exist.')
        for i in range(2**(self.size-1)):
            I = 2*i-i%(2**m)
            J = I+2**m
            a = s*self.state[I] + s*self.state[J]
            b = s*self.state[I] - s*self.state[J]
            self.state[I] = a
            self.state[J] = b

    def x(self, m):
        """Apply the X Pauli Gate on the m'th qubit.
        Args.
            m (int): the qubit we apply our gate on.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        for i in range(2**(self.size-1)):
            I = 2*i-i%(2**m)
            J = I+2**m
            a = self.state[I]
            self.state[I] = self.state[J]
            self.state[J] = a

    def y(self, m):
        """Apply the Y Pauli Gate on the m'th qubit.
        Args.
            m (int): the qubit we apply our gate on.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        for i in range(2**(self.size-1)):
            I = 2*i -i%(2**m)
            J = I+2**m
            a = -1.j * self.state[I]
            self.state[I] = 1.j*self.state[J]
            self.state[J] = a

    def z(self, m):
        """Apply the Z Pauli Gate on the m'th qubit.
        Args.
            m (int): the qubit we apply our gate on.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        for i in range(2**(self.size-1)):
            J = 2*i - i%(2**m) + 2**m
            self.state[J] *= -1

    def s(self, m):
        """Apply the Phase Gate on the m'th qubit.
        Args.
            m (int): the qubit we apply our gate on.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        for i in range(2**(self.size-1)):
            J = 2*i - i%(2**m) + 2**m
            self.state[J] *= 1.j
                
    def t(self, m):
        """Apply the pi/8 Gate on the m'th qubit.
        Args.
            m (int): the qubit we apply our gate on.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        aux = cmath.exp(0.25j*math.pi)
        for i in range(2**(self.size-1)):
            J = 2*i - i%(2**m) + 2**m
            self.state[J] *= aux
               
    def rx(self, m, th):
        """Apply a x-rotation on the m'th qubit.
        Args.
            m (int): the qubit we apply our gate on.
            th (float): angle we rotate.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        th2 = 0.5*th
        c = math.cos(th2)
        s = -1.j * math.sin(th2) # beware of conventions
        for i in range(2**(self.size-1)):
            I = 2*i - i%2**m
            J = I + 2**m
            a = c*self.state[I] + s*self.state[J]
            b = s*self.state[I] + c*self.state[J]
            self.state[I] = a
            self.state[J] = b
        
    def ry(self, m, th):
        """Apply a y-rotation on the m'th qubit.
        Args.
            m (int): the qubit we apply our gate on.
            th (float): angle we rotate.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        th2 = 0.5*th
        c = math.cos(th2)
        s = math.sin(th2) # beware of conventions
        for i in range(2**(self.size-1)):
            I = 2*i - i%2**m
            J = I + 2**m
            a = c*self.state[I] - s*self.state[J]
            b = s*self.state[I] + c*self.state[J]
            self.state[I] = a
            self.state[J] = b
        
    def rz(self, m, th):
        """Apply a z-rotation on the m'th qubit.
        Args.
            m (int): the qubit we apply our gate on.
            th (float): angle we rotate.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        aux1 = cmath.exp(0.5j*th)
        aux2 = cmath.exp(-0.5j*th)
        for i in range(2**(self.size-1)):
            I = 2*i - i%2**m
            J = I + 2**m
            self.state[I] *= aux1
            self.state[J] *= aux2
        
                
    #######################################
    # 2-Qubit Gates, Entanglement
    #######################################
            
    def cnot(self, c, t):
        """Apply a Controlled-NOT gate.
        Args.
            c (int): control qubit.
            t (int): target qubit.
        """
        if c>=self.size: raise ValueError('Control does not exist.')
        if t>=self.size: raise ValueError('Target does not exist.')
        if c==t: raise ValueError('Control and Target cannot be the same.')
        for i in range(2**(self.size-2)):
            I = (2**c + i%2**c + ((i-i%2**c)*2)%2**t + 2*((i-i%2**c)*2 -
                 ((2*(i-i%2**c))%2**t)))
            J = I + 2**t
            self.state[I], self.state[J] = self.state[J], self.state[I]
                
    def cz(self, c, t):
        """Apply a Controlled-Z gate.
        Args.
            c (int): control qubit.
            t (int): target qubit.
        """
        if c>=self.size: raise ValueError('Control does not exist.')
        if t>=self.size: raise ValueError('Target does not exist.')
        if c==t: raise ValueError('Control and Target cannot be the same.')
        if t<c: t,c = c,t
        for i in range(2**(self.size-2)):
            I = (2**c + i%2**c + ((i-i%2**c)*2)%2**t + 2*((i-i%2**c)*2 -
                 ((2*(i-i%2**c))%2**t)) + 2**t)
            self.state[I] *= -1

    def swap(self, m, n):
        """Apply a SWAP gate.
        Args.
            m (int): first qubit.
            n (int): second qubit.
        """
        if m>=self.size: raise ValueError('First Qubit does not exist.')
        if n>=self.size: raise ValueError('Second Qubit does not exist.')
        if m==n: raise ValueError('Both Qubits cannot be the same.')
        for i in range(2**(self.size-2)):
            I = (i%2**m + ((i-i%2**m)*2)%2**n + 2*((i-i%2**m)*2 - 
                 ((2*(i-i%2**m))%2**n)) + 2**n)
            J = I + 2**m - 2**n
            self.state[I], self.state[J] = self.state[J], self.state[I]
    ############################################
    # Circuits
    ############################################
   
    # The following are intended to be used with 1-qubit circuits.
    def unitary(self, m, theta, phi, lamb):
        """Apply an arbitrary unitary gate on the m'th qubit.
        Every unitary gate is characterized by three angles.
        Args.
            m (int): qubit the gate is applied on.
            theta (float): first angle.
            phi (float): second angle.
            lamb (float): third angle.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        c = math.cos(0.5*theta)
        s = math.sin(0.5*theta)
        ephi = cmath.exp(1j*phi)
        elamb = cmath.exp(1j*lamb)
        for i in range(2**(self.size-1)):
            I = 2*i -i%(2**m)
            J = I+2**m
            a = c*self.state[I] - s*elamb*self.state[J]
            b = s*ephi*self.state[I] + c*ephi*elamb*self.state[J]
            self.state[I] = a
            self.state[J] = b

    def transunit(self, m, theta, phi, lamb, vector):
        """Transpose of the unitary gate right above.
        Args.
            m (int): qubit the gate is applied on.
            theta (float): first angle.
            phi (float): second angle.
            lamb (float): third angle.
            vector (array complex): what the operator acts on.
        Ret.
            vector (array complex): image of the operation.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        c = math.cos(0.5*theta)
        s = math.sin(0.5*theta)
        ephi = cmath.exp(1j*phi)
        elamb = cmath.exp(1j*lamb)
        for i in range(2**(self.size-1)):
            I = 2*i -i%(2**m)
            J = I+2**m
            a = c*self.state[I] + s*ephi*self.state[J]
            b = -s*elamb*self.state[I] + c*ephi*elamb*self.state[J]
            self.state[I], self.state[J] = a, b
        return vector

    def difunit1(self, m, theta, phi, lamb, vector):
        """Unitary gate differentiated with respect to theta on m'th qubit.
        Not unitary anymore.
        Args.
            m (int): qubit the operator is applied on.
            theta (float): first angle.
            phi (float): second angle.
            lamb (float): third angle.
            vector (array complex): what the operator acts on.
        Ret.
            vector (array complex): image of the operation.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        c = 0.5*math.cos(0.5*theta)
        s = 0.5*math.sin(0.5*theta)
        ephi = cmath.exp(1j*phi)
        elamb = cmath.exp(1j*lamb)
        for i in range(2**(self.size-1)):
            I = 2*i-i%(2**m)
            J = I+2**m
            a = -s*vector[I] -c*elamb*vector[J]
            b = c*ephi*vector[I] -s*ephi*elamb*vector[J]
            vector[I], vector[J] = a, b
        return vector

    def difunit2(self, m, theta, phi, lamb, vector):
        """Unitary operator differentiated with respect to phi on m'th qubit.
        Not unitary anymore.
        Args.
            m (int): qubit the operator is applied on.
            theta (float): first angle.
            phi (float): second angle.
            lamb (float): third angle.
            vector (array complex): what the operator acts on.
        Ret.
            vector (array complex): image of the operation.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        c = math.cos(0.5*theta)
        s = math.sin(0.5*theta)
        ephi = 1j*cmath.exp(1j*phi)
        elamb = cmath.exp(1j*lamb)
        for i in range(2**(self.size-1)):
            I = 2*i-i%(2**m)
            J = I+2**m
            vector[J] = s*ephi*vector[I]+c*ephi*elamb*vector[J]
            vector[I] = 0
        return vector

    def difunit3(self, m, theta, phi, lamb, vector):
        """Unitary operator differentiated with respect to lamb on m'th qubit.
        Not unitary anymore.
        Args.
            m (int): qubit the operator is applied on.
            theta (float): first angle.
            phi (float): second angle.
            lamb (float): third angle.
            vector (array complex): what the operator acts on.
        Ret.
            vector (array complex): image of the operation.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        c = math.cos(0.5*theta)
        s = math.sin(0.5*theta)
        ephi = cmath.exp(1j*phi)
        elamb = 1j*cmath.exp(1j*lamb)
        for i in range(2**(self.size-1)):
            I = 2*i-i%(2**m)
            J = I+2**m
            vector[I] = -s*elamb*vector[J]
            vector[J] = c*ephi*elamb*vector[J]
        return vector


    def block(self, m, point, angles, style=1):
        """Apply a learning block on the m'th qubit.
        Args.
            m (int): qubit the block is applied on.
            point (dim=2 float): coordinates of input.
            angles (dim=3 float): angles that determine a unitary gate.
            style (int): customizes the block.
        """
        if m>=self.size: raise ValueError('Qubit does not exist.')
        if style:
            self.unitary(m, point[0]+angles[0], point[1]+angles[1], angles[2])
        else:
            self.ry(m, point[0]*0.5*Pi)
            self.rz(m, (1+point[1])*Pi)
            self.unitary(m, angles[0], angles[1], angles[2])

    def run(self, point, parameters, style=1):
        """Runs the circuit and restores to the initial state.
        Args.
            point (dim=2 float): coordinates of input.
            style (int): customizes the block.
        Ret.
            p0 (float): probability of the final state being |0>.
        """
        for layer in parameters:
            self.block(0, point, layer, style)
        p0 = self.state[0]
        p0 = p0*np.conj(p0)
        self.initialize()
        return p0

    def Cp(self, point, label, parameters):
        """Computes the cost value of a single input.
        Args.
            point(dim=2 float): coordinates of input point.
            label(int): what class does point belong to.
            parameters(array float): NEEDS BE GOTTEN RID OF.
        Ret.
            cp (float): cost value for a single input.
        """
        cp = np.linalg.norm(label - 3*self.run(point, parameters))
        cp = 0.5*cp*cp
        return cp


# The day shall come, when we don't need parameters to be an argument
# of the cost. That is the day when we will know how to analytically
# compute the gradient. Until that day comes, MasterCard.
    def C(self, data, parameters):
        """Computes the cost function over the whole data set.
        Args.
            data (array float): set of input points.
            parameters (array float): NEEDS BE GOTTEN RID OF.
        Ret.
            c (float): cost value for the entire set.
        """
        c = 0
        n = len(data[0])
        for i in range(n):
            c+=self.Cp(data[0][i], data[1][i], parameters)
        c/=n
        return c

    def gradC(self, data, parameters, step):
        """Computes the gradient of the cost function.
        Args.
            data (array float): set of input points.
            parameters (array float): NEEDS BE GOTTEN RID OF.
            step (float): nummerical differentiation step.
        Ret.
            nabla (array float): gradient vector.
        """
        nabla = ([np.zeros(a.shape) for a in parameters])
        dif = np.asarray([np.zeros(a.shape) for a in parameters])
        half=0.5*step
        inv=1/step
        for i in range(len(dif)):
            for j in range(len(dif[i])):
                dif[i][j] = half
                nabla[i][j] = (self.C(data, parameters+dif)-
                               self.C(data, parameters-dif))*inv
                dif[i][j] = 0
        return nabla

    #########################################
    # Quantum Backpropagation
    #########################################

    def SGD(self, training_data, epochs, mini_batch_size, learning_rate,
            test_data):
        """Train a variational circuit using Stochastic Gradient Descent.
        Args.
            training_data (array float): set of training input points.
            epochs (int): number of learning epochs performed.
            mini_batch_size (int): size of the learning batches.
            learning_rate (float): distance moved on every step.
            test_data (array float): set of test input points.
        """
        n = len(training_data)
        if test_data:
            test = []
            n_test = len(test_data)
            for i in range(n_test):
                test.append([test_data[0][i],test_data[1][i]])
        for j in range(epochs):
            comb = list(zip(training_data[0],training_data[1]))
            random.shuffle(comb)
            #training_data[0][:], training_data[1][:]=zip(*comb)
            mini_batches = [
                training_data[k:k+mini_batch_size] for k in
                range(0,n,mini_batch_size)
            ]
            for mini_batch in mini_batches:
                self.update_mini_batch(mini_batch, learning_rate)
            print('after epoch',j,'score is: ', self.accuracy(test_data))
            print('\tcost = ',self.C(test_data,self.angles))

    def update_mini_batch(self, mini_batch, learning_rate):
        """Propose a new set of parameters for an input batch.
        Args.
            mini_batch (array float): set of input points.
            learning_rate (float): distance moved along the gradient.
        """
        new_angles = [np.zeros(a.shape) for a in self.angles]

        for x,y in zip(mini_batch[0], mini_batch[1]):
            delta_new_angles = self.backpropagate(x,y)
            new_angles = [na + dna for na, dna in zip(new_angles,
                                                      delta_new_angles)]
        self.angles = [a - (learning_rate/len(mini_batch)) * na 
                       for a, na in zip(self.angles, new_angles)]

    def backpropagate(self, x, y):
        """Propose a new set of angles using a backpropagation algorithm.
        Args.
            x (dim=2 float): coordinates of input point.
            y (int): desired output for input point.
        Ret.
            new_angles (array float): proposal of new angles for one input.
        """
        # CAREFUL! the gradient is not perfectly calculated
        # and so the derivatives dC/dth are complex, 
        # which does not make much sense. This leads to poor performance
        # since the parameters update is kind of unkown to us right now.
        new_angles = [np.zeros(a.shape) for a in self.angles]
        #---------------FeedForward------------
        activations = [self.state.copy()]
        #activations.append(list(self.state))
        for a in self.angles:
            self.unitary(0,a[0],a[1],a[2])
            activations.append(list(self.state))
        #---------------Backward pass----------
            # First step
        delta = [-3*np.conj(activations[-1][0])*
                 (y-3*activations[-1][0]*np.conj(activations[-1][0])),0]
        dif_act = self.difunit1(0, self.angles[-1][0]+x[0],
                                self.angles[-1][1]+x[1],
                                self.angles[-1][2], activations[-2].copy())
        new_angles[-1][0] = np.dot(delta,dif_act)
        dif_act = self.difunit2(0, self.angles[-1][0]+x[0],
                                self.angles[-1][1]+x[1],
                                self.angles[-1][2], activations[-2].copy())
        new_angles[-1][1] = np.dot(delta,dif_act)
        dif_act = self.difunit3(0, self.angles[-1][0]+x[0],
                                self.angles[-1][1]+x[1],
                                self.angles[-1][2], activations[-2].copy())
        new_angles[-1][2] = np.dot(delta,dif_act)
            # Recursive steps
        for l in range(1, self.depth):
            psi = activations[-l-2].copy()
            delta = self.transunit(0, self.angles[-l][0]+x[0],
                                   self.angles[-l][1]+x[1],
                                   self.angles[-l][2], delta)
            dif_act = self.difunit1(0, self.angles[-l-1][0]+x[0],
                                    self.angles[-l-1][1]+x[1],
                                    self.angles[-l-1][2],psi)
            new_angles[-l-1][0] = np.dot(delta,dif_act)
            dif_act = self.difunit2(0, self.angles[-l-1][0]+x[0],
                                    self.angles[-l-1][1]+x[1],
                                    self.angles[-l-1][2],psi)
            new_angles[-l-1][1] = np.dot(delta,dif_act)
            dif_act = self.difunit3(0, self.angles[-l-1][0]+x[0],
                                    self.angles[-l-1][1]+x[1],
                                    self.angles[-l-1][2],psi)
            new_angles[-l-1][2] = np.dot(delta,dif_act)
        self.initialize()
        return new_angles

    def accuracy(self, data):
        """Computes the classification accuracy.
        Args.
            data (array float): input data set.
        Ret.
            score (int): number of inputs correctly classified.
        """
        results = []
        for (x,y) in zip(data[0],data[1]):
            p0 = self.run(x,self.angles)
            if p0<0.25: p=0
            elif p0<0.5: p=1
            elif p0<0.75: p=2
            else: p=3
            results.append([p,y])
        score=sum(int(x==y) for (x,y) in results)
        return score

    ######################################
    # Quantum Backpropagation v2
    ######################################
    def run2(self, point, style=1):
        """Runs the circuit and restores to the initial state.
        Args.
            point (dim=2 float): coordinates of input.
            style (int): customizes the block.
        Ret.
            final_state (dim=2 float): wavefunction at the end.
        """
        for layer in self.angles:
            self.block(0, point, layer, style)
        final_state = self.state.copy()
        self.initialize()
        return final_state

    def fp(self, point, label):
        """Computes the fidelity of a single input.
        Args.
            point(dim=2 float): coordinates of input point.
            label(int): what class does point belong to.
        Ret.
            fp (float): fidelity for a single input.
        """
        fp = np.dot(self.run2(point),label)
        fp = (np.conj(fp)*fp).real
        return fp

    def C2(self, data):
        """Computes the cost function over thw whole data set.
        Args.
            data (array float): set of input points.
        Ret.
            cost (float): cost value for the entire set.
        """
        cost = 0
        n = len(data[0])
        for point, label in zip(data[0],data[1]):
            cost+=self.fp(point,label)
        cost/=n
        return cost

    def JISGD(self, training_data, epochs, mini_batch_size, learning_rate,
               test_data):
        """Train a variational circuit using Stochastic Gradient Descent.
        Args.
            training_data (array float): set of training input points.
            epochs (int): number of learning epochs performed.
            mini_batch_size (int): size of the learning batches.
            learning_rate (float): distance moved on every step.
            test_data (array float): set of test input points.
        """
        n = len(training_data[0])
        sq3 = 1/math.sqrt(3)
        sq23 = math.sqrt(2/3)
        exp23 = cmath.exp(2j*Pi/3)
        if test_data:
            new_test = []
            test = []
            n_test = len(test_data[0])
            for i in range(n_test):
                if test_data[1][i] == 0:
                    new_test.append([1,0])
                elif test_data[1][i] == 1:
                     new_test.append([sq3,sq23])
                elif test_data[1][i] == 2:
                   new_test.append([sq3,sq23*exp23])
                elif test_data[1][i] == 3:
                    new_test.append([sq3,sq23/exp23])
                test.append([test_data[0][i],new_test[i]])

        # Debugging needed
        new_labels = [] 
        for i in range(n):
            if training_data[1][i] == 0:
                new_labels.append([1,0])
            elif training_data[1][i] == 1:
                new_labels.append([sq3,sq23])
            elif training_data[1][i] == 2:
                new_labels.append([sq3,sq23*exp23])
            elif training_data[1][i] == 3:
                new_labels.append([sq3,sq23/exp23])
        
        for j in range(epochs):
            comb = list(zip(training_data[0],new_labels))
            random.shuffle(comb)
            
            mini_batches = [
                comb[k:k+mini_batch_size] for k in
                range(0,n,mini_batch_size)
            ]
            for mini_batch in mini_batches:
                self.update_mini_batch2(mini_batch, learning_rate)
#            print('after epoch',j,'score is: ', self.accuracy(test_data))
            print('\tcost = ',self.C2(test))
            
    def update_mini_batch2(self, mini_batch, learning_rate):
        """Propose a new set of parameters for an input batch.
        Args.
            mini_batch (array float): set of input points.
            learning_rate (float): distance moved along the gradient.
        """
        new_angles = [np.zeros(a.shape) for a in self.angles]

        for x,y in zip(mini_batch[0], mini_batch[1]):
            delta_new_angles = self.backpropagate2(x,y)
            new_angles = [na + dna for na, dna in zip(new_angles,
                                                      delta_new_angles)]
        self.angles = [a + (learning_rate/len(mini_batch)) * na 
                       for a, na in zip(self.angles, new_angles)]

    def backpropagate2(self, x, y):
        """Propose a new set of angles using a backpropagation algorithm.
        Args.
            x (dim=2 float): coordinates of input point.
            y (dim=2 complex): desired output state for input point.
        Ret.
            new_angles (array float): proposal of new angles for one input.
        """
        new_angles = [np.zeros(a.shape) for a in self.angles]
        #-----------------FeedForward----------------
        activations = [self.state.copy()]
        for a in self.angles:
            self.unitary(0,a[0],a[1],a[2])
            activations.append(self.state.copy())
        fid = np.dot(np.conj(activations[-1]),y)
        self.initialize()
        #-----------------Backward pass--------------
        # First step

        # ENHANCEMENT:
        # We should be coding this taking advantadge of the vector 
        # nature of angles, activations, and all.
        delta = y
        dif_act = self.difunit1(0, self.angles[-1][0]+x[0],
                               self.angles[-1][1]+x[1],
                               self.angles[-1][2], activations[-2].copy())
        new_angles[-1][0] = 2*(fid*np.dot(np.conj(delta),dif_act)).real
        dif_act = self.difunit2(0, self.angles[-1][0]+x[0],
                               self.angles[-1][1]+x[1],
                               self.angles[-1][2], activations[-2].copy())
        new_angles[-1][1] = 2*(fid*np.dot(np.conj(delta),dif_act)).real
        dif_act = self.difunit3(0, self.angles[-1][0]+x[0],
                               self.angles[-1][1]+x[1],
                               self.angles[-1][2], activations[-2].copy())
        new_angles[-1][2] = 2*(fid*np.dot(np.conj(delta),dif_act)).real
        # Recursive steps
        for l in range(1, self.depth):
            delta = self.transunit(0, self.angles[-l][0]+x[0],
                                   -self.angles[-l][1]-x[1],
                                   -self.angles[-l][2], delta)
            psi = activations[-l-2].copy()
            dif_act = self.difunit1(0, self.angles[-l-1][0]+x[0],
                                    self.angles[-l-1][1]+x[1],
                                    self.angles[-l-1][2], psi)
            new_angles[-l-1][0] = 2*(fid*np.dot(np.conj(delta),dif_act)).real
            dif_act = self.difunit2(0, self.angles[-l-1][0]+x[0],
                                    self.angles[-l-1][1]+x[1],
                                    self.angles[-l-1][2], psi)
            new_angles[-l-1][1] = 2*(fid*np.dot(np.conj(delta),dif_act)).real
            dif_act = self.difunit3(0, self.angles[-l-1][0]+x[0],
                                    self.angles[-l-1][1]+x[1],
                                    self.angles[-l-1][2], psi)
            new_angles[-l-1][2] = 2*(fid*np.dot(np.conj(delta),dif_act)).real


        return new_angles

# Choose a way of classifying with the states!
    def accuracy2(self, data):
        """Computes the classification accuracy.
        Args.
            data (array float): input data set.
        Ret.
            score (int): number of inputs correctly classified.
        """
        results = []
        for (x,y) in zip(data[0], data[1]):
            psi = self.run2(x):

import datagen

training_data, test_data = datagen.read('../data/data3.txt', 50, 50)
qlass = QC(1,1)
qlass.JISGD(training_data, 30, 5, 0.5, test_data)
#print('backprop 1: ', qlass.backpropagate([0.05,0.95],0))
#print('backprop 2: ', qlass.backpropagate2([0.05,0.95],np.asarray([1,0])))
#print('gradC : ', qlass.gradC([[[0.05,0.95]],[0]],qlass.angles,0.01))
