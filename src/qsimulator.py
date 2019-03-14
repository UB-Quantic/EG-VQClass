from itertools import product as prod
import numpy as np
import math
import cmath
import matplotlib.pyplot as plt

Pi = math.pi

def ii(i):
    """Computes the double index.
    
    Args.
        i (vector): array of bits with value 1 or 0.

    Ret.
        ii (int): base 10 representation of the bits array.
    """
    l = len(i)
    return np.dot(i,[2**(l-j-1) for j in range(l)])

 


class QC(object):

    def __init__(self, qubits):
        self.size = qubits
        """
        The quantum state is initialized with all qubits at 0.
        """
        self.state = [0]*2**self.size
        self.state[0] = 1.

    def initialize(self):
        """Brings the state vector back to its initial state.
        """
        self.state = [0]*2**self.size
        self.state[0] = 1.


    
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
        for i in prod([0,1], repeat=self.size):
            j = np.array(i)
            j[m] ^= 1
            if(i[m]):
                self.state[self.ii(i)], self.state[self.ii(j)] = (
                    self.state[self.ii(j)],
                    -self.state[self.ii(i)]
                )
        self.state *= 1.j
        
    def z(self, m):
        """Apply the Z Pauli Gate on the m'th qubit.
        
        Args.
            m (int): the qubit we apply our gate on.
        """
        for i in prod([0,1], repeat=self.size):
            if(i[m]):
                self.state[self.ii(i)] *= -1
    
    def s(self, m):
        """Apply the Phase Gate on the m'th qubit.
        
        Args.
            m (int): the qubit we apply our gate on.
        """
        for i in prod([0,1], repeat=self.size):
            if(i[m]):
                self.state[self.ii(i)] *= 1.j
                
    def t(self, m):
        """Apply the pi/8 Gate on the m'th qubit.
        
        Args.
            m (int): the qubit we apply our gate on.
        """
        aux = cmath.exp(0.25j*math.pi)
        for i in prod([0,1], repeat=self.size):
            if(i[m]):
                self.state[self.ii(i)] *= aux
                
    def rx(self, m, th):
        """Apply a x-rotation on the m'th qubit.
        
        Args.
            m (int): the qubit we apply our gate on.
            th (float): angle we rotate.
        """
        c = math.cos(th / 2)
        s = -1.j * math.sin(th / 2) # beware of conventions
        for i in prod([0,1], repeat=self.size):
            j = np.array(i)
            j[m] ^= 1
            if(i[m]):
                self.state[self.ii(i)],self.state[self.ii(j)] = (
                    s*self.state[self.ii(j)] + c*self.state[self.ii(i)],
                    c*self.state[self.ii(j)] + s*self.state[self.ii(i)]
                )
    
    def ry(self, m, th):
        """Apply a y-rotation on the m'th qubit.
        
        Args.
            m (int): the qubit we apply our gate on.
            th (float): angle we rotate.
        """
        c = math.cos(th / 2)
        s = math.sin(th / 2) 
        for i in prod([0,1], repeat=self.size):
            j = np.array(i)
            j[m] ^= 1
            if(i[m]):
                self.state[self.ii(i)],self.state[self.ii(j)] = (
                    s*self.state[self.ii(j)] + c*self.state[self.ii(i)],
                    c*self.state[self.ii(j)] - s*self.state[self.ii(i)]
                )
                # beware of conventions
                
    def rz(self, m, th):
        """Apply a z-rotation on the m'th qubit.
        
        Args.
            m (int): the qubit we apply our gate on.
            th (float): angle we rotate.
        """
        aux = cmath.exp(0.5j*th)
        for i in prod([0,1], repeat=self.size):
            if(i[m]):
                self.state[self.ii(i)] *= aux
            if(not i[m]):
                self.state[self.ii(i)] /= aux
                
        
                
    #######################################
    # 2-Qubit Gates, Entanglement
    #######################################
            
    def cnot(self, c, t):
        """Apply a Controlled-NOT gate.
        
        Args.
            c (int): control qubit.
            t (int): target qubit.
        """
        if(c==t):
            print('Error: control cannot be target')
            return 
        for i in prod([0,1], repeat=self.size):
            j = np.array(i)
            j[t] ^= 1
            if(i[c] and i[t]):
                self.state[self.ii(i)],self.state[self.ii(j)]=(
                    self.state[self.ii(j)],
                    self.state[self.ii(i)]
                )
                
    def cz(self, c, t):
        """Apply a Controlled-Z gate.
        
        Args.
            c (int): control qubit.
            t (int): target qubit.
        """
        if(c==t):
            print('Error: control cannot be target')
            return 
        for i in prod([0,1], repeat=self.size):
            if(i[c] and i[t]):
                self.state[self.ii(i)] *= -1
                
    ############################################
    # Circuits
    ############################################
    
    def encode(self, point):
        """Creates the encoding layer.
        
        Args.
            point (dim=2 float): coordinates of one input point.
        """
        for i in range(self.size):
            self.h(i)
            self.rz(i, point[i%2])
                
    def blocka(self, angles, qubits=[0,1,2,3]):
        """Adds a block of type a.
        
        Args.
            angles (dim=8 float): rotation angles for each gate .
            qubits (dim=4 int): qubits the block acts on.
        """
        for i in range(4):
            self.rx(qubits[i], angles[i])
        self.cz(qubits[0], qubits[1])
        self.cz(qubits[2], qubits[3])
        for i in range(4):
            self.ry(qubits[i], angles[4+i])
        self.cz(qubits[1], qubits[2])
        self.cz(qubits[0], qubits[3])
        
    def blockb(self, angles, qubits=[0,1,2,3]):
        """Adds a block of type b.
        
        Args.
            angles (dim=8 float): rotation angles for each gate.
            qubits (dim=4 int): qubits the block acts on.
        """
        for i in range(4):
            self.ry(qubits[i], angles[i]) 
        self.cz(qubits[0], qubits[1])
        self.cz(qubits[2], qubits[3])
        for i in range(4):
            self.rx(qubits[i], angles[4+i])
        self.cz(qubits[1], qubits[2])
        self.cz(qubits[0], qubits[3])
        
    def blockc(self, angles, qubits=[0,1,2,3]):
        """Adds a block of type c.
        
        Args.
        angles (dim=8 float): rotation angles for each gate.
        qubits (dim=4 int): qubits the block acts on.
        """
        self.rx(qubits[0], angles[0])
        self.ry(qubits[1], angles[1])
        self.rx(qubits[2], angles[2])
        self.ry(qubits[3], angles[3])
        self.cz(qubits[0], qubits[1])
        self.cz(qubits[2], qubits[3])
        self.ry(qubits[0], angles[4])
        self.rx(qubits[1], angles[5])
        self.ry(qubits[2], angles[6])
        self.rx(qubits[3], angles[7])
        self.cz(qubits[1], qubits[2])
        self.cz(qubits[0], qubits[3])
        
    def blockd(self, angles, qubits=[0,1,2,3]):
        """Adds a block of type d.
        
        Args.
            angles (dim=8 float): rotation angles for each gate.
            qubits (dim=4 int): qubits the block acts on.
        """
        self.rx(qubits[0], angles[0])
        self.ry(qubits[1], angles[1])
        self.rx(qubits[2], angles[2])
        self.ry(qubits[3], angles[3])
        self.cz(qubits[0], qubits[1])
        self.cz(qubits[2], qubits[3])
        self.rx(qubits[0], angles[4])
        self.ry(qubits[1], angles[5])
        self.rx(qubits[2], angles[6])
        self.ry(qubits[3], angles[7])
        self.cz(qubits[1], qubits[2])
        self.cz(qubits[0], qubits[3])
        
    def blockx(self, angles, qubits=[0,1,2,3]):
        """Adds a block of type x.
        
        Args.
            angles (dim=8 float): rotation angles for each gate.
            qubits (dim=4 int): qubits the block acts on.
        """
        for i in range(4):
            self.rx(qubits[i], angles[i])
        self.cz(qubits[0], qubits[1])
        self.cz(qubits[2], qubits[3])
        for i in range(4):
            self.rx(qubits[i], angles[4+i])
        self.cz(qubits[1], qubits[2])
        self.cz(qubits[0], qubits[3])
        
    def blocky(self, angles, qubits=[0,1,2,3]):
        """Adds a block of type y.
        
        Args.
            angles(dim=8 float): rotation angles for each gate.
            qubits (dim=4 int): qubits the block acts on.
        """
        for i in range(4):
            self.ry(qubits[i], angles[i])
        self.cz(qubits[0], qubits[1])
        self.cz(qubits[2], qubits[3])
        for i in range(4):
            self.ry(qubits[i], angles[4+i])
        self.cz(qubits[1], qubits[2])
        self.cz(qubits[0], qubits[3])
        
    def add(self, typ, angles, qubits=[0,1,2,3]):
        """Adds a block of a certain type in a given position.
        
        Args.
            typ (char): type of circuit 'a', 'b', 'c' or 'd'.
            angles (dim=8 float): rotation angles for each gate.
            qubits (dim=4 int): which qubits the block acts on.
            
        Rets.
            success (int): indicates whether some error flag was raised.
        """
       
        if(typ not in 'abcdxy'):
            print("Wrong key for type.")
            return 1
        
        return {
            'a': self.blocka(angles, qubits),
            'b': self.blockb(angles, qubits),
            'c': self.blockc(angles, qubits),
            'd': self.blockd(angles, qubits),
            'x': self.blockx(angles, qubits),
            'y': self.blocky(angles, qubits)
        }.get(typ, 1)


import time

hola = QC(4)
inici1 = time.time()
for i in range(4):
    hola.x(i)
print(hola.state)
final1 = time.time()
print("temps 1 = ", final1-inici1)
hola.initialize()
inici2 = time.time()
for i in range(4):
    hola.x2(i)
print(hola.state)
final2 = time.time()
print("temps 2 = ", final2-inici2)
