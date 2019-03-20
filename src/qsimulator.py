import numpy as np
import math
import cmath

Pi = math.pi


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
