from enum import Enum 

class PauliGate(Enum):
    I=0
    X=1
    Y=2
    Z=3

    def __mul__(self, other):
        if isinstance(other, PauliGate):
            return PauliString(1.0,str(self)+str(other))
        elif isinstance(other, float):
            return PauliString(other, str(self))
        raise TypeError("Cannot multiply PauliGate with {}".format(type(other)))

    __rmul__=__mul__

    def __str__(self):
        return str(self.name)

I=PauliGate.I
X=PauliGate.X
Y=PauliGate.Y
Z=PauliGate.Z

class PauliString:
    def __init__(self, coef, gates):
        self.coef=coef
        self.gates=gates

    def __mul__(self, other):
        if isinstance(other, PauliString):
            print(self, other)
            print(type(self), "  ", type(other))
            self.coef*=other.coef
            self.gates+=other.gates
            return self
        
        elif isinstance(other, PauliGate):
            self.gates+=str(other)
            return self
        
        raise TypeError("Cannot multiply PauliString with {}".format(type(other)))

    def __eq__(self, other):
        return self.gates==other.gates
    
    def __str__(self):
        return str(self.coef)+"*"+self.gates

# a0 III + a1 IIX + a2 XYZ + ...
class PauliStringExpr:
    def __init__(self, data):
        pass