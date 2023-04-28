from enum import Enum 

class PauliGate(Enum):
    I=0
    X=1
    Y=2
    Z=3

I=PauliGate.I
X=PauliGate.X
Y=PauliGate.Y
Z=PauliGate.Z

class PauliString:
    def __init__(self, coef, gates):
        self.coef=coef
        self.gates=gates

    def __eq__(self, other):
        return self.gates==other.gates
    
# a0 III + a1 IIX + a2 XYZ + ...
class PauliStringExpr:
    def __init__(self, data):
        self.expr=data

    def __add__(self, expr):
        if expr.gates in self.expr:
            self.expr[expr.gates]+=expr.coef
            if self.expr[expr.gates]==0:
                self.expr.pop(expr.gates)
        else:
            self.expr[expr.gates]=expr.coef