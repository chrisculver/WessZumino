
class PauliGate(Enum):
    I=0
    X=1
    Y=2
    Z=3

class PauliString:
    def __init__(self, gates):
        self.gates=gates

    def __eq__(self, other):
        return self.gates==other.gates
    
class PauliStringExpr:
    def __init__(self, data):
        self.expr=data

    def __add__(self, expr):
        pass