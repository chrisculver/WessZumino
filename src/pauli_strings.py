from enum import Enum 

# This could be a nice standalone package, 
# TODOS:
#
# 1.  Users shouldn't specify PS like X*X + Y*X*I, arbitrary I padding on first expression needed
# 2.  Should the reduction of 0 valued PS be done always?  Is it faster?
#     With string overhead, maybe
#     With ints.... maybe not... 
# 3.  Storing main data as a string...  probably better way to sort?
#     PS as base 4 numbers? i.e. IIII=0000=0  ZZZZ=3333=63?
#     then they can be sorted as ints (FAST!), and strings can easily be 
#     obtained from a dictionary of max size 4^NQ?
#
#     DEF the above should be faster and straightforward to implement
#
#     At that point why use a dictionary and not an array of size 4^NQ, 
#     and never worry about reducing coefs, default is zero...?  I guess
#     the exponential scaling is bad if you think the final result will 
#     be sparse in the types of pauli gates, I don't have an idea at 
#     the moment if random PSE's are like this
#     


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
    
    def __add__(self,other):
        return PauliStringExpr({self.gates: self.coef, other.gates: other.coef})

    def __str__(self):
        return str(self.coef)+"*"+self.gates

# a0 III + a1 IIX + a2 XYZ + ...
class PauliStringExpr:
    def __init__(self, data):
        self.data=data
    
    def __add__(self,other):
        #print("Adding {} and {}".format(self,other))
        if isinstance(other, PauliString):
            if other.gates in self.data:
                self.data[other.gates]+=other.coef 
                
                if self.data[other.gates]==0:
                    self.data.pop(other.gates)
        
            else:
                self.data[other.gates]=other.coef

            return self

        elif isinstance(other, PauliStringExpr):
            for gates,coef in other.data.items():
                if gates in self.data:
                    self.data[gates]+=coef
                    if self.data[gates]==0:
                        self.data.pop(gates)
                else:
                    self.data[gates]=coef

            return self
    
        else: 
            raise TypeError("Cannot add PauliStringExpr with {}".format(type(other)))

    def __sub__(self,other):
        #print("Subbing {} and {}".format(self,other))
        if isinstance(other, PauliString):
            other.coef*=-1
            return self+other
        
        elif isinstance(other, PauliStringExpr):
            for gates in other.data:
                other.data[gates]*=-1
            return self+other 
        
        else:
            raise TypeError("Cannot sub PauliStringExpr with {}".format(type(other)))
        
    def __mul__(self, other):
        if isinstance(other, (int, float, complex)):
            for gates in self.data:
                self.data[gates]*=other 
            return self 
    
        elif isinstance(other, PauliStringExpr):
            for gates, coef in other.items():
                pass
            return self

        else: 
            raise TypeError("Cannot mul PauliStringExpr with {}".format(type(other)))

    # TODO: Is this correct?
    __rmul__=__mul__

    def __str__(self):
        s=""
        for gates,coef in self.data.items():
            s+=str(coef)+"*"+str(gates)
            s+="+"
        return s[:-1]
    

