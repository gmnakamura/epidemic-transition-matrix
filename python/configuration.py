

import gmpy2

class Configuration(object):
    """
    represents an unique network configuration.
    * Important: in this representation, states are limited 
                 up to maximum of 10 states.
    """
    dimension=1
    base     =2

    def __init__(self,val=''):
        """ initialize the class. val may assume either string or
        integer representation. Integer decomposition is performed
        by gmpy2.
        """
        if isinstance(val,str):
            self.label = val     # must be a string '0101110001'
        else:
            self.label =gmpy2.digits(val,Configuration.base).zfill(Configuration.dimension)            


    def init_globals(base,N):
        """ define global class variables such as base and dimension"""
        Configuration.base      = base
        Configuration.dimension = N
    
    def get(self):
        """ return the configuration representation"""
        return self.label
    
    def get_state(self,k):
        """ return the state at node/position k held by configuration
        Ex: 
        10100  -> [0, 0, 1, 0, 1]. state[2] = 1 and state[4]=1
        """
        return self.label[k] # does not check if  0 < k < dimension - 1
    
    def get_dimension(self):
        """ returns number of nodes/elements for each representation"""
        return Configuration.dimension
    
    def set(self,conf):
        self.label = conf
        return None
    
    def set_state(self,k,new_state):
        """ change state[k] to new_state"""
        self.label[k] = new_state
        return self.label    
    
    def get_integer(self):
        """ returns the corresponding integer with base self.base"""
        return int(self.label,base)
    
    def get_count(self,state):
        """ returns the state sum over all nodes: sum(state[k],k=0..N-1)"""
        #return sum( int(x) for x in self.label )
        return self.label.count(state)

    def get_repetition(self):
        return None
    
    def get_representative(state):
        return state,int(state,Configuration.base),1

    def get_basis():
        """ returns the vector space"""
        return [Configuration(x) for x in range(Configuration.base**Configuration.dimension)]

    def __repr__(self):
        return '<Configuration %r>' % (self.label)

##==================================================



if __name__ == '__main__':

    N=3
    base=3
    Configuration.init_globals(base,N)
    for i in range(base**(N)):
        print(Configuration(i))

