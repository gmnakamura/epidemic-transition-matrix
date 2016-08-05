
import gmpy2


# import functools
# 
# def memoize(obj):
#     cache=obj.cache={}
#     @functools.wraps(obj)
#     def memoizer(*args,**kwargs):
#         if args not in cache:
#             cache[args]=obj(*args,**kwargs)
#         return cache[args]
#     return memoizer

class Memoize(object):
    def __init__(self,foo):
        self.foo = foo
        self.memo= {}
    def __call__(self, *args):
        if args not in self.memo:
            self.memo[args] = self.foo(*args)
        return self.memo[args]

class SymConf(object):
    """
    represents an unique network configuration.
    * Important: in this representation, states are limited 
                 up to maximum of 10 states.
    """
    dimension=1
    base     =2
    basemax  =2
    
    def __init__(self,val=0):
        """ initialize the class. val may assume either string or
        integer representation. Integer decomposition is performed
        by gmpy2.
        """
        tmp=val
        if isinstance(val,str):
            tmp = int(val,SymConf.base)     
        self.label,self.label_int,self.repetition = SymConf.get_representative(tmp)
    
    def get(self):
        """ return the representative configuration label"""
        return self.label
    def get_state(self,k):
        """ return the state at node/position k held by configuration
        Ex: 
        10100  -> [0, 0, 1, 0, 1]. state[2] = 1 and state[4]=1
        """
        return self.label[k] # does not check if  0 < k < dimension - 1
    def get_dimension(self):
        """ returns number of nodes/elements for each representation"""
        return SymConf.dimension
    
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

    def get_repetition(self):
        """ returns the number of times the configuration repeats itself
        after N cyclic permutations"""
        return self.repetition
    
    def get_count(self,state):
        """ returns the state sum over all nodes: sum(state[k],k=0..N-1)"""
        #return sum( int(x) for x in self.label )
        return self.label.count(state)
    

    def init_globals(base,N):
        """ define global class variables such as base and dimension"""
        SymConf.base      = base
        SymConf.dimension = N
        SymConf.basemax   = base**N

    
    @Memoize
    def get_representative(val,p=0):
        """ return the representative configuration of val
        with eigennumber p.

        This function explores the correspondence between 
        finite groups {O} and permutation operations. 
        Given a single state S (for instance '0101'), the
        linear combination with eigenvalue p is obtained 
        by (1/Normalization) sum_k  ([p Ô]^k S.

        We define the representative configuration as the 
        state with lowest integer 

        Ex: (base 2)

        state    representative (p=0 and N=4)
        1000     0001
        0100     0001
        
        0101     0101
        1010     0101

        """

        #first let us consider only p=0 (no null norm)
        if isinstance(val,str):
            val = int(val,SymConf.base)
        rep=val                
        current = val
        repetition=1
        for k in range(SymConf.dimension-1):            
            new  = current * SymConf.base
            shift,current = divmod(new,SymConf.basemax)
            current = current + shift
            if not (current > rep):
                repetition = repetition + max(0, 1-abs(rep-current))
                rep = current

        return gmpy2.digits(rep,SymConf.base).zfill(SymConf.dimension),rep,repetition

    # def get_all_representatives(p=0):
    #     lookup_int = {}
    #     lookup_str = {}
    #     for k in range(SymConf.basemax):
    #         x=SymConf.get_representative(k,p)
    #         lookup_str[k] = x[0]
    #         lookup_int[k] = x[1]
    #     reps_str=set(lookup_str.values())
    #     reps_int=set(lookup_int.values())        
    #     return lookup_str,reps_str,lookup_int,reps_int

    def get_basis(p=0):
        """ return the vector space for p-sector"""
        lookup = []
        for k in range(SymConf.basemax):
            lookup.append( SymConf.get_representative(k,p)[1])
        return [SymConf(x) for x in set(lookup)]

    def __repr__(self):
        return '<SymConf %r>' % (self.label)


##==================================================



if __name__ == '__main__':

    N=4
    base=2
    SymConf.init_globals(base,N)

    from timeit import default_timer as timer

    t=timer()
    for i in range(base**(N)):
        #print(SymConf(i))
        SymConf(i)
    print('Elapsed time without memoization %s' % (timer()-t))
    
    t=timer()
    for i in range(base**(N)):
        #print(SymConf(i))
        SymConf(i)
    print('Elapsed time with    memoization %s' % (timer()-t))    

    import sys
    print(sys.getsizeof(SymConf))
    print(sys.getsizeof(SymConf(0)))

    for i in range(base**N):
        print(i, SymConf(i).label_int )


    a=SymConf.get_basis()

    
