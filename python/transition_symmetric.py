from math import *
import scipy
from scipy import linalg

class Transition(object):
    """ Transition matrix.
    
    Describes the temporal evolution of probability vectors
    given one and two-two body transition rules. 
    """
    def __init__(self,_class,rules={}):
        self.rules=rules
        self._class=_class # define the class of vector object (symconf, configuration)
        #self.elements={ n:{} for n in range(_class.dimension*(_class.base-1)+1) }
        self.elements={}
        self.symmetries={}

    def set_rules(self,rules=None):
        """
        rules variable is a dictionary containing transitions.
        each key describes the current state(s) while the value
        express the respective (outcome state and coupling) 
        via tuple.

        rules may include only non-diagonal contributions.

        Ex. (SIS/SI)
           rules={'1':[('0',beta)], '01':[('11',alpha)], '10':[('11',alpha)]}
        """
        self.rules=rules # depends on problem


    # non-diagonal one body operator action
    # must be called before two-body
    def onebody(self,conf,output={}):
        """ returns a dictionary with all transitions from
        conf via one-body transformations."""
        cconf=conf.get()
        repc =conf.get_repetition() or 1.0
        repc =1.0/repc
        for k in range(conf.get_dimension()):
            for rule in (self.rules.get(cconf[k]) or []):
                trial =cconf[:k]+rule[0]+cconf[k+1:]
                (trial,repres,rep)=self._class.get_representative(trial)                
                output[trial]=rule[1]*sqrt(rep*repc)+(output.get(trial) or 0)
        return output

    def twobody(self,conf,A,output={}):
        """ returns a dictionary with all transitions from
        conf via two-body transformations.
        
        Nodes connections are accounted via adjacency matrix A (NxN)
        """
        cconf=conf.get()
        size =conf.get_dimension()
        repc =conf.get_repetition() or 1.0
        repc =1.0/repc
        for j in range(size):
            sj     =cconf[j]
            trialj =cconf[:j]            
            # run over all remaining nodes
            for i in range((j+1),size): 
                # pair = two-node state string (key for self.rules).
                pair = sj+cconf[i]
                # check if pair produces non-diagonal transitions
                for rule in (self.rules.get(pair) or []):
                    # new configuration
                    trial=trialj+rule[0][0]+cconf[j+1:i]+rule[0][1]+cconf[i+1:]
                    a=trial
                    (trial,repres,rep)=self._class.get_representative(trial)

                    # if trial not in output: # ***bottleneck***                        
                    #     output[trial]=0
                    # # update coupling
                    # output[trial] += rule[1]*A[j][i]*sqrt(rep*repc)

                    output[trial]=rule[1]*A[j][i]*sqrt(rep*repc)+(output.get(trial) or 0)
        return output

    def diagonal(self,conf,A):
        contrib=1.0
        cconf=conf.get()
        size=conf.get_dimension()
        for k in range(size):
            for rule in (self.rules.get(cconf[k]) or []):
                contrib=contrib - rule[1]
                print(conf,rule,' one',contrib)
                
        for j in range(size):
            sj     =cconf[j]
            for i in range((j+1),size): 
                pair = sj+cconf[i]
                for rule in (self.rules.get(pair) or []):
                    contrib=contrib - rule[1]*A[j][i]
                    print(conf,rule,' two',contrib)
        return contrib
    
    def compute_column(self,conf,A):
        """ returns all available transitions for a given
        state conf with adjacency matrix A.

        Obs: onebody and twobody may be computed independently.
        This suggests a simple parallel solution.
        1) Carolina task #1
        """
        output={}
        #non-diagonal contributions
        output= self.onebody(conf,output)   #1) use output1 = ..
        output= self.twobody(conf,A,output) #2) use output2 = ..  merge 
        #output[conf.get()]=1-sum( output.values() )
        output[conf.get()]=self.diagonal(conf,A)+(output.get(conf.get()) or 0.0)
        return output


    def compute(self,A):
        """ returns all matrix elements. Data is organized as follows:
        
        T.elements[n = eigennumber] [starting configuration] = 
                 ={ out_conf_0:coup_0, out_conf_1 : coup_1,...}
        """        
        for element in (self._class.get_basis() or []):            
            # n=element.get_count('1')
            # self.elements[n][element.get()]=self.compute_column(element,A)
            self.elements[element.get()]=self.compute_column(element,A)
        return None

    def get_matrix(self):
        """ convert self.elements into numpy array"""
        n=len(self.elements.keys())
        matrix = scipy.zeros(shape=(n,n))
        index  = {x:y for x,y in zip(self.elements.keys(),range(n))}
        for state in self.elements.keys():
            j = index[state]
            for out in self.elements[state]:
                i = index[out]
                matrix[j][i]=self.elements[state][out]
        return matrix.transpose(),index

    def get_eigenvals(self):
        return scipy.linalg.eigvals(self.get_matrix()[0])

if __name__ == "__main__":
    
    import sys
    from timeit import default_timer as timer
    from configuration import Configuration

    beta =0.001
    alpha=0.100
    rules={'1':[('0',beta)], '01':[('11',alpha)], '10':[('11',alpha)]}
    base =2
    N    =4

    
    Nmax=20
    Nmin=2
    size=[0]*Nmax
    time=[0]*Nmax
    

    import resource
    print('Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

    t1=timer()
    Configuration.init_globals(base,N)
    adjacency=[ [1]*N for i in range(N) ]
    for k in range(N):
        adjacency[k][k]=0
    T=Transition(Configuration,rules)
    T.compute(adjacency)
    print('Elapsed time: %s' % (timer()-t1))

    print('Memory usage: %s (kb)' % resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    print(sys.getsizeof(T.elements))
    print(sys.getsizeof(T))


    from symmetric_configuration import *
    SymConf.init_globals(base,N)
    TS=Transition(SymConf,rules)
    TS.compute(adjacency)
    
