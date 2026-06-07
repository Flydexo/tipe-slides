class Sponge():
    def __init__(self, rate, capacity, permutation, r_f, r_p, rc):
        self.r = rate
        self.c = capacity
        self.pi = permutation
        self.state = None
        self.r_f = r_f
        self.r_p = r_p
        self.mds = generate_cauchy_matrix(rate+capacity)
        self.rc = ARC
        self.trace = []

    def apply(self,data):
        self.state = np.array([FieldElement(0) for i in range(0,self.c+self.r)])
        self.trace = [self.state]
        #remaining = self.pad(remaining)
        chunks = [data[self.r*i:self.r*(i+1)] for i in range(0,len(data)//self.r)]
        i = 0
        while(len(chunks) > 0):
            i+=1
            self.state, trace = self.pi(self.state + np.pad(chunks[0], (0,self.c), mode='constant'), self.r_f, self.r_p, self.mds, self.rc)
            chunks = chunks[1:]
            self.trace.extend(trace)
        return self.state[0:self.c], self.trace     


def generate_cauchy_matrix(t):
    x = [FieldElement(e) for e in range(0,t)]
    y = [FieldElement(e) for e in range(t, 2*t)]
    return np.array([[FieldElement(1)/(x[i]-y[j]) for j in range(0,len(y))] for i in range(0,len(x))])

def arc(state, rc):
    return np.array(state) + rc

def mix(state, mds):
    return mds @ np.array(state)    

def full_s_box(state):
    return np.array(state) ** 5

def partial_s_box(state):
    return [state[0] ** 5, *state[1:]]

def hades_round_permutation(state, r_f, r_p, mds, rc):
    trace = []
    for r in range(0,r_f+r_p):
        state = arc(state, rc[r])
        if(r_f/2 <= r < r_f/2+r_p):
            state = partial_s_box(state)
        else:
            state = full_s_box(state)
        state = mix(state, mds)
        trace.append(state)
    return state, trace

class Poseidon():
    def __init__(self):
        self.sponge = Sponge(7,1,hades_round_permutation,8,57, ARC)

    def hash(self, data):
        return self.sponge.apply(data)
