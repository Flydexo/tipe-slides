#word = open("text.txt", "r").read()
word = "I think STARKS are fun to learn"

import pandas as pd
from itertools import product
import time
from poseidon import Poseidon, ARC, generate_cauchy_matrix
from field import FieldElement
import numpy as np
from polynomial import Polynomial, interpolate_poly
from channel import Channel
from merkle import MerkleTree, verify_decommitment
from functools import reduce
from constants import N_QUERY, r_f, r_p, blowup_factor, T
import constants
from utils import Automata, primes
from tqdm import tqdm

start_all = time.time()

automata = Automata.from_json(open("automata.json", "r").read())
trace = automata.generate_trace(word)

mds = generate_cauchy_matrix(T)

hash_input = [FieldElement(x) for x in np.pad([automata.caracter_map[x] for x in word], (0, (-len(word)) % (T-1)))]

phash, ptrace = Poseidon().hash(hash_input)

X = Polynomial.X()

def next_power_of_two(n: int) -> int:
    if n <= 0:
        return 1
    # If n is already a power of two, return n
    if (n & (n - 1)) == 0:
        return n
    power = 1
    while power < n:
        power <<= 1
    return power

def hash_input_index(i,j):
    if (i % (r_f+r_p) == 0 and j < 7 and i//(r_f+r_p) < len(hash_input)//7):
        return hash_input[7*(i//(r_f+r_p))+j]
    else:
        return FieldElement(0)

ptrace_domain_size = next_power_of_two(len(ptrace))
gp = FieldElement.generator() ** (3*2**30/(ptrace_domain_size))
x = [gp**i for i in range(len(trace['i']))]
z = [gp**i for i in range(len(ptrace))]

state_polynomial = interpolate_poly(x, [FieldElement(automata.state_map[s]) for s in trace['s']])
caracter_polynomial = interpolate_poly(x, [FieldElement(automata.caracter_map[c]) for c in trace['c']])
transi_x = [FieldElement(automata.state_map[s]*automata.caracter_map[c]) for s,c in tqdm(product(automata.states, automata.alphabet))]

transition_polynomial = interpolate_poly(transi_x, [FieldElement(automata.state_map[automata.transition[s][c]]) for (s,c) in tqdm(product(automata.states, automata.alphabet))])

step_polynomial = state_polynomial(gp*X) - transition_polynomial(state_polynomial(X)*caracter_polynomial(X))

ptrace_polynomials = [interpolate_poly(z, [row[i] for row in ptrace]) for i in range(0,len(ptrace[0]))]
arc_polys = [interpolate_poly(z, [FieldElement(ARC[j%(r_f+r_p)][i]) for j in range(len(ptrace))]) for i in range(len(ARC[0]))]
hash_input_poly = [interpolate_poly(z, [hash_input_index(i,j) for i in range(len(z))]) for j in range(len(ptrace[0]))]

w = FieldElement.generator()
h = FieldElement.generator() ** (3*2**30/(ptrace_domain_size*blowup_factor))
H = [h**i for i in range(ptrace_domain_size*blowup_factor)]
eval_domain = [w*h for h in H]

state_eval = [state_polynomial(d) for d in tqdm(eval_domain)]
caracter_eval = [caracter_polynomial(d) for d in tqdm(eval_domain)]
#transition_eval = [transition_polynomial(state*caracter) for state,caracter in tqdm(product(state_eval,caracter_eval))]
transition_eval = [transition_polynomial(state*caracter) for state,caracter in zip(state_eval,caracter_eval)]
ptrace_evals = [[p(x) for x in tqdm(eval_domain)] for p in ptrace_polynomials]
hash_evals = [[p(x) for x in tqdm(eval_domain)] for p in hash_input_poly]
arc_polys_eval = [[p(x) for x in tqdm(eval_domain)] for p in arc_polys]

state_merkle = MerkleTree(state_eval)
caracter_merkle = MerkleTree(caracter_eval)
transition_merkle = MerkleTree(transition_eval)
ptrace_merkle = MerkleTree([e for p in ptrace_evals for e in p])
hash_merkle = MerkleTree([e for p in hash_evals for e in p])
arc_merkle = MerkleTree([e for p in arc_polys_eval for e in p])

initial_constraint_p = (state_polynomial-automata.state_map['q0'])/(X-1)
Z_G = Polynomial([FieldElement.one()])
for i in range(len(x)):
  Z_G = Z_G * (X-x[i])
ap = Polynomial([FieldElement.one()])
for c in automata.alphabet:
    ap = ap * (caracter_polynomial - automata.caracter_map[c])
caracter_constraint_p = ap/Z_G
ap = Polynomial([FieldElement.one()])
for f in automata.accept_states:
  ap = ap * (state_polynomial - automata.state_map[f])
accepting_constraint_p = ap/(X-x[-1])
step_constraint_p = step_polynomial/(Z_G/(X-x[-1]))
poseidon_full_round_polynomials = mds@[(ptrace_polynomials[j](X) + hash_input_poly[j](X) + arc_polys[j](X))**5 for j in range(T)]
poseidon_partial_round_polynomials = mds@[(ptrace_polynomials[j](X) + hash_input_poly[j](X) + arc_polys[j](X))**(5 if j == 0 else 1) for j in range(T)]
X_G = Polynomial([FieldElement.one()])
Z_G = Polynomial([FieldElement.one()])
for i in range(len(ptrace)-1):
    r = i%(r_f+r_p)
    if(r < r_f/2 or r >= r_f/2+r_p):
      X_G = X_G * (X-z[i])
    else:
        Z_G = Z_G * (X-z[i])
poseidon_constraint_full_round = [(ptrace_polynomials[j](gp*X) - poseidon_full_round_polynomials[j](X))/X_G for j in range(len(ptrace_polynomials))]
poseidon_constraint_partial_round = [(ptrace_polynomials[j](gp*X) - poseidon_partial_round_polynomials[j](X))/Z_G for j in range(len(ptrace_polynomials))]
poseidon_constraints_initial = [p/(X-1) for p in ptrace_polynomials]
poseidon_output_constraint = (ptrace_polynomials[0](X)-phash[0])/(X-z[-1])

channel = Channel()
channel.send(str(ptrace_domain_size))
channel.send(str(len(ptrace)))
channel.send(str(len(eval_domain)))
channel.send(str(len(trace['i'])))
channel.send(str(phash[0]))
channel.send(state_merkle.root)
channel.send(caracter_merkle.root)
channel.send(transition_merkle.root)
channel.send(ptrace_merkle.root)
channel.send(hash_merkle.root)
channel.send(arc_merkle.root)

ps = [poseidon_output_constraint, initial_constraint_p, caracter_constraint_p, accepting_constraint_p, step_constraint_p, *poseidon_constraints_initial, *poseidon_constraint_full_round, *poseidon_constraint_partial_round]
alphas = [channel.receive_random_field_element() for _ in range(len(ps))]
cp = reduce(lambda x, y: x+y,map(lambda a, p: a*p, alphas, ps),FieldElement(0))

cp_eval = [cp(x) for x in eval_domain]
cp_merkle = MerkleTree(cp_eval)
channel.send(cp_merkle.root)

def fold_domain(domain):
  return [x**2 for x in domain[:len(domain)//2]]

def fold(p, beta):
  odd = p.poly[1::2]
  even = p.poly[::2]
  o = Polynomial(odd)
  e = Polynomial(even)
  return e + beta*o

def next_fri_layer(poly, domain, beta):
    next_poly = fold(poly, beta)
    next_domain = fold_domain(domain)
    next_layer = [next_poly(x) for x in next_domain]
    return next_poly, next_domain, next_layer

cp.degree(), len(cp_eval)

def FriCommit(cp, domain, cp_eval, cp_merkle, channel):
    fri_polys = [cp]
    fri_domains = [domain]
    fri_layers = [cp_eval]
    fri_merkles = [cp_merkle]
    betas = []
    i = 0
    while fri_polys[-1].degree() > 0:
        beta = channel.receive_random_field_element()
        betas.append(beta)
        next_poly, next_domain, next_layer = next_fri_layer(fri_polys[-1], fri_domains[-1], beta)
        fri_polys.append(next_poly)
        fri_domains.append(next_domain)
        fri_layers.append(next_layer)
        fri_merkles.append(MerkleTree(next_layer))
        channel.send(fri_merkles[-1].root)
        i+=1
    channel.send("FINISHED_FRI")
    channel.send(str(fri_polys[-1].poly[0]))
    return fri_polys, fri_domains, fri_layers, fri_merkles, betas

fri_polys, fri_domains, fri_layers, fri_merkles, betas = FriCommit(cp, eval_domain,cp_eval,cp_merkle,channel)

len(fri_layers[-1])

def decommit_on_fri_layers(idx, channel):
    prev_idx = None
    prev_sibidx = None
    p_idx = None
    i = 0
    for layer, merkle in zip(fri_layers[:-1], fri_merkles[:-1]):
        length = len(layer)
        idx = idx % length
        sib_idx = (idx + length // 2) % length
        channel.send(str(layer[idx]))
        channel.send(str(merkle.get_authentication_path(idx)))
        channel.send(str(layer[sib_idx]))
        channel.send(str(merkle.get_authentication_path(sib_idx)))
    channel.send(str(fri_layers[-1][0]))

def decommit_on_query(idx, channel): 
    assert idx + blowup_factor < len(eval_domain), f'query index: {idx} is out of range. Length of layer: {len(eval_domain)}.'
    channel.send(str(state_eval[idx%len(eval_domain)])) # f(x).
    channel.send(str(state_merkle.get_authentication_path(idx%len(eval_domain)))) # auth path for f(x).
    assert(state_eval[(idx + blowup_factor)%len(eval_domain)] == state_polynomial(gp*eval_domain[idx%len(eval_domain)]))
    channel.send(str(state_eval[(idx + blowup_factor)%len(eval_domain)])) # f(gx).
    channel.send(str(state_merkle.get_authentication_path((idx + blowup_factor)%len(eval_domain)))) # auth path for f(gx).
    assert (idx + blowup_factor)%len(eval_domain) < len(caracter_eval), f'query index: {idx%len(eval_domain)} is out of range. Length of layer: {len(caracter_eval)}.'
    channel.send(str(caracter_eval[idx%len(eval_domain)])) # f(x).
    channel.send(str(caracter_merkle.get_authentication_path(idx%len(eval_domain)))) # auth path for f(x).
    assert (idx + blowup_factor)%len(eval_domain) < len(transition_eval), f'query index: {idx%len(eval_domain)} is out of range. Length of layer: {len(transition_eval)}.'
    transi_x = transition_polynomial(caracter_eval[idx%len(eval_domain)]*state_eval[idx%len(eval_domain)])
    channel.send(str(transi_x)) # f(x).
    transi_idx = (len(eval_domain)*idx+idx)%len(eval_domain)
    channel.send(str(transition_merkle.get_authentication_path(transi_idx))) # auth path for f(x).


    #### POSEIDON ####
    for row in range(len(ptrace[0])):
        trace_x = ptrace_evals[row][idx%len(eval_domain)]
        channel.send(str(trace_x))
        channel.send(str(ptrace_merkle.get_authentication_path((len(eval_domain)*row+(idx%len(eval_domain))))))
        assert(ptrace_polynomials[row](eval_domain[idx%len(eval_domain)]) == trace_x)

        next_trace_x = ptrace_evals[row][(idx+blowup_factor)%len(eval_domain)]
        channel.send(str(next_trace_x))
        channel.send(str(ptrace_merkle.get_authentication_path((len(eval_domain)*row+((idx+blowup_factor)%len(eval_domain))))))
        assert(ptrace_polynomials[row](gp*eval_domain[idx%len(eval_domain)]) == next_trace_x)

        hash_x = hash_evals[row][idx%len(eval_domain)]
        channel.send(str(hash_x))
        channel.send(str(hash_merkle.get_authentication_path((len(eval_domain)*row+(idx%len(eval_domain))))))
        
        arc_x = arc_polys_eval[row][idx%len(eval_domain)]
        channel.send(str(arc_x))
        channel.send(str(arc_merkle.get_authentication_path((len(eval_domain)*row+(idx%len(eval_domain))))))
        
    decommit_on_fri_layers(idx, channel)

def decommit_fri(channel):
    for query in range(N_QUERY):
        decommit_on_query(channel.receive_random_int(0, len(eval_domain)), channel)

decommit_fri(channel)

proof_file = open("proof.txt", "w")
proof_file.write(str(channel.proof))

print_prover_success_banner(start_all, channel, len(ptrace), len(state_eval))
