import ast
import constants
from constants import blowup_factor
from channel import Channel
from field import FieldElement
from merkle import verify_decommitment
from poseidon import generate_cauchy_matrix
from functools import reduce
from utils import Automata

proof_file = open("proof.txt", "r")
proof = ast.literal_eval(proof_file.read())
automata = Automata.from_json(open("automata.json", "r").read())

def parse_proof(method, parse = None):
    global v_current_channel_idx
    global v_proof
    #print("channel",proof[current_channel_idx], "method", method)
    v_current_channel_idx+=1
    if parse == None:
        return v_proof[v_current_channel_idx-1][len(method)+1:]
    else:
        return parse(v_proof[v_current_channel_idx-1][len(method)+1:])

channel = Channel(proof)
v_proof = channel.proof
v_c = Channel()
v_current_channel_idx = 0
PTRACE_DOMAIN_SIZE=parse_proof("send", lambda x: int(x))
POSEIDON_TRACE_EVAL_LENGTH=parse_proof("send", lambda x: int(x))
EVAL_DOMAIN_LENGTH=parse_proof("send", lambda x: int(x))
AUTOMATA_TRACE_LENGTH=parse_proof("send", lambda x: int(x))
v_hash = v_proof[v_current_channel_idx][5:]
v_c.send(v_hash)
v_current_channel_idx+=1
v_state_merkle = v_proof[v_current_channel_idx][5:]
v_c.send(v_state_merkle)
v_current_channel_idx+=1
v_caracter_merkle = v_proof[v_current_channel_idx][5:]
v_c.send(v_caracter_merkle)
v_current_channel_idx+=1
v_transition_merkle = v_proof[v_current_channel_idx][5:]
v_c.send(v_transition_merkle)
v_current_channel_idx+=1
v_ptrace_merkle = v_proof[v_current_channel_idx][5:]
v_c.send(v_ptrace_merkle)
v_current_channel_idx+=1
v_hash_merkle = v_proof[v_current_channel_idx][5:]
v_c.send(v_hash_merkle)
v_current_channel_idx+=1
v_arc_merkle = v_proof[v_current_channel_idx][5:]
v_c.send(v_arc_merkle)
v_current_channel_idx+=1
v_alphas = []
for i in range(constants.CONSTRAINTS_LENGTH):
    v_a = FieldElement(int(v_proof[v_current_channel_idx][len('receive_random_field_element:'):]))
    v_current_channel_idx+=1
    v_alphas.append(v_a)
#a0 = FieldElement(int(proof[current_channel_idx][len('receive_random_field_element:'):]))
#assert(c.receive_random_field_element() == a0)
#current_channel_idx+=1
#a1 = FieldElement(int(proof[current_channel_idx][len('receive_random_field_element:'):]))
#assert(c.receive_random_field_element() == a1)
#current_channel_idx+=1
#a2 = FieldElement(int(proof[current_channel_idx][len('receive_random_field_element:'):]))
#assert(c.receive_random_field_element() == a2)
#current_channel_idx+=1
#a3 = FieldElement(int(proof[current_channel_idx][len('receive_random_field_element:'):]))
#assert(c.receive_random_field_element() == a3)
#current_channel_idx+=1
#a4 = FieldElement(int(proof[current_channel_idx][len('receive_random_field_element:'):]))
#assert(c.receive_random_field_element() == a4)
#current_channel_idx+=1
cp_merkle_root = v_proof[v_current_channel_idx][len('send:'):]
v_current_channel_idx+=1

w = FieldElement.generator()
h = FieldElement.generator() ** (3*2**30/(PTRACE_DOMAIN_SIZE*blowup_factor))
gp = FieldElement.generator() ** (3*2**30/(PTRACE_DOMAIN_SIZE))

#def x(i):
#    return gp**i

def x(i):
    return gp**(i%AUTOMATA_TRACE_LENGTH)

def z(i):
    return gp**(i%POSEIDON_TRACE_EVAL_LENGTH)

#x = [gp**i for i in range(AUTOMATA_TRACE_LENGTH)]
#z = [gp**i for i in range(POSEIDON_TRACE_EVAL_LENGTH)]
mds = generate_cauchy_matrix(constants.T)

def verify_FRI_commitment():
    global v_proof
    global v_current_channel_idx
    v_betas = []
    v_roots = [cp_merkle_root]
    v_l = 0
    while v_proof[v_current_channel_idx] != "send:FINISHED_FRI":
        v_betas.append(FieldElement(int(v_proof[v_current_channel_idx][len('receive_random_field_element:'):])))
        v_current_channel_idx+=1
        v_roots.append(v_proof[v_current_channel_idx][len('send:'):])
        v_current_channel_idx+=1
        v_l+=1
    v_current_channel_idx+=1
    v_poly0 = int(v_proof[v_current_channel_idx][len('send:'):])
    v_current_channel_idx+=1
    return v_betas,v_roots,v_poly0,v_l

v_betas,v_roots,v_poly0,v_l = verify_FRI_commitment()

print(v_l), 2**10

def verify_query_decommitment(l):
    global current_channel_idx
    for query in range(constants.N_QUERY):
        # r = int(proof[current_channel_idx][len('receive_random_field_element:'):])
        r = parse_proof("receive_random_int", int)
        # current_channel_idx+=1
        state_x = parse_proof("send", lambda x: FieldElement(int(x)))
        #state_x = FieldElement(int(proof[current_channel_idx][len('send:'):]))
        #current_channel_idx+=1
        state_x_path = parse_proof("send", lambda x: ast.literal_eval(x))
        state_x2 = parse_proof("send", lambda x: FieldElement(int(x)))
        #state_x = FieldElement(int(proof[current_channel_idx][len('send:'):]))
        #current_channel_idx+=1
        state_x2_path = parse_proof("send", lambda x: ast.literal_eval(x))
        #state_x_path = ast.literal_eval(proof[current_channel_idx][len('send:'):])
        #current_channel_idx+=1
        assert(verify_decommitment(r,state_x,state_x_path,v_state_merkle))
        assert(verify_decommitment((r+blowup_factor),state_x2,state_x2_path,v_state_merkle))

        #caracter_x = FieldElement(int(proof[current_channel_idx][len('send:'):]))
        #current_channel_idx+=1
        caracter_x = parse_proof("send", lambda x: FieldElement(int(x)))
        caracter_x_path = parse_proof("send", lambda x: ast.literal_eval(x))
        #caracter_x_path = proof[current_channel_idx][len('send:'):]
        #current_channel_idx+=1
        assert(verify_decommitment(r,caracter_x,caracter_x_path,v_caracter_merkle))

        transition_x = parse_proof("send", lambda x: FieldElement(int(x)))
        transition_x_path = parse_proof("send", lambda x: ast.literal_eval(x))
        assert(verify_decommitment((EVAL_DOMAIN_LENGTH*r+r),transition_x,transition_x_path,v_transition_merkle))

        poseidon_x = []
        poseidon_path = []
        poseidon_next_x = []
        poseidon_next_path = []
        hash_x = []
        hash_path = []
        arc_x = []
        arc_path = []
        
        for i in range(constants.T):
            poseidon_x.append(parse_proof("send", lambda x: FieldElement(int(x))))
            poseidon_path.append(parse_proof("send", lambda x: ast.literal_eval(x)))
            assert(verify_decommitment((i*EVAL_DOMAIN_LENGTH + (r)),poseidon_x[-1],poseidon_path[-1],v_ptrace_merkle))

            poseidon_next_x.append(parse_proof("send", lambda x: FieldElement(int(x))))
            poseidon_next_path.append(parse_proof("send", lambda x: ast.literal_eval(x)))
            assert(verify_decommitment((i*EVAL_DOMAIN_LENGTH  + ((r + blowup_factor))),poseidon_next_x[-1],poseidon_next_path[-1],v_ptrace_merkle))

            hash_x.append(parse_proof("send", lambda x: FieldElement(int(x))))
            hash_path.append(parse_proof("send", lambda x: ast.literal_eval(x)))
            assert(verify_decommitment((i*EVAL_DOMAIN_LENGTH  + (r)),hash_x[-1],hash_path[-1],v_hash_merkle))

            arc_x.append(parse_proof("send", lambda x: FieldElement(int(x))))
            arc_path.append(parse_proof("send", lambda x: ast.literal_eval(x)))
            assert(verify_decommitment((i*EVAL_DOMAIN_LENGTH  + (r)),arc_x[-1],arc_path[-1],v_arc_merkle))
        
        #transition_x = FieldElement(int(proof[current_channel_idx][len('send:'):]))
        #current_channel_idx+=1
        #transition_x_path = proof[current_channel_idx][len('send:'):]
        #current_channel_idx+=1

        layers = []
        layer_paths = []

        #length = 2**(v_l+1)
        length = EVAL_DOMAIN_LENGTH
        cps = []
            
        for j in range(0,l):
            idx = r % length
            sib_idx = (r + length // 2) % length
            layer_x = parse_proof("send", lambda x: FieldElement(int(x)))
            #layer_x = int(proof[current_channel_idx][len('send:'):])
            #current_channel_idx+=1
            layer_x_path = parse_proof("send",  lambda x: ast.literal_eval(x))
            assert(verify_decommitment(idx,layer_x,layer_x_path,v_roots[j]))
            #layer_x_path = proof[current_channel_idx][len('send:'):]
            #current_channel_idx+=1
            layer_sibx = parse_proof("send", lambda x: FieldElement(int(x)))
            #layer_sibx = int(proof[current_channel_idx][len('send:'):])
            #current_channel_idx+=1
            #layer_sibx_path = proof[current_channel_idx][len('send:'):]
            layer_sibx_path = parse_proof("send",  lambda x: ast.literal_eval(x))
            assert(verify_decommitment(sib_idx,layer_sibx,layer_sibx_path,v_roots[j]))
            e = (w*h**idx)**(2**j)
            
            g_x = (layer_x + layer_sibx)/FieldElement(2)
            h_x = (layer_x - layer_sibx)/(2*e) 
            cp_ip1 = g_x + v_betas[j]*h_x
            cps.append(cp_ip1)
            
            if j == 0:
                initial_constraint = (state_x-automata.state_map['q0'])/(e-1)
                Z_G = FieldElement(1)
                for i in range(AUTOMATA_TRACE_LENGTH):
                  Z_G = Z_G * (e-x(i))
                ap = FieldElement(1)
                for c in automata.alphabet:
                    ap = ap * (caracter_x - automata.caracter_map[c])
                caracter_constraint = ap/Z_G
                ap = FieldElement(1)
                for f in automata.accept_states:
                  ap = ap * (state_x - automata.state_map[f])
                accepting_constraint = ap/(e-x(-1))
                step_constraint = (state_x2 - transition_x)/(Z_G/(e-x(-1)))
                X_G = FieldElement(1)
                Y_G = FieldElement(1)
                for i in range(POSEIDON_TRACE_EVAL_LENGTH-1):
                    round_x = i%(constants.r_f+constants.r_p)
                    if(round_x < constants.r_f/2 or round_x >= constants.r_f/2+constants.r_p):
                      X_G = X_G * (e-z(i))
                    else:
                      Y_G = Y_G * (e-z(i))
                v_poseidon_full_round = mds@[(poseidon_x[j] + hash_x[j] + arc_x[j])**5 for j in range(constants.T)]
                v_poseidon_partial_round = mds@[(poseidon_x[j] + hash_x[j] + arc_x[j])**(5 if j == 0 else 1) for j in range(constants.T)]
                poseidon_full_round_cr = [(poseidon_next_x[j] - v_poseidon_full_round[j])/X_G for j in range(constants.T)]
                poseidon_partial_round_cr = [(poseidon_next_x[j] - v_poseidon_partial_round[j])/Y_G for j in range(constants.T)]
                poseidon_initial_cr = [poseidon_x[j]/(e-1) for j in range(constants.T)]
                poseidon_output_cr = (poseidon_x[0]-FieldElement(int(v_hash)))/(e-z(-1))

                v_ps = [poseidon_output_cr, initial_constraint, caracter_constraint, accepting_constraint, step_constraint, *poseidon_initial_cr, *poseidon_full_round_cr, *poseidon_partial_round_cr]
                cp_e = reduce(lambda x, y: x+y,map(lambda a, p: a*p, v_alphas, v_ps),FieldElement(0))
                assert(cp_e == layer_x)
            else:
                assert(cps[j-1] == layer_x)
            
            layers.append((layer_x,layer_sibx))
            layer_paths.append((layer_x_path, layer_sibx_path))
            length //= 2

        last_layer_x = parse_proof("send", int)
        assert(cps[-1] == last_layer_x)

verify_query_decommitment(v_l)

print_success_banner(verification_start, channel, EVAL_DOMAIN_LENGTH, POSEIDON_TRACE_EVAL_LENGTH)
