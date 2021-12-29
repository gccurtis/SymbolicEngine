import sympy
from itertools import combinations
sympy.init_printing(use_unicode=True)

# The goal of this script is to generate the 
# matrix form of a Q[r_1, r_2, ... r_n] with the respective commutative rules

# basis = ["i","j", ...] This will be the variables that represent r_1 ... r_n

# power_rules = [{"i":"1-1i","j":"-1"},{"k":"i"}, ...] 
# power_rules represent the equations that r_i solves
# each dictionary represents a polynomial equation whose order is two more than its index in the list
# example based on the definition: i^2 = 1-1*i, j^2 = -1, k^3 = i

# commutative_rules = {"ji" : -ij} this represents the transformation from a noncanonical term to a canonical one

# create_matrix_from_rules(["a","b","c","d","e","f","g"],commutative_rules, power_rules, basis)
# This will create the matrix using ["a","b","c","d",...] as the list of variables the matrix will use

# This will delineate 'string' by the return values of 'fun' evaluated on each character
# The delineation is done by grouping the chars that have the same return value
def string_groupby(string, fun):
    if string == "":
        return [""]
    return_list = []
    prev = fun(string[0])
    current_word = ""
    for char in string:
        val = fun(char)
        if val == prev:
            current_word += char
        else:
            prev = val
            return_list.append(current_word)
            current_word = char
    return_list.append(current_word)
    return return_list

# The same as isnumeric but returns true for negative numbers
def string_isnumeric(word):
    return word.isnumeric() or (len(word) > 1 and word[0] == "-" and word[1:].isnumeric())

# Returns the string 'string' but without the first occurence of the substring 'sub' 
def string_remove(string, sub):
    return string[:string.find(sub)]+string[string.find(sub)+len(sub):]

# Returns the string 'string' but without the substring 'sub'
def string_removeall(string, sub):
    s = string
    for i in range(string.count(sub)):
        s = string_remove(s, sub)
    return s

# This attempts to return a list of lists representing the mathematical terms of the string
# example: "-1" -> "-1", "1-i" = ["1","-i"], "1+(2+2j)" -> ["1",["2", "2j"]]
def separate_terms(string):
    terms = []
    current_word = ""
    i = 0
    p_level = 0
    for char in string:
        if char == "+" and p_level == 0 :
            terms.append(current_word)
            current_word = ""
        elif char == "-" and p_level == 0 :
            if current_word != "":
                terms.append(current_word)
            current_word = "-"
        elif char == "(" and p_level == 0:
            current_word = "" #Not handling non-simplified expressions
            terms.append(separate_terms(string[i+1:]))
            p_level += 1
        elif char == "(" and p_level > 0:
            p_level += 1
        elif char == ")" and p_level == 0:
            break
        elif char == ")":
            p_level -= 1
            current_word = ""
        elif p_level == 0 and char not in ["+",")","(","-","*"]:
            current_word += char
        i += 1
    if current_word != "":
        terms.append(current_word)
    return terms

# This will take a list of terms 'li' and separates the vector and scalars 
# before simplifying the the expression 
def separate_terms2(li, basis_vector):
    terms = []
    i = 0
    while i < len(li):
        if type(li[i]) == str:
            terms.append(separate_special_terms(li[i], basis_vector))
        elif type(li[i]) == list:
            i += 1
            for scalar in li[i-1]:
                terms.append(separate_special_terms(scalar+li[i], basis_vector))
        i += 1
    return combine_term_pairs(terms)

# Takes a term 'string' and separates it into it vectors and scalars 
# "2ij" -> [["2"],"ij"], "(1+7)i" -> [["1+7"],"i"]
def separate_special_terms(string, basis_vectors):
    scalar_term = ""
    vector_term = ""
    for char in string:
        if char in basis_vectors:
            vector_term += char
        else:
            scalar_term += char
    if scalar_term == string:
        return [[scalar_term], "1"]
    else:
        return [[scalar_term], vector_term]

def apply_commutative_rules_to_vector_terms(vector_terms, commutative_rules, basis_vectors):
    vector_terms_so_far = ""
    svector_terms = vector_terms.replace("1","")
    for i in range(len(svector_terms)-1):
        if svector_terms[i:i+2] in commutative_rules:
            new_replacement_terms = separate_terms2(separate_terms(commutative_rules[svector_terms[i:i+2]]), basis_vectors)
            new_replacement_terms = [[new_term[0], vector_terms_so_far + new_term[1] + svector_terms[i+2:]] for new_term in new_replacement_terms]
            return new_replacement_terms
        else:
            vector_terms_so_far += svector_terms[i]
    if svector_terms == "":
        return [[["1"], "1"]]
    else:
        return [[["1"], svector_terms]]

# This applies the power rules
def apply_special_rules_to_vector_terms(vector_terms, all_special_rules, basis_vectors):
    vector_terms_so_far = ""
    svector_terms = vector_terms.replace("1","")
    for i in range(len(svector_terms)-1):
        k = 2
        for special_rules in all_special_rules:
            if svector_terms[i]*k == svector_terms[i:i+k] and svector_terms[i] in special_rules:
                new_replacement_terms = separate_terms2(separate_terms(special_rules[svector_terms[i]]), basis_vectors)
                new_replacement_terms = [[new_term[0], vector_terms_so_far + new_term[1] + svector_terms[i+k:]] for new_term in new_replacement_terms]
                return new_replacement_terms
            k += 1
        vector_terms_so_far += svector_terms[i]
        
    if svector_terms == "":
        return [[["1"], "1"]]
    else:
        return [[["1"], svector_terms]]

# This takes a list of scalar and vector parts of a term 'term_pair' and returns a it scaled by 'scalar'
# example: 2 ,[["2","j"],"i"] -> [["4","2j"],"i"]
def scale_term_pair(scalar, term_pair):
    return_terms = []
    for term1 in term_pair[0]:
        for term2 in scalar:
            scalar_vars = ""
            num = 1
            isnumeric = lambda x: (x.isnumeric() or x=="-")
            for word in string_groupby(term1, isnumeric) + string_groupby(term2, isnumeric):
                if string_isnumeric(word):
                    num *= int(word)
                
                else:
                    scalar_vars += word
            return_terms.append(str(num) + scalar_vars)
    return [return_terms, term_pair[1]]
            
# Returns a list  will simplify term_pairs that have the same vectors
# [ [["2","j"],"i"], [["1"],"i"] ] -> [ [["3","j"],"i"] ]
def combine_term_pairs(list_of_term_pairs):
    combined_dict = dict()
    combined_list = []
    for term_pair in list_of_term_pairs:
        if term_pair[1] in combined_dict:
            if term_pair[0][0] == "0":
                continue
            if string_isnumeric(term_pair[0][0]):
                if string_isnumeric(combined_dict[term_pair[1]][0]):
                    combined_dict[term_pair[1]][0] = str(int(combined_dict[term_pair[1]][0])+int(term_pair[0][0]))
                    combined_dict[term_pair[1]].extend(term_pair[0][1:])
                else:
                    combined_dict[term_pair[1]].reverse()
                    combined_dict[term_pair[1]].append(term_pair[0][0])
                    combined_dict[term_pair[1]].reverse()
                    combined_dict[term_pair[1]].extend(term_pair[0][1:])
            else:
                combined_dict[term_pair[1]].extend(term_pair[0])
        else:
            combined_dict[term_pair[1]] = term_pair[0]
    for key, val in combined_dict.items():
        combined_list.append([val, key])
    return combined_list

# Returns a list of term pairs that are the sum of term_pairs in both lists 
def sum_term_pairs(list_of_term_pairs1, list_of_term_pairs2):
    return combine_term_pairs(list_of_term_pairs1 + list_of_term_pairs2)

# Returns a list of term pairs that are the distributed term pairs of the 2 lists 
def distribute_term_pairs(list_of_term_pairs1, list_of_term_pairs2):
    return_list = []
    for term_pair1 in list_of_term_pairs1:
        for term_pair2 in list_of_term_pairs2:
            return_list.append([scale_term_pair(term_pair1[0], term_pair2)[0], term_pair1[1]+term_pair2[1]])
    return return_list

# Returns a list of term pairs that is the product of the lists of term pairs inputed following the commutative and power rules
# vector terms are defined in basis_vectors
def multiply_term_pairs(list_of_term_pairs1, list_of_term_pairs2, commutative_rules, power_rules, basis_vectors):
    distributed_list = distribute_term_pairs(list_of_term_pairs1, list_of_term_pairs2)
    simplified_list = []
    loop_list = distributed_list
    while True:
        simplified_list = []
        for term_pair in loop_list:
            #tmp_list is a list of term pairs after applying commutative rules
            tmp_list = [scale_term_pair(term_pair[0], x) for x in apply_commutative_rules_to_vector_terms(term_pair[1], commutative_rules, basis_vectors)]
            for tmp_term_pair in tmp_list:
                simplified_list.extend([scale_term_pair(tmp_term_pair[0], x) for x in apply_special_rules_to_vector_terms(tmp_term_pair[1], power_rules, basis_vectors)])
        simplified_list = combine_term_pairs(simplified_list)
        if loop_list == simplified_list:
            break
        else:
            loop_list = simplified_list
    return simplified_list

# Creates a list of term pairs that represent a canonical form element of the ring created by commutative_rules and power_rules 
def generate_multi_vector_from_rules(coefs, commutative_rules, power_rules, basis_vectors):
    basis = ["1"]
    for vector in basis_vectors:
        k = 1
        for rule in power_rules:
            basis.append(vector*k)
            k += 1
            if vector in rule:
                break
    basis_pairs = list(combinations(basis,2))
    mvector = [[[coefs[0]],"1"]]
    i = 1
    for pair in basis_pairs:
        if pair[0] == "1":
            vector_term = pair[1]
        elif pair[1] == "1":
            vector_term = pair[0]
        elif pair[0].count(pair[1])+pair[1].count(pair[0]) == 0:
            vector_term = pair[0]+pair[1]
        else:
            continue
        mvector.append([[coefs[i]], vector_term])
        i += 1
    return mvector

# Puts a * between units in a term
def preprocess_string_for_sympify(s):
    ret = ""
    for c in s:
        if c.isnumeric() or c.isalpha():
            ret += c+"*"
        else:
            ret = ret[:-1]+c
    if ret[-1] == "*":
        ret = ret[:-1]
    return ret
 
def create_matrix_from_rules(coefs, commutative_rules, power_rules, basis_vectors):
    apply_vector_coefs = ["x"*(i+1) for i in range(len(coefs))]
    origin_mvector = generate_multi_vector_from_rules(coefs, commutative_rules, power_rules, basis_vectors)
    apply_mvector = generate_multi_vector_from_rules(apply_vector_coefs, commutative_rules, power_rules, basis_vectors)
    origin_term_pairs = multiply_term_pairs(origin_mvector,apply_mvector, commutative_rules, power_rules, basis_vectors)
    print(apply_mvector)
    print(origin_term_pairs)
    organized_matrix = []
    for row in origin_term_pairs:
        tmp = [""]*len(origin_mvector)
        for c in row[0]:
            tmp[c.count("x")-1] += string_removeall(c,"x")+"+"
        print("tmp:",tmp)
        organized_matrix.append(list(map(lambda x: x[:-1], tmp)))
    #print(organized_matrix)
    sympy_matrix = []
    for r in organized_matrix:
        tmp = []
        for e in r:
            tmp.append(sympy.sympify(preprocess_string_for_sympify(e)))
        sympy_matrix.append(tmp)
    return sympy.Matrix(sympy_matrix)         

basis = ["i"]
power_rules = [{},{"i":"u"}]
commutative_rules = {}
m = create_matrix_from_rules(["a","b","c","d","e","f","g"],commutative_rules, power_rules, basis)
#m_alt = create_matrix_from_rules(["q","w","e","r"],commutative_rules, power_rules, basis)
print(m)

'''
em = sympy.exp(m)
dim = int(len(em)**(1/2))
#genexp = sympy.Matrix([em[i*dim].simplify().subs({sympy.Symbol("a"): sympy.Symbol("x"), sympy.Symbol("b"): sympy.Symbol("y")}) for i in range(dim)])
genexp = sympy.Matrix([em[i*dim].simplify().subs(
    {sympy.Symbol("a"): sympy.Symbol("x"),
     sympy.Symbol("b"): sympy.Symbol("y"),
     sympy.Symbol("c"): sympy.Symbol("z"),})
                       for i in range(dim)])
unit_scalar = m.det().subs({sympy.Symbol("a"): genexp[0],
                        #sympy.Symbol("b"): genexp[1]}).simplify()     
                        sympy.Symbol("b"): genexp[1],
                        sympy.Symbol("c"): genexp[2]}).simplify()
'''


tosolve_vars = [sympy.Symbol("x"),sympy.Symbol("y"),sympy.Symbol("z")]
tosolve = sympy.Matrix(tosolve_vars)
tosolve = m*tosolve
result = sympy.Matrix([m.det(), *[0 for i in range(len(tosolve_vars)-1)]])
eqs = [sympy.Eq(tosolve[i],result[i]) for i in range(len(tosolve_vars))]
conj = sympy.solve(eqs,*tosolve_vars)
conj = sympy.Matrix([conj[s] for s in tosolve_vars])

