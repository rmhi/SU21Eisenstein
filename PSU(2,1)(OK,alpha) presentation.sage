"""
This code was written by Richard M. Hill,
Department of Mathematics, University College London.
It is available at https://github.com/rmhi/SU21Eisenstein
It runs on sagemath version 9.0.
"""


K.<omega>=NumberField(x^2+x+1)
OK = K.ring_of_integers()
omega = OK(omega)     # omega is a primitive cube root of unity
alpha = 2*omega+1     # alpha is sqrt(-3)
mu_3 = [1,omega,omega^2]
mu_6 = [(1+omega)^n for n in range(6)]

MatOK = MatrixSpace(OK,3,3)
VOK = OK^3
VK = K^3
F3.<one_mod,omega_mod> = OK.quotient_ring(OK.ideal([alpha])) # the finite field OK/(alpha) with 3 elements. 
MatF3 = MatrixSpace(F3,3,3)

def K2CC(zK):         # convert elements of K into complex numbers.
    a,b = list(K(zK))
    return CC(a+b*exp(2*pi*i/3))

J = MatOK([[0,0,1],[0,1,0],[1,0,0]])

def Form(v1,v2):      # the Hermitian form corresponding to J.
    return (v1.transpose()*J*v2.conjugate())[0,0]

def N(z,x):     # construct an upper triangular element in SU(2,1)(OK)
    if z in OK and x in ZZ and (OK(z).norm()+x)%2==0:
        return MatOK([[1,z,(-OK(z).norm()+x*alpha)/2],[0,1,-OK(z).conjugate()],[0,0,1]])
    else:
        print("Error in N(z,x)")
        print("You entered z=',z','and x=",x)
        raise ValueError
        
def Nt(z,x):     # construct a lower triangular element in SU(2,1)(OK)
    if z in OK and x in ZZ and (OK(z).norm()+x)%2==0:
        return N(z,x).transpose()
    else:
        print("Error in Nt(z,x)")
        print("You entered z=',z','and x=",x)
        raise ValueError

def in_SU(g):    # check whether a matrix is in SU(2,1)(OK)
    if g in MatOK:
        return MatOK(g).transpose()*J*MatOK(g).conjugate()==J and MatOK(g).det()==1
    else:
        print("Error in function in_SU.")
        print("The argument should be a matrix")
        print("Your argument is ",g," of type ", g.parent())
        raise ValueError

def in_SUalpha(g):
    return in_SU(g) and MatF3(g)==MatF3(1)

def in_U(g):    # check whether a matrix is in U(2,1)(OK)
    if g in MatOK:
        return g.transpose()*J*g.conjugate()==J
    else:
        print("Error in function in_SU.")
        print("The argument should be a matrix")
        print("Your argument is ",g," of type ", g.parent())
        raise ValueError

    

# The following is a set of generators of PSU(2,1)(OK,alpha).
# we shall always identify this group with a subgroup of SU(2,1)(OK,alpha), as
# SU(2,1)(OK,alpha) is a direct sum of this group with its centre <Z>
n1=N(alpha,1)
n2=N(omega*alpha,1)
n3=N(0,2)
n1t = n1.transpose()
n2t = n2.transpose()
n3t = n3.transpose()
generatorsPSUalpha = [n1,n2,n3,n1t,n2t,n3t]


# This is a projection homomorphism from SU(2,1)(OK,alpha) to PSU(2,1)(OK,alpha), regarded as a subgroup
def SUalpha2PSUalpha(g):
    temp = g[0][0]
    for zeta in mu_3:
        x,y = list(temp-zeta)
        if x %3==0 and y%3==0:
            return MatOK(zeta)^-1 * g
    print('Error in SUalpha2PSUalpha. No cube root of unity is congruent to g[0][0] modulo 3.')
    return None


# Falbel - Parker generators and relations.
# These are generators of PU(2,1)(OK), i.e. U(2,1)(OK) modulo its centre <Z>.
FalPar_P = diagonal_matrix([1,omega,1])*N(1,1)
FalPar_Q = diagonal_matrix([1,-1,1])*N(1,1)
FalPar_R = diagonal_matrix([1,-1,1])*J
FalPar_generators= {'P':FalPar_P, 'Q':FalPar_Q, 'R':FalPar_R}
FalPar_relations = [[FalPar_R,FalPar_R],
                    6*[FalPar_Q,FalPar_P^-1],
                    [FalPar_P,FalPar_Q^-1,FalPar_R,FalPar_Q,FalPar_P^-1,FalPar_R],
                   [FalPar_P,FalPar_P,FalPar_P,FalPar_Q^-1,FalPar_Q^-1],
                   3*[FalPar_R,FalPar_P]]





"""
This code finds a set of coset reps for U(2,1)(OK) / SU(2,1)(OK,alpha).
More precisely, it finds a Schreier transversal with respect to the Falbel-Parker generators.
Note that this quotient is isomorphic to PU(2,1)(OK) / PSU(2,1)(OK,alpha), and
SU(2,1)(OK,alpha) = PSU(2,1)(OK,alpha) x mu_3 (direct sum)

The data is stored as a list Transversal,
whose elements are triples
      (r,h,zeta),
where r is the coset representative,
      h is its reduction of r modulo alpha,
      zeta = det(r).
"""
print('Finding a Schreier transversal for PSU(2,1)(OK,alpha) in PU(2,1)(OK)')
print('with respect to the Falbel-Parker generators....')
Transversal=[(MatOK(1),MatF3(1),1)]
temp = [(MatF3(1),1)]
unfinished = True
while unfinished:
    unfinished=False
    for g in FalPar_generators.values():
        for elt in Transversal:
            test = MatOK(elt[0])*g
            test_mod = (MatF3(test),test.det())
            if not test_mod in temp:
                temp.append(test_mod)
                Transversal.append((test,MatF3(test),test.det()))
                unfinished=True
print('Found Schreier transversal of size', len(Transversal))

Transversal_dict = {str((y,z)):x for x,y,z in Transversal}


"""
Given an element g in U(2,1)(OK), this returns a pair (r,h),
where r is in Transversal,
      h is in SU(2,1)(OK)(alpha),
      g=h * r.
"""
def Schreier_decomp(g):
    g_mod = MatF3(g)
    g_det = g.det()
    key = str((g_mod,g_det))
    coset_rep = Transversal_dict[key]
    return [coset_rep,  g*coset_rep^-1]



#  Find generators in SU(2,1)(OK,alpha) and PSU(OK,alpha)
#  using the Schreier transversal.
print('Finding generators of PSU(2,1)(OK,alpha) using Schreier\'s lemma....')
RSgens_PSUalpha = []
RSgens_PSUalpha2 = {}
RSindices_PSUalpha = {}
current_index = 1
for g in FalPar_generators.values():
    for r in Transversal:
        test = SUalpha2PSUalpha(Schreier_decomp(r[0]*g)[1])
        if test != MatOK(1):
            if not(test in RSgens_PSUalpha or test^-1 in RSgens_PSUalpha):
                test_inv = test^-1
                RSgens_PSUalpha.append(test)
                RSgens_PSUalpha2[str(test)]=current_index
                RSgens_PSUalpha2[str(test_inv)]=-current_index
                RSindices_PSUalpha[current_index]= test
                RSindices_PSUalpha[-current_index]= test_inv
                current_index += 1
print("Generators of PSU(2,1)(OK,alpha) found from the Schreier transversal:",current_index-1)



def RSgen2index(g):
    if str(g) in RSgens_PSUalpha2:
        return RSgens_PSUalpha2[str(g)]
    else:
        print('Error in RSgen2index(g).')
        print('That element and its inverse are not in the list of generators.')
        pretty_print('g=',g)
        raise ValueError


# cancel an RS-word, so that it does not contain [x,-x]
# An RS-word is a list of integers, which are labels for RS-generators
# and their inverses.
# NOTE : CAN BE SPEEDED UP IF NECESSARY.
def reduce_RSword(symbolic_word):
    pruned_word = symbolic_word
    unfinished = True
    while unfinished:
        unfinished=False
        for x in range(len(pruned_word)-1):
            if pruned_word[x] == -pruned_word[x+1]:
                pruned_word = pruned_word[:x]+ pruned_word[x+2:]
                unfinished = True
                break
    return pruned_word



"""
Find relations in PSU(OK,alpha)
This is the Reidemeister - Schreier algorithm
"""
print('Finding relations in PSU(2,1)(OK,alpha) using the Reidemeister - Schreier algorithm....')
RSrelations = []
#symbolic_relns = []
RSsymbolic_relns = []
for r in FalPar_relations:
    #pretty_print(r)
    for starting_point in Transversal:
        current_point = starting_point[0]
        subgp_word=[]
        #symbolic_word=[]
        symbolic_word2=[]
        for g in r:
            current_point,skip = Schreier_decomp(current_point*g)
            PSUskip = SUalpha2PSUalpha(skip)
            if PSUskip != MatOK(1):
                subgp_word.append(PSUskip)
                #symbolic_word.append(generator2indexPSUalpha(PSUskip))
                symbolic_word2.append(RSgen2index(PSUskip))
        #if prod(r)!=MatOK(1):
        #    subgp_word.append(prod(r)^-1)
        # TO DO: ADD TO THE SYMBOLIC WORD HERE
        symbolic_word2 = reduce_RSword(symbolic_word2)
        if len(subgp_word)>2:
            RSrelations.append(subgp_word)
            #symbolic_relns.append(symbolic_word)
            RSsymbolic_relns.append(symbolic_word2)
            #print(symbolic_word)
        #print(prod(subgp_word))
print('Found',len(RSrelations),'relations.')



def inverse_word(symbolic_word):
    return [-x for x in reversed(symbolic_word)]



    

"""
return the nearest Eisenstein integer to q.
The difference z-q will be in a regular hexagon centred at 0
with two vertical sides distance 1 apart. This has corners at zeta/alpha,
where zeta runs through the 6th roots of unity.
"""
def NearestOK(q):
    try:
        y = (q-K(q).conjugate())/alpha
    except ValueError:
        print('ValueError in NearestOK(q).\n    q=',q,q.parent())
        raise ValueError
    try:
        y = QQ(y)
    except ValueError:
        print('ValueError in NearestOK(q).\n    y=',y,y.parent())
        raise ValueError
    y_int = round(y)
    temp = q - y_int * omega
    x_int = round(temp.trace()/2)
    z = x_int+y_int*omega
    for test in [z+omega,z+omega+1, z-omega,z-omega-1]:
        if (q-test).norm() < (q-z).norm():
            return test
    return z




# The function decomp2_generator_PSUalpha(g) takes an element g,
# which is either an upper or lower triangular matrix in PSU(2,1)(OK,alpha)
# and returns a word in the Euclidean generators.
# The word is a list of integers in {1,...,6,-1,...,-6}
def decomp2_generator_PSUalpha(g):
    #pretty_print(g)
    if g[0,2]==0:
        #print("lower triangular")
        z = g[1][0]/alpha
        a,b = z
        #print("z=",z,"a=",a,"b=",b, (a+omega* b)*alpha)
        h = n2t^-b*n1t^-a * g
        #pretty_print("h=",h)
        c = ZZ(h[2][0]/alpha)
        #pretty_print(h == Nt(0,2)^ZZ(c))
        #pretty_print( n1t^a * n2t^b * n3t^c *== g)
        symbolic_word=[]
        if a>0:
            symbolic_word += a*[4]
        elif a<0:
            symbolic_word += (-a)*[-4]
        if b>0:
            symbolic_word += b*[5]
        elif b<0:
            symbolic_word += (-b)*[-5]
        if c>0:
            symbolic_word += c*[6]
        elif c<0:
            symbolic_word += (-c)*[-6]
        return symbolic_word
    elif g[2,0]==0:
        #print("upper triangular")
        z = g[0][1]/alpha
        a,b = z
        #print("z=",z,"a=",a,"b=",b, (a+omega* b)*alpha)
        h = n2^-b*n1^-a * g
        c = ZZ(h[0][2]/alpha)
        #pretty_print(h==N(0,2)^ZZ(c))
        #pretty_print(n1^a * n2^b * n3^c == g)
        symbolic_word=[]
        if a>0:
            symbolic_word += a*[1]
        elif a<0:
            symbolic_word += (-a)*[-1]
        if b>0:
            symbolic_word += b*[2]
        elif b<0:
            symbolic_word += (-b)*[-2]
        if c>0:
            symbolic_word += c*[3]
        elif c<0:
            symbolic_word += (-c)*[-3]
        return symbolic_word
    else:
        print('Error in decomp2_generator_PSUalpha(g). g is neither upper triangular nor lower triangular.')



# The function decomp2_PSUalpha(g) takes an element g in PSU(2,1)(OK,alpha)
# and returns a word in the Euclidean generators.
# The word is a list of integers in {1,...,6,-1,...,-6}
def decomp2_PSUalpha(g):
    #pretty_print("decomposing ",g)
    if not in_SUalpha(g):
        print('Error in decomp_level_alpha(g). g is not in SU(2,1)(OK,alpha)')
        return None
    else:
        if g[2][0]==0:
            if g==MatOK(1):
                return []
            else:
                return decomp2_generator_PSUalpha(g)
        elif g[2][0].norm() > g[0][0].norm():
            #print("top-left corner is smaller")
            quo = g[1][0] / g[0][0]
            z = NearestOK(quo/alpha)
            #print("quo=",quo,"z=",z, "N(quo-z*alpha)=",(quo-alpha*z).norm())
            temp1 = Nt(z*alpha, z.norm()%2)
            g2 = temp1^-1 * g
            quo2 = g2[2][0]/g2[0][0]
            x = round(QQ((quo2-quo2.conjugate())/(2*alpha)))
            temp2 = temp1*Nt(0,2*x)
            return decomp2_generator_PSUalpha(temp2)+decomp2_PSUalpha(temp2^-1*g)
        else:
            #print("bottom-left corner is smaller")
            quo = g[1][0]/g[2][0]
            z = NearestOK(quo/alpha)
            #print("z =",z)
            #print("N(quo-z*alpha)=",RR((quo-z*alpha).norm()))
            temp1 = N(z.conjugate()*alpha,z.norm()%2)
            g2 = temp1^-1 * g
            #pretty_print("g2=",g2)
            quo2 = g2[0][0]/g2[2][0]
            x = round(QQ((quo2-quo2.conjugate())/(2*alpha)))
            #print("quo2 =",quo2)
            #print("x =",x)
            #print((quo2-x*alpha).norm())
            temp2 = temp1*N(0,2*x)
            return decomp2_generator_PSUalpha(temp2) + decomp2_PSUalpha(temp2^-1 * g)



print('Translating each Reidemeister-Schreier generator into a Euclidean word....',end='')
RSgen2Eword={}
for g in RSgens_PSUalpha:
    #print(g)
    ind = RSgens_PSUalpha.index(g)+1
    Eword =decomp2_PSUalpha(MatOK(g))
    Eword_inv = inverse_word(Eword)
    RSgen2Eword[ind] = Eword
    RSgen2Eword[-ind] = Eword_inv
    #print(Eword_inv)
print('done.')


def RSword2Eword(RSword):
    Eword=[]
    for g in RSword:
        Eword =Eword+ RSgen2Eword[g]
    return Eword



print('Rewriting relations as Euclidean words....',end='')
E_relns=[RSword2Eword(RSword) for RSword in RSsymbolic_relns]
print('done.\n')



E_gens = {
    1:n1,
    2:n2,
    3:n3,
    4:n1t,
    5:n2t,
    6:n3t,
    -1:n1^-1,
    -2:n2^-1,
    -3:n3^-1,
    -4:n1t^-1,
    -5:n2t^-1,
    -6:n3t^-1
}

FP_gens ={
    1:FalPar_P,
    2:FalPar_Q,
    3:FalPar_R,
    -1:FalPar_P^-1,
    -2:FalPar_Q^-1,
    -3:FalPar_R^-1
}
def FPword2elt(word):
    return prod([FP_gens[x] for x in word])
def RSword2elt(word):
    return prod([RSindices_PSUalpha[x] for x in word])
def Eword2elt(word):
    return prod([E_gens[x] for x in word])






print('For each Euclidean generator n, we need to add another relation.\n'
    'This relation consists of first writing n in terms of the Reidemeister Schreier generators\n'
    'and then replacing each RS generator by its decomposition into E-generators.\n'
    'The first step (decomposition into RS generators is achieved by decomposing n\n'
    'into Falbel-Parker generators, and then using the RS-tree.\n') 


print('Writing n1 as a Falbel-Parker word ...',end='')
n1_FPword = [2,-1,2,-1,-1,2,-1,2,-1,-1,-1,2,-1,2,-1,2,2]
print(FPword2elt(n1_FPword)==n1)

print('Rewriting n1 as an RS word ...',end='')
current_point = MatOK(1)
n1_RSword=[]
for x in n1_FPword:
    current_point,skip = Schreier_decomp(current_point * FP_gens[x])
    if skip != 1:
        n1_RSword.append(RSgen2index(skip))
print(n1_RSword)
print('Checking ...',RSword2elt(n1_RSword)==n1, current_point==1)

print('Rewriting n1 as an E word ...',end='')
n1_Eword = RSword2Eword(n1_RSword)
print(Eword2elt(n1_Eword)==n1)

print('Creating and testing extra relation for n1 ...',end='')
extra_reln_n1 = n1_Eword + [-1]
print(Eword2elt(extra_reln_n1)==1, '\n')

E_relns.append(extra_reln_n1)


print('Rewriting n2 as an FP word ...',end='')
n2_FPword =[-1,2,-1,2,1,-2,1,-2]
print(FPword2elt(n2_FPword)==n2)

print('Rewriting n2 as an RS word ...',end='')
current_point = MatOK(1)
n2_RSword=[]
for x in n2_FPword:
    current_point,skip = Schreier_decomp(current_point * FP_gens[x])
    if skip != 1:
        n2_RSword.append(RSgen2index(skip))
print(n2_RSword)
print('Checking ...',RSword2elt(n2_RSword)==n2, current_point==1)

print('Rewriting n2 as an E word ...',end='')
n2_Eword = RSword2Eword(n2_RSword)
print(Eword2elt(n2_Eword)==n2)

print('Creating and testing extra relation for n2 ...',end='')
extra_reln_n2 = n2_Eword + [-2]
print(Eword2elt(extra_reln_n2)==1,'\n')

E_relns.append(extra_reln_n2)


print('Rewriting n3 as an FP word ...',end='')
n3_FPword =[2,2]
print(FPword2elt(n3_FPword)==n3)

print('Rewriting n3 as an RS word ...',end='')
current_point = MatOK(1)
n3_RSword=[]
for x in n3_FPword:
    current_point,skip = Schreier_decomp(current_point * FP_gens[x])
    if skip != 1:
        n3_RSword.append(RSgen2index(skip))
print(n3_RSword)
print('Checking ...', RSword2elt(n3_RSword)==n3, current_point==1)

print('Rewriting n3 as an E word ...',end='')
n3_Eword = RSword2Eword(n3_RSword)
print(Eword2elt(n3_Eword)==n3)

print('Creating and testing extra relation for n3 ...',end='')
extra_reln_n3 = n3_Eword + [-3]
print(Eword2elt(extra_reln_n3)==1,'\n')

E_relns.append(extra_reln_n3)


print('Rewriting n1t as an FP word ...',end='')
n1t_FPword= [3]+inverse_word(n1_FPword)+[2,2,3]
print(FPword2elt(n1t_FPword)==n1t)

print('Rewriting n1t as an RS word ...',end='')
current_point = MatOK(1)
n1t_RSword=[]
for x in n1t_FPword:
    current_point,skip = Schreier_decomp(current_point * FP_gens[x])
    if skip != 1:
        n1t_RSword.append(RSgen2index(skip))
print(n1t_RSword)
print('Checking ...',RSword2elt(n1t_RSword)==n1t, current_point==1)

print('Rewriting n1t as an E word ...',end='')
n1t_Eword = RSword2Eword(n1t_RSword)
print(Eword2elt(n1t_Eword)==n1t)

print('Creating and testing extra relation for n1t ...',end='')
extra_reln_n1t = n1t_Eword + [-4]
print(Eword2elt(extra_reln_n1t)==1,'\n')

E_relns.append(extra_reln_n1t)


print('Rewriting n2t as an FP word ...',end='')
n2t_FPword = [3]+n2_FPword + n1_FPword + n3_FPword +[3]
print(FPword2elt(n2t_FPword)==n2t)

print('Rewriting n2t as an RS word',end='')
current_point = MatOK(1)
n2t_RSword=[]
for x in n2t_FPword:
    current_point,skip = Schreier_decomp(current_point * FP_gens[x])
    if skip != 1:
        n2t_RSword.append(RSgen2index(SUalpha2PSUalpha(skip)))
print(n2t_RSword)
print('Checking ... ',RSword2elt(n2t_RSword)==n2t, current_point==1)

print('Rewriting n2t as an E word ...',end='')
n2t_Eword = RSword2Eword(n2t_RSword)
print(Eword2elt(n2t_Eword)==n2t)

print('Creating and testing extra relation for n2t',end='')
extra_reln_n2t = n2t_Eword + [-5]
print(Eword2elt(extra_reln_n2t)==1,'\n')

E_relns.append(extra_reln_n2t)


print('Rewriting n3t as an FP word ...',end='')
n3t_FPword = [3,2,2,3]
print(FPword2elt(n3t_FPword)==n3t)

print('Rewriting n3t as an RS word',end='')
current_point = MatOK(1)
n3t_RSword=[]
for x in n3t_FPword:
    current_point,skip = Schreier_decomp(current_point * FP_gens[x])
    if skip != 1:
        n3t_RSword.append(RSgen2index(SUalpha2PSUalpha(skip)))
print(n3t_RSword)
print('Checking ...',RSword2elt(n3t_RSword)==n3t, current_point==1)

print('Rewriting n3t as an E word ...',end='')
n3t_Eword = RSword2Eword(n3t_RSword)
print(Eword2elt(n3t_Eword)==n3t)

print('Creating and testing extra relation for n3t ...')
extra_reln_n3t = n3t_Eword + [-6]
print(Eword2elt(extra_reln_n3t)==1,'\n')

E_relns.append(extra_reln_n3t)

print('We now have a presentation for PSU(2,1)(OK,alpha) with the six Euclidean generators\n'
    'n1,n2,n3,n1t,n2t,n3t\n'
    'and',len(E_relns),'relations.\n'
    'Our next step is to simplify this presentation using GAP.')

G_free.<N1,N2,N3,Nt1,Nt2,Nt3> = FreeGroup(6)
abstract_relns = [G_free(Eword) for Eword in E_relns]
Gamma = G_free.quotient(abstract_relns)

GammaGap=Gamma.gap()
print('Number of relations remaining ...', end = '')
print(len(GammaGap.RelatorsOfFpGroup()), end=' ')
for x in range(50):
    GammaGap=GammaGap.SimplifiedFpGroup()
    print(len(GammaGap.RelatorsOfFpGroup()), end=' ')

print('\n\n',GammaGap,'with relations:')
for reln in GammaGap.RelatorsOfFpGroup():
    print('    ',Gamma(reln))
print('Abelian invariants',GammaGap.AbelianInvariants(),'\n')

simplified_E_relns = [G_free(reln).Tietze() for reln in GammaGap.RelatorsOfFpGroup()]

save(simplified_E_relns,'Simplified_Euclidean_Relns.sobj')


# Find a short relation which expresses Nt2 in terms of the other generators.
# This code should produce the relation
#
#      Nt2 = N3^-1*N1*Nt1*N1*N3^-2*N2
#
print('Finding a short relation which allows us to eliminate n2t from the list of Euclidean generators')
for reln in abstract_relns:
    word = reln.Tietze()
    if word.count(5)==1 and word.count(-5)==0 and len(word)<=8:
        print('    ',reln^-1,'= 1')
    if word.count(-5)==1 and word.count(5)==0 and len(word)<=8:
        print('    ',reln,'= 1')
# Check the relation found above
print('Checking the expression for n2t ...',end='')
print(n2t==n3^-1*n1*n1t*n1*n3^-2*n2,'\n')



# This is the weight 1 multiplier system c*tau+d.
# It's assumed here that tau is in the form (a,b,1)^t, where Tr(a)+N(b)<0.
def J1(g,tau):
    return (g*tau)[2]


# action(g,tau) is the action of g on tau.
# This returns the vector in Span(g*tau), normalized so that its final coordinate is 1.
def action(g,tau):
    temp = g*tau
    return temp[2]^-1 *temp


# argj(g,tau) returns the a tau-continuous branch of arg(j(g,tau)).
# If g is in the Borel subgrop then this is just a constant depending on g.
# If not then j(g,tau) is in a*H, where a =i*g[2,0] and H is the upper half-plane.
# to guarantee continuity, we add the constant arg(a) to the continuous function arg(j(g,tau)/a).
# Note that arg is continuous away from the negative reals, an in particular on H.
# It's assumed here that tau is in the form (a,b,1)^t, where Tr(a)+N(b)<0.
# It's also assumes that g is in SU(2,1)(K).
# This is an approximation - it gives a real value.
def argj(g,tau):
    g_matrix = MatOK(g)
    #pretty_print(g_matrix,tau)
    a= alpha*g_matrix[2,0]
    if a==0:
        return arg( K2CC(g_matrix[2,2]))
    else:
        a = alpha*g[2,0]
        return arg(K2CC(a)) + arg(K2CC(J1(g,tau)/a))

# sigma(g,h) is a 2-cocycle on SU(2,1) taking values in ZZ. It corresponds to the univarsal cover of SU(2,1).
# It uses the approximations in argj(g,tau). However, it only uses three values, and sums them.
# The errors would need to sum to at least pi for this to give the wrong answer.
def sigma(g,h):
    tau0 = VK([-1,0,1])
    winding_num = (argj(g*h,tau0) - argj(g,action(h,tau0)) - argj(h,tau0))/(2*pi)
    return(ZZ(round(winding_num)))
    


def E_reln_winding_number(reln):
    current_product=1
    current_winding_number=0
    for x in reln:
        g = E_gens[x]
        current_winding_number += sigma(current_product,g)
        if x<0:
            current_winding_number -= sigma(g,g^-1)
        current_product = current_product*g
    if current_product==MatOK(1):
        #print(current_winding_number)
        return current_winding_number
    else:
        print('Error in RS_reln_winding_num(reln). Your reln is not a relation.')
        print(reln)
        raise ValueError



print('For each relation r, we find the corresponding relation in the '
    'pre-image of our group in the universal cover.\n'
    'This pre-image is generated by lifts of the generators n1,n2,n3,n1t,n3t,'
    'together with a generator of the centre.\n'
    'We then find, for each lifted relation, the row vector in Z^6 of multiplicities of generators.\n'
    'The first five columns correspond to lifts of the Euclidean generators n1,n2,n3,n1t,n3t.\n'
    'The last column corresponds to a generator of the centre of the universal cover.')
row_list=[]
for reln in GammaGap.RelatorsOfFpGroup():
    gt=list(G_free(reln).Tietze())
    gt_correct=[]
    for x in gt:
        if x==5:
            gt_correct.append(6)
        elif x==-5:
            gt_correct.append(-6)
        else:
            gt_correct.append(x)
    gt=gt_correct
    row_list.append( [gt.count(x) -gt.count(-x) for x in [1,2,3,4,6]] + [-E_reln_winding_number(gt)])
mx = matrix(row_list)
print(mx,'\n')

print('A 2-cycle is a formal Z-linear combination of relations,\n'
    'such that the total multiplicity of each generator is 0.\n'
    'Finding a basis for the 2-cocycles (i.e. the kernel of the matrix above) ...')
cycles = mx[0:,:5].kernel().gens()
for Z in cycles:
    print('    ',Z)
    #print(Z*mx)


print('\nIn the lifted relation for any 2-cycle, the total multiplicity of each of the\n'
    'lifted generators is zero. We define the winding number of the cycle to be the\n'
    'multiplicity of the generator of the centre in the lifted relation.\n'
    'The weight denominator is the highest common factor of these winding numbers.\n')
for Z in cycles:
    #print('\n',Z)
    cycle=G_free(1)
    T_cycle = []
    for x in range(13):
        reln = G_free(GammaGap.RelatorsOfFpGroup()[x])
        T_reln = reln.Tietze()
        power = Z[x]
        if power>0:
            T_cycle += power* T_reln
        elif power<0:
            T_cycle +=(-power)* inverse_word(T_reln)
    # need to correct T_cycle by changing 5 to 6 and -5 to -6 ...
    T_cycle_correct = []
    for letter in T_cycle:
        if letter == 5:
            T_cycle_correct.append(6)
        elif letter == -5:
            T_cycle_correct.append(-6)
        else:
            T_cycle_correct.append(letter)
    #print(T_cycle_correct)
    #print(G_free(T_cycle_correct))
    winding_number = (Z*mx[0:,-1:])[0]
    #print([T_cycle.count(x) -T_cycle.count(-x) for x in [1,2,3,4,5]] + [E_reln_winding_number(T_cycle)])
    if winding_number!=0:
        print('The following cycle has winding number ',winding_number)
        pretty_print(G_free(T_cycle_correct))

print('\nTherefore the weight denominator of the group PSU(2,1)(OK,alpha) is 1.\n')

