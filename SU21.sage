run_tests=True

print('Defining the following for SU(2,1) over the Eisenstein integers:\n')

K.<omega>=NumberField(x^2+x+1)
OK = K.ring_of_integers()
omega = OK(omega)     # omega is a primitive cube root of unity
alpha = 2*omega+1     # alpha is sqrt(-3)
mu_3 = [1,omega,omega^2]
mu_6 = [(1+omega)^n for n in range(6)]
F3.<one_mod,omega_mod> = OK.quotient_ring(OK.ideal([alpha])) # the finite field OK/(alpha) with 3 elements. 
R9.<one_mod,omega_mod> = OK.quotient_ring(OK.ideal([3])) # the finite ring OK/(3) with 9 elements. 
def K2CC(zK):         # convert elements of K into complex numbers.
    a,b = list(K(zK))
    return CC(a+b*exp(2*pi*i/3))
"""
return the nearest Eisenstein integer to q.
The difference z-q will be in a regular hexagon centred at 0
with two vertical sides distance 1 apart. This has corners at zeta/alpha,
where zeta runs through the 6th roots of unity.
"""
def NearestOK(q):
    y=QQ((q-K(q).conjugate())/(alpha))
    y_int = round(y)
    temp = q - y_int * omega
    x_int = round(temp.trace()/2)
    z = x_int+y_int*omega
    for test in [z+omega,z+omega+1, z-omega,z-omega-1]:
        if (q-test).norm() < (q-z).norm():
            return test
    return z


print('omega             a primitive cube root of unity \n'
	  'alpha=2*omega+1   a square root of -3 \n'
	  'K=Q(omega),\n'
      'OK=Z[omega] \n'
	  'mu_3 and mu_6     the groups of cube and 6th roots of unity in K \n'
	  'F3 and R9         the quotient rings OK/alpha and OK/3 \n'
      'K2CC(zK)          an ambedding of K into the complex numbers.\n'
      'NearestOK(q)      the nearest element of OK to an element q in K.')


MatOK = MatrixSpace(OK,3,3)
VOK = OK^3
VK = K^3
MatF3 = MatrixSpace(F3,3,3)
MatR9 = MatrixSpace(R9,3,3)
print('MatOK             the space of 3x3 matrices over OK.'
	  'MatF3, MatR9      the spaces of 3x3 matrices over F3 and R9')


J = MatOK([[0,0,1],[0,1,0],[1,0,0]])
print('J                 the antidiagonal matrix of 1s')

def Form(v1,v2):      # the Hermitian form corresponding to J.
    return (v1.transpose()*J*v2.conjugate())[0,0]
print('Form(v1,v2)       The Hermitian form on OK^3 corresponding to J.')

def N(z,x):     # construct an upper triangular element in SU(2,1)(OK)
    if z in OK and x in ZZ and (OK(z).norm()+x)%2==0:
        return MatOK([[1,z,(-OK(z).norm()+x*alpha)/2],[0,1,-OK(z).conjugate()],[0,0,1]])
    else:
        print("Error in N(z,x)")
        print("You entered z=',z','and x=",x)
        return None
print('N(z,x)            an upper triangular unipotent matrix in SU(2,1).')

def Nt(z,x):     # construct a lower triangular element in SU(2,1)(OK)
    if z in OK and x in ZZ and (OK(z).norm()+x)%2==0:
        return N(z,x).transpose()
    else:
        print("Error in Nt(z,x)")
        print("You entered z=',z','and x=",x)
        return None
print('Nt(z,x)           the transpose of N(z,x)')

def in_U(g):    # check whether a matrix is in U(2,1)(OK)
    if g in MatOK:
        return MatOK(g).transpose()*J*MatOK(g).conjugate()==J
    else:
        print("Error in function in_SU.")
        print("The argument should be a matrix")
        print("Your argument is ",g," of type ", g.parent())
        return None
print('in_U(g)           check whether g is in the unitary group U(2,1)(OK)')

def in_SU(g):    # check whether a matrix is in SU(2,1)(OK)
    if g in MatOK:
        return in_U(g) and MatOK(g).det()==1
    else:
        print("Error in function in_SU(g).")
        print("The argument should be a matrix")
        print("Your argument is ",g," of type ", g.parent())
        return None
print('in_SU(g)          check whether g is in the special unitary group SU(2,1)(OK)')

def in_SUalpha(g):
    return in_SU(g) and MatF3(g)==MatF3(1)
print('in_SUalpha(g)     check whether g is in the level alpha principal congruence subgroup SU(2,1)(OK,3)')

def in_SU3(g):  # check whether g is in the subgroup SU(2,1)(OK,3) of U(2,1)(OK)
    return in_SU(g) and MatR9(g)==MatR9(1)
print('in_SU3(g)         check whether g is in the level 3 principal congruence subgroup')


# check whether g is in the subgroup Z.SU(2,1)(OK,3) of SU(2,1)(OK), where Z is the centre
# (a cyclic group of order 3). Note that coset reps for Z.SU(2,1)(OK,3) in U(2,1)(OK,3) are
# the same as those of SU(2,1)(OK,3) in PSU(2,1)(OK,3).
def in_ZSU3(g):
    return in_SU(g) and (MatR9(g)==MatR9(1) or MatR9(g)==MatR9(omaga) or MatR9(g)==MatR9(omaga^2))
print('in_ZSU3(g)        check whether g is in Z.SU(2,1)(OK,3), where Z is the centre of SU(2,1)')


"""
generate random elements of PSU(2,1)(OK,alpha) for testing purposes.
"""
def random_PSUalpha():
    temp = MatOK(1)
    for counter in range(10):
        z1 = 3*randint(-10,10) + randint(-10,10)*alpha
        x1 = 6*randint(-10,10) + (z1.norm()%2)
        z2 = 3*randint(-10,10) + randint(-10,10)*alpha
        x2 = 6*randint(-10,10) + (z2.norm()%2)
        temp *= N(z1,x1)*Nt(z2,x2)
    return temp
print('random_PSUalpha() a random element of PSU(2,1)(OK,alpha).')
    
# The following is a set of generators for SU(2,1)(OK).    
#gamma1 = N(1,1)
#gamma2 = N(omega,1)
#gamma3 = -J
#generatorsSU1 = [gamma1,gamma2,gamma3,Z]

# The following is a set of generators of PSU(2,1)(OK,alpha).
# we shall always identify this group with a subgroup of SU(2,1)(OK,alpha), as
# SU(2,1)(OK,alpha) is a direct sum of this group with its centre <Z>
n1=N(alpha,1)
n2=N(omega*alpha,1)
n3=N(0,2)
n1t = n1.transpose()
n2t = n2.transpose()
n3t = n3.transpose()
print('n1,n2,n3,n1t,n3t  generators of SU(2,1)(OK,alpha).')

#generatorsSUalpha = [n1,n2,n3,n1t,n2t,n3t,Z]
#generatorsPSUalpha = [n1,n2,n3,n1t,n3t]
gensPSUalpha = {
    1:n1,
    2:n2,
    3:n3,
    4:n1t,
    5:n3t,
    -1:n1^-1,
    -2:n2^-1,
    -3:n3^-1,
    -4:n1t^-1,
    -5:n3t^-1
}
labelsPSUalpha = { str(g): g_key for g_key,g in gensPSUalpha.items()}
print('gensPSUalpha      a dictionary to these generators and their inverses \n'
	  '                  with keys 1,...,5, -1,...,-5')


def gen2Tz(g,gen_dict = gensPSUalpha):
    for key,h in gen_dict.items():
        if MatOK(g)==MatOK(h):
            return key
    print('Error in gen2Tz(g).\n'
          'g should be a generator or the inverse of a generator.\n'
          'you have entered g=',g)
    return None
def word2Tzword(word,gen_dict = gensPSUalpha):
    return [gen2Tz(g,gen_dict = gen_dict) for g in word]


relationsPSUalpha =[
    [n2^-1,n3^-1,n2,n3],
    [n1^-1,n3^-1,n1,n3],
    [n1t^-1,n3t^-1,n1t,n3t],
    3*[n3,n3t],
    [n3,n2,n1,n2^-1,n3,n1^-1,n3],
    3*[n1^-1,n3,n1t^-1],
    [n3t^-1,n2,n3t,n1t^-1,n1^-1,n2^-1,n3,n1t,n3^-1,n1],
    [n1t^-1,n1^-1,n3,n3t,n2,n1,n3t^-1,n1t,n2^-1,n3^-1],
    [n3t^-1,n1t,n1,n3t,n3^-1,n1^-1,n2^-1,n1t^-1,n3^-1,n3^-1,n2],
    [n3t^-1,n2,n1,n3t^-1,n1t,n1,n1t,n1,n3t^-1,n1t,n2^-1],
    [n3,n3t,n1,n1t,n3t^-1,n2^-1,n1t,n3^-1,n1,n1t,n2,n1],
    [n3^-1,n1,n1t,n2,n3,n1,n3t^-1,n1^-1,n1t^-1,n3t,n1^-1,n2^-1],
    [n1t^-1,n3^-1,n3t,n3,n1^-1,n1t^-1,n2,n1,n3^-1,n1t,n1,n3t^-1,n1t,n2^-1,n1^-1,n3]
]


"""relationsPSUalpha =[
    [n2^-1,n3^-1,n2,n3],
    [n1^-1,n3^-1,n1,n3],
    [n1t^-1,n3t^-1,n1t,n3t],
    3*[n3,n3t],
    [n3,n2^-1,n3,n1^-1,n2,n3,n1],
    [n1t^-1,n1^-1,n3]+2*[n1t^-1,n3,n1^-1],
    [n1t^-1]+2*[n1^-1,n1t^-1,n3t]+[n1^-1,n3t],
    [n2^-1,n3t^-1,n1,n3^-1,n1t,n3,n2,n3t,n1t^-1,n1^-1],
    [n3t,n1,n1t,n3t^-1,n2^-1,n3,n1^-1,n1t^-1,n2,n3^-1],
    [n1^-1,n3,n1t^-1,n3^-1,n2,n1,n1t,n3t^-1,n2^-1,n3t],
    [n2^-1,n1t^-1,n3t^-1,n2,n1t^-1,n3t,n1^-1,n3t^-1,n1^-1,n1t^-1,n3t,n1^-1]
]
"""
Tz_relns_level_alpha = [word2Tzword(word) for word in relationsPSUalpha]



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

#pretty_print(FalPar_generators)
#for g in FalPar_generators:
#    print(in_U(FalPar_generators[g]))
#for r in FalPar_relations:
#    print(prod(r))


"""
The following functions convert between matrices in PSU(2,1)(alpha) and words with respect to
the generating set reduced_gens.

This is the set of Euclidean generators with n2t removed,
since it may be expressed in terms of the other generators.
"""
def word2elt(word,
             generator_dict = gensPSUalpha):
    return prod([generator_dict[letter] for letter in word])


def upper_tri2word_PSUalpha(n_matrix):
    z = n_matrix[0,1]
    a,b = z/alpha
    t = (n_matrix[0,2]/alpha).trace()
    c = (t - 3*a*b-a-b)/2
    #pretty_print(n_matrix, N(z,t), n1^a*n2^b*n3^c)
    if a>0:
        word = a*[1]
    else:
        word = (-a)*[-1]
    if b>0:
        word += b*[2]
    else:
        word += (-b)*[-2]
    if c>0:
        word += c*[3]
    else:
        word += (-c)*[-3]
    return word
    
def lower_tri2word_PSUalpha(nt_matrix):
    n2t_word = [-3,1,4,1,-3,-3,2]
    n2t_inv_word = [-2,3,3,-1,-4,-1,3]
    #print(word2PSUalpha(n2t_word) == n2t)
    #print(word2PSUalpha(n2t_inv_word) == n2t^-1)
    z = nt_matrix[1,0]
    t = (nt_matrix[2,0]/alpha).trace()
    a,b = z/alpha
    c = (t - 3*a*b-a-b)/2
    #pretty_print(nt_matrix, Nt(z,t), n2t^b*n1t^a*n3t^c)
    if b>0:
        word = b*n2t_word
    else:
        word = (-b)*n2t_inv_word
    if a>0:
        word += a*[4]
    else:
        word += (-a)*[-4]
    if c>0:
        word += c*[5]
    else:
        word += (-c)*[-5]
    return word

def elt2word_PSUalpha(g):
    #pretty_print("decomposing ",g)
    if not in_SUalpha(g):
        print('Error in elt2word_PSUalpha(g). g is not in SU(2,1)(OK,alpha)')
        pretty_print('g=',g)
        return None
    else:
        if g[2][0]==0:
            return upper_tri2word_PSUalpha(g)
        elif g[2][0].norm() > g[0][0].norm():
            quo = g[1][0] / g[0][0]
            z = NearestOK(quo/alpha)
            temp1 = Nt(z*alpha, z.norm()%2)
            g2 = temp1^-1 * g
            quo2 = g2[2][0]/g2[0][0]
            x = round(QQ((quo2-quo2.conjugate())/(2*alpha)))
            temp2 = temp1*Nt(0,2*x)
            return lower_tri2word_PSUalpha(temp2)+elt2word_PSUalpha(temp2^-1*g)
        else:
            quo = g[1][0]/g[2][0]
            z = NearestOK(quo/alpha)
            temp1 = N(z.conjugate()*alpha,z.norm()%2)
            g2 = temp1^-1 * g
            quo2 = g2[0][0]/g2[2][0]
            x = round(QQ((quo2-quo2.conjugate())/(2*alpha)))
            temp2 = temp1*N(0,2*x)
            return upper_tri2word_PSUalpha(temp2) + elt2word_PSUalpha(temp2^-1 * g)


        
if run_tests:
    print('Testing the functions converting between elements of PSU(2,1)(alpha) and Tietze words.\n')
    for counter in range(10):
        r = randint(-5,5)
        s = randint(-5,5)
        t = 2*randint(-5,5) + ((r+r*s+s)%2)
        z = (r+s*omega)*alpha
        test_mat = N(z,t)
        print(test_mat == word2elt(upper_tri2word_PSUalpha(test_mat)),end='')
        print(test_mat.transpose()==word2elt(lower_tri2word_PSUalpha(test_mat.transpose())),end='')
    for counter in range(10):
        test_mat = random_PSUalpha()
        print(word2elt(elt2word_PSUalpha(test_mat))==test_mat, end='')
else:
    print('NOT RUNNING TESTS.')


"""def reduce_keyword(keyword):
    reduced_word = keyword
    unfinished = True
    while unfinished:
        unfinished=False
        for x in range(len(reduced_word)-1):
            if reduced_word[x] == -reduced_word[x+1]:
                reduced_word = reduced_word[:x]+ reduced_word[x+2:]
                unfinished = True
                break
    return reduced_word
"""


print('Checking relations in PSU(2,1)(OK,alpha)')
for reln in relationsPSUalpha:
    print(prod([g for g in reln])==1,end=' ')
for reln in Tz_relns_level_alpha:
    print(prod([gensPSUalpha[label] for label in reln])==1, end=' ')
print('')

# This is the weight 1 multiplier system c*tau+d.
# It's assumed here that tau is in the form (a,b,1)^t, where Tr(a)+N(b)<0.
def J1(g,tau):
    return (g*tau)[2]
print('J1(g,tau)        The weight 1 multiplier system.\n'
     '                 where g is in SU(2,1)(OK) and tau=(a,b,1)^t has \n'
     '                 negative Hermitian norm.')

# action(g,tau) is the action of g on tau.
# This returns the vector in Span(g*tau), normalized so that its final coordinate is 1.
def action(g,tau):
    temp = g*tau
    return temp[2]^-1 *temp
print('action(g,tau)    The linear fractional action of g on tau.')

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
print('argj(g,tau)      a branch of the argument of J1(g,tau),\n'
      '                 which is continuous as a function of tau.')
# sigma(g,h) is a 2-cocycle on SU(2,1) taking values in ZZ. It corresponds to the univarsal cover of SU(2,1).
# It uses the approximations in argj(g,tau). However, it only uses three values, and sums them.
# The errors would need to sum to at least pi for this to give the wrong answer.
def sigma(g,h):
    tau0 = VK([-1,0,1])
    winding_num = (argj(g*h,tau0) - argj(g,action(h,tau0)) - argj(h,tau0))/(2*pi)
    return(ZZ(round(winding_num)))
print('sigma(g,h)       a Z-valued homogeneous 2-cocycle corresponding to\n'
      '                 the universal cover of SU(2,1).')




"""
 This is a homomorphism phi from PSU(2,1)(OK,alpha) to OK.
 The image is (alpha).
 Pre-images of subgroups of (alpha) will be arithmetic subgroups, some of them non-congruence subgroups.
 
 The phi is defined on generators. It takes Nt(z,x) to z and N(z,x) to z.conjugate().
 Here z and x must be chosen so that these elements are in the level alpha subgroup,
 ie. z is a multiple of alpha and x is Nz + a multiple of 6.
"""

def phi(g):
    #pretty_print("decomposing ",g)
    if not in_SUalpha(g):
        print('Error in phi(g). g is not in SU(2,1)(OK,alpha)')
        return None
    else:
        if g[2][0]==0:
            return -(OK(g[0][1])/alpha).conjugate()
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
            return z+phi(temp2^-1*g)
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
            return -z + phi(temp2^-1 * g)
print('phi(g)           a homomorphism from SU(2,1)(OK,alpha) to OK')

if run_tests:
    print('\ntesting to confirm that phi is a homomorphism:')
    for counter in range(10):
        g1 = random_PSUalpha()
        g2 = random_PSUalpha()
        print(phi(g1*g2)== phi(g1)+ phi(g2), end = ' ')
        print(phi(g1^-1)==-phi(g1), end=' ')
    print()
else:
    print('NOT RUNNING TESTS.')


"""
A class whose objects are elements of the universal cover of SU(2,1)(RR).
"""
print('Defining a class covering_group_element, wose objects '
      'are elements of the universal cover of SU(2,1)(RR).\n')
class covering_group_element:
    def __init__(self,mat,central_value=0):
        self._mat = MatOK(mat)
        self._central_value = ZZ(central_value)
    def __str__(self):
        return str((self._mat,self._central_value))
    def _latex_(self):
        return latex((self._mat,self._central_value))
    def __repr__(self):
        return str('The element '+str((self._mat,self._central_value))+' in the universal cover of SU(2,1)(RR).')
    # the image of self in SU(2,1)(RR)
    def __eq__(self,other):
        return self._mat==other._mat and self._central_value==other._central_value
    def projection(self):
        return self._mat
    # multiply two elements of the covering group
    def __mul__(self,other):
        product_mat = self._mat * other._mat
        central_value = self._central_value + other._central_value + sigma(self._mat,other._mat)
        return covering_group_element(product_mat,central_value)
    def inverse(self):
        inverse_mat = self._mat^-1
        central_value = -self._central_value - sigma(self._mat,inverse_mat)
        return covering_group_element(inverse_mat,central_value)
    # overloading the power operator
    def __pow__(self,other):
        if other not in ZZ:
            print('Error in taking power of covering group element. The power must be an integer.')
            raise ValueError
            return None
        if other > 0:
            return prod([self for counter in range(other)])
        if other < 0:
            inv = self.inverse()
            return prod([inv for counter in range(-other)])
        return covering_group_element(1)
        
if run_tests:
    print('Testing multiplication and power laws int the covering group.')
    for counter in range(10):
        g1= covering_group_element(random_PSUalpha(), randint(-10,10))
        g2= covering_group_element(random_PSUalpha(), randint(-10,10))
        g3= covering_group_element(random_PSUalpha(), randint(-10,10))
        a = randint(-10,10)
        b = randint(-10,10)
        print((g1*g2)*g3 == (g1*(g2*g3)),end='')
        print(g1^a * g1^b == g1^(a+b),end='')
        print(g1.projection()*g2.projection()==(g1*g2).projection(),end='')
else:
    print('Not testing.')


"""
A class whose objects are arithmetic subgroups of PSU(2,1)(alpha),
regarded as a subgroup of SU(2,1)(RR).
"""
class ArithGp:
    def __init__(self,gp_name,coset_key,latex_name=None):
        self.filename ='SU(2,1)(OK,alpha) subgroup data/'+gp_name+'.sobj'
        try:
            data = load(self.filename)
            self._name = gp_name
            self._coset_key = coset_key
            self._transversal = data[0]
            self._RSgens = data[1]
            self._RSedges = data[2]
            self._index = data[3]
            self._RSrelations = data[4]
            self._no_of_gens = data[5]
            self._no_of_relns = data[6]
            self._reln_matrix = data[7]
            self._extended_reln_matrix = data[8]
            self._gap_group = None
            self._gap_gens = [9]
            self._gap_relns = [10]
            self._no_of_gap_gens = [11]
            self._no_of_gap_relns = [12]
        except FileNotFoundError:
            self._name = gp_name
            self._coset_key = coset_key
            self._transversal = None
            self._RSgens = None
            self._RSedges = None
            self._index = None
            self._RSrelations = None
            self._no_of_gens = None
            self._no_of_relns = None
            self._reln_matrix = None
            self._extended_reln_matrix = None
            self._gap_group = None
            self._gap_gens = None
            self._gap_relns = None
            self._no_of_gap_gens = None
            self._no_of_gap_relns = None
        self._H_1 = None
        self._reln_matrix_smith = None
        self._H_1_elementary_divisors = None
        self._H_1_met = None
        self._extended_reln_matrix_smith = None
        self._H_1_met_elementary_divisors = None
        if latex_name==None:
            self._latex_name = LatexExpr(self._name)
        else:
            self._latex_name = latex_name
        self.__save()
    
    def __str__(self):
        return self._name
    def __repr__(self):
        return 'The arithmetic subgroup '+self._name+' of SU(2,1)(alpha)'
    def _latex_(self):
        return self._latex_name
    def __contains__(self,elt):
        try:
            return self._coset_key(elt) == self._coset_key(MatOK(1))
        except:
            return False
    def __le__(self,other):    # self is a subgroup of other
        gens = self.generator_list()
        for elt in gens:
            if not elt in other:
                return False
        return True    
    def __ge__(self,other):     # self is a supergroup of other
        return other <= self
    def __eq__(self,other):     # self is the same group as other (although they may have quite different internal data)
        return self <= other and other <= self
    def __ne__(self,other):
        return not(self==other)
    def __lt__(self,other):
        return self <= other and not self >= other
    def __gt__(self,other):
        return self >= other and not self <= other
    def __and__(self,other):   # The intersection of two subgroups of PSU(2,1)(OK,alpha) "G & H"
        new_name = self._name + "\\cap " + other._name
        new_latex_name = self._latex_name + LatexExpr('\\cap ')+other._latex_name
        def new_coset_key(elt):
            return str([self._coset_key(elt),other._coset_key(elt)])
        return ArithGp(new_name,new_coset_key,latex_name=new_latex_name)
    
    def normalized_by(self,g):
        gens = self.generator_list()
        for gen in gens:
            test = g*gen*g^-1
            if SUalpha2PSUalpha(test)!=test:
                return False
            if not g*gen*g^-1 in self:
                return False
        return True
    
    """ To do:
    add an attribute with a list of known supergroups, to be checked and updated by __le__
    """
        
    def __save(self):
        data = [self._transversal,
                self._RSgens,
                self._RSedges,
                self._index,
                self._RSrelations,
                self._no_of_gens,
                self._no_of_relns,
                self._reln_matrix,
                self._extended_reln_matrix,
                self._gap_gens,
                self._gap_relns,
                self._no_of_gap_gens,
                self._no_of_gap_relns,
                self._reln_matrix_smith,
                self._H_1,
                self._H_1_elementary_divisors
               ]
        save(data,self.filename)

    
    def SchreierTransversal(self,
                            generator_dict = gensPSUalpha,
                            verbose = False,
                            check = False
                           ):
        if self._transversal != None:
            if verbose:
                print('Already have the Schreier Transversal.')
            return self._transversal
        gp_name = self._name
        coset_key = self._coset_key
        if verbose:
            print('Creating Reidemeister - Schreier graph for '+gp_name+' in PSU(2,1)(alpha) ....')
        Transversal = {coset_key(MatOK(1)):MatOK(1)}
        todolist = [MatOK(1)]  # a list of coset reps not yet operated on by the generators.
        RS_gens = {0:MatOK(1)} # a dictionary to Reidemeister-Schreier generators
        RS_gen_keys = {str(MatOK(1)):0}
        next_gen_label = 1     # the label of the next RS generator to be found
        RS_edges = {}          # the edges of the RS graph in the form rep_key,g_label: (new_rep_key,RS_gen_label)
        while todolist != []:
            rep = todolist.pop(0)
            rep_key = coset_key(rep)
            for g_label,g in generator_dict.items():
                test = rep*g
                test_key = coset_key(test)
                try:
                    test_rep = Transversal[test_key]
                    RS_gen = test * test_rep^-1
                    RS_gen_key = str(RS_gen)
                    try:
                        RS_gen_label = RS_gen_keys[RS_gen_key]
                        RS_edges[rep_key,g_label] = (test_key,RS_gen_label)
                    except KeyError:  # in this case RS_gen is a new generator.
                        RS_gen_inv = RS_gen^-1
                        RS_gens[next_gen_label] = RS_gen
                        RS_gens[-next_gen_label] = RS_gen_inv
                        RS_gen_keys[RS_gen_key] = next_gen_label
                        RS_gen_keys[str(RS_gen_inv)] = -next_gen_label
                        RS_edges[rep_key,g_label] = (test_key,next_gen_label)
                        next_gen_label += 1                   
                except KeyError:
                    Transversal[test_key]=test
                    todolist.append(test)
                    RS_edges[rep_key,g_label] = (test_key,0)
                    if verbose and len(Transversal)%100==0:
                        print(len(Transversal), end='\r')
        self._transversal = Transversal
        self._RSgens = RS_gens
        self._RSedges = RS_edges
        self._index = len(Transversal)
        self._no_of_gens = next_gen_label-1
        if verbose:
            print('Found Schreier transversal of size',len(Transversal),'=',factor(len(Transversal)))
            print('Number of generators:',self._no_of_gens)
        if check:
            print('Testing the Reidemeister-Schreier graph ...', end='')
            for rep_key, rep in Transversal.items():
                for g_label,g in generator_dict.items():
                    new_rep_key , skip_label = RS_edges[rep_key,g_label]
                    new_rep = Transversal[new_rep_key]
                    skip = RS_gens[skip_label]
                    if rep * g != skip * new_rep and coset_key(skip)==coset_key(MatOK(1)):
                        print('TEST FAILED!')
                        return None
            print('all ok.')
        self.__save()
        return self._transversal

    def generators(self, verbose=False):
        if self._RSgens==None:
            if verbose:
                print('Getting Schreier transversal...')
            self.SchreierTransversal(verbose=verbose)
        return self._RSgens
    def generator_list(self, verbose=False):
        self.generators(verbose=verbose)
        return [self._RSgens[counter] for counter in range(1,self._no_of_gens+1)]
    
    def no_of_generators(self,verbose=False):
        if self._no_of_gens==None:
            if verbose:
                print('Getting Schreier transversal...')
            self.SchreierTransversal(verbose=verbose)
        return self._no_of_gens
    
    def index(self, verbose=False):
        if self._index == None:
            if verbose:
                print('getting Schreier transversal...')
            self.SchreierTransversal(verbose=verbose)
        return self._index

    # add a letter to a word in Tietze format.
    # this automatically cancels if n is added to a word ending in -n.
    def word_add(word,letter):
        if letter !=0:
            try:
                if letter == - word[-1]:
                    word.pop()
                else:
                    word.append(letter)
            except IndexError:
                word.append(letter)
        return None

    # find a single Reidemeister-Schreier relation in a subgroup, given
    # a relation in the group,
    # a starting point in the transversal,
    # and the Reidemeister-Schreier graph.
    # This defaults to subgroups pf PSU(2,1)(OK,alpha)
    def RS_relation(self,
                    reln,    # a relation in the larger group
                    rep_key, # a transversal key
                    generator_dict = gensPSUalpha):
        current_point = rep_key
        RS_reln = []
        for g in reln:
            current_point, skip = self._RSedges[current_point,g]
            ArithGp.word_add(RS_reln,skip)
        if current_point==rep_key:
            return RS_reln
        else:
            print('Error in RS_relation. Your word does not generate an element of the subgroup.')
            print('relation =',reln)
            print('starting point =',rep_key)
            print('finished at',current_point)
            return None

    # find all the Reidemeister-Schreier relations in a subgroup
    def RS_relations(self,
                     reln_list = Tz_relns_level_alpha,
                     verbose = False):
        if self._RSrelations != None:
            if verbose:
                print('Already have this.')
            return self._RSrelations
        if self._transversal == None:
            self.SchreierTransversal(verbose=verbose)
        if verbose:
            print('Finding Reidemeister-Schreier relations....')
        self._RSrelations = []
        for reln in reln_list:
            for rep_key in self._transversal.keys():
                RS_reln = self.RS_relation(reln,rep_key,gensPSUalpha)
                if RS_reln !=[]:
                    self._RSrelations.append(RS_reln)
        self._no_of_relns = len(self._RSrelations)
        if verbose:
            print('Number of relations:',self._no_of_relns)
        self.__save()
        return self._RSrelations



    # convert relations to a matrix, whose row reduction gives the abelianization

    def reln_matrix(self,verbose=False):
        if self._reln_matrix != None:
            if verbose:
                print('Already have this.')
            return self._reln_matrix
        if self._RSrelations == None:
            self.RS_relations(verbose=verbose)
        no_of_gens = self._no_of_gens
        no_of_relns = self._no_of_relns
        if verbose:
            print('converting list a of',no_of_relns,'relations in',no_of_gens,'generators to a matrix ...')
        self._reln_matrix = zero_matrix(ZZ,
                                       no_of_relns,
                                       no_of_gens,
                                       sparse=True)
        for current_row in range(no_of_relns):
            if verbose:
                print(current_row,end='\r')
            reln = self._RSrelations[current_row]
            for g_label in reln:
                if g_label > 0:
                    self._reln_matrix[current_row,g_label-1]+=1
                else:
                    self._reln_matrix[current_row,-g_label-1]-=1
        if verbose:
            print('done.')
        self.__save()
        return self._reln_matrix

    # convert relations to a matrix, whose row reduction gives the abelianization
    # of the pre-image in a covering group defined by a cocycle.
    def extended_reln_matrix(self, verbose = False, cocycle = sigma):
        if self._extended_reln_matrix !=None:
            if verbose:
                print('Already have this.')
            return self._extended_reln_matrix
        if self._RSrelations == None:
            self.RS_relations(verbose=verbose)
        no_of_gens = self._no_of_gens
        no_of_relns = self._no_of_relns
        if verbose:
            print('converting a list of',no_of_relns,'relations in',no_of_gens,'generators to a matrix ...')
        extended_matrix = zero_matrix(ZZ,no_of_relns,no_of_gens+1,sparse=True)
        for current_row in range(no_of_relns):
            if verbose:
                print(current_row,end='\r')
            reln = self._RSrelations[current_row]
            current_posn = MatOK(1)
            for g_label in reln:
                g = self._RSgens[g_label]
                if g_label > 0:
                    extended_matrix[current_row,g_label-1]+=1
                else:
                    extended_matrix[current_row,-g_label-1]-=1
                    extended_matrix[current_row,-1] -=cocycle(g^-1,g)
                extended_matrix[current_row,-1] += cocycle(current_posn,g)
                current_posn = current_posn * g
            if current_posn !=1:
                print('Error in reln2extended_row. Not a relation.')
                return None
        self._extended_reln_matrix=extended_matrix
        if self._reln_matrix==None:
            self._reln_matrix= extended_matrix[:,:-1]
        if verbose:
            print('done.')
        self.__save()
        return self._extended_reln_matrix
    
    def gap_group(self,verbose=False):
        if self._gap_group != None:
            return self._gap_group
        if self._RSrelations == None:
            self.RS_relations(verbose=verbose)
        Gamma = FreeGroup(self._no_of_gens) / self._RSrelations
        GammaGap = Gamma.gap()
        presentation = GammaGap.PresentationFpGroup( 0 )
        presentation.TzInitGeneratorImages()
        for counter in range(10):
            presentation.SimplifyPresentation()
        NewGammaGap = presentation.FpGroupPresentation()
        gap_gens ={0:MatOK(1)}
        self._no_of_gap_gens =0
        for g in presentation.TzPreImagesNewGens():
            self._no_of_gap_gens += 1
            g_matrix = word2elt(Gamma(g).Tietze(),self._RSgens)
            gap_gens[self._no_of_gap_gens] = g_matrix
            gap_gens[-self._no_of_gap_gens] = g_matrix^-1
        self._gap_gens = gap_gens
        self._gap_relns = [Gamma(reln).Tietze() for reln in NewGammaGap.RelatorsOfFpGroup()]
        self._gap_group = NewGammaGap
        self._no_of_gap_relns = len(NewGammaGap.RelatorsOfFpGroup())
        self.__save()
        return self._gap_group
    
    """
    find a word in the RS generators, which evaluates to the element g.
    Such a word is not unique. The algorithm first uses the Euclidean algorithm to
    write g as a word in the Euclidean generators of PSU(2,1)(alpha). It then walks around
    the Reidemeister-Schreier graph of the quotient PSU(2,1)(alpha)/self using this word.
    """
    def elt2RSword(self,g):
        sualpha_word = elt2word_PSUalpha(g)
        current_point = self._coset_key(MatOK(1))
        RS_word = []
        for letter in sualpha_word:
            current_point, skip = self._RSedges[current_point,letter]
            ArithGp.word_add(RS_word,skip)
        if current_point==self._coset_key(MatOK(1)):
            return RS_word
    
    """
    convert a word in the Reidemeister-Schreier generators to
    an elenment of the group.
    """
    def RSword2elt(self,RS_word):
        return word2elt(RS_word,generator_dict=self.generators())

    """
    given a word in RS generators, this evaluates the product of the lifts of the corresponding
    generators, giving an element of the pre-image of self in the univarsal cover of SU(2,1)(RR).
    """
    def RSword2met(self,RS_word):
        return prod([covering_group_element(self._RSgens[letter]) for letter in RS_word])
    """
    Old version
    
    def RSword2met(self,RS_word):
        self.generators()
        central_value=0
        current_elt = MatOK(1)
        for letter in RS_word:
            gen = self._RSgens[letter]
            central_value += sigma(current_elt,gen)
            current_elt = current_elt * gen
        return (current_elt,central_value)

    """
    
    """
    Given an element (elt,central_value), this returns a pair (RS_word,new_central_value),
    which evaluates to (elt,central_value).
    """
    def elt2RSword_met(self,elt,central_value=0):
        RS_word = self.elt2RSword(elt)
        new_central_value = central_value - self.RSword2met(RS_word)[1]
        return (RS_word,new_central_value)

    
    
    def reduced_reln_matrix(self,verbose=False):
        if verbose:
            print('Getting the relation matrix.')
        mat = self.reln_matrix(verbose=verbose)
        if mat.is_mutable():
            if verbose:
                print('Row reducing the matrix. This could take some time, or might not finish if the matrix is big.')
            mat.echelonize()
            mat=mat[:mat.rank()]
            self._reln_matrix = mat
            self.__save()
        elif verbose:
            print('Relation matrix is immutable, so it\'s probably already reduced. Doing no more reduction.')
        return self._reln_matrix

    def reduced_extended_reln_matrix(self,verbose=False):
        if verbose:
            print('Getting the extended relation matrix.')
        mat = self.extended_reln_matrix(verbose=verbose)
        if mat.is_mutable():
            if verbose:
                print('Row reducing the matrix. This could take some time, or might not finish if the matrix is big.')
            mat.echelonize()
            mat=mat[:mat.rank()]
            self._extended_reln_matrix = mat
            self.__save()
        elif verbose:
            print('Extended relation matrix is immutable, so it\'s probably already reduced. Doing no more reduction.')
        return self._extended_reln_matrix

    def weight_denominator(self,verbose=False):
        return self.reduced_extended_reln_matrix(verbose=verbose)[-1,-1]
        
    def H_1(self,verbose=False):    # first homology group
        if self._H_1 == None:
            mat = self.reduced_reln_matrix(verbose=verbose)
            rank = len(mat.columns())-len(mat.rows()) # then rank of H_1
            self._reln_matrix_smith = mat.smith_form()
            self._H_1_elementary_divisors = self._reln_matrix_smith[0].T.elementary_divisors()
            self._H_1_reduced_elementary_divisors = [x for x in self._H_1_elementary_divisors if x !=1]
            self._H_1 = AdditiveAbelianGroup(self._H_1_reduced_elementary_divisors,remember_generators=True)
        return self._H_1

    def H_1_met(self,verbose=False):    # first homology of pre-image in covering group
        if self._H_1_met == None:
            mat = self.reduced_extended_reln_matrix(verbose=verbose)
            self._extended_reln_matrix_smith = mat.smith_form()
            self._H_1_met_elementary_divisors = self._extended_reln_matrix_smith[0].T.elementary_divisors()
            self._H_1_met = AdditiveAbelianGroup([x for x in self._H_1_met_elementary_divisors if x !=1])
        return self._H_1_met
    
    """
    convert a group element to an element of H_1.
    """
    def elt2H_1(self,g,verbose=False):
        g_word = self.elt2RSword(g)
        #print('word=',g_word)
        g_multiplicities = zero_vector(ZZ,self._no_of_gens)
        for letter in g_word:
            if letter >0:
                g_multiplicities[letter-1] +=1
            else:
                g_multiplicities[-letter-1] -=1
        #print('multiplicities:',g_multiplicities)
        self.H_1(verbose=verbose)
        g_smith = g_multiplicities * self._reln_matrix_smith[2]
        #print('Smith form:    ',g_smith)
        start_column = self._H_1_elementary_divisors.count(1)
        #print('removing nulls:',g_smith[start_column:])
        return self._H_1(g_smith[start_column:])
    
    """
    find a lift to the group of an element of H_1.
    """
    def H_1_lift(self,h_1_elt):
        self.H_1()
        new_h_1_elt = self._H_1(h_1_elt)
        start_column = self._H_1_elementary_divisors.count(1)
        g_smith = vector(ZZ,start_column*[0]+new_h_1_elt.lift().list())
        #print('     g_smith=',g_smith)
        g_multiplicities = g_smith * self._reln_matrix_smith[2]^-1
        #print('multilicities',g_multiplicities)
        return prod([self._RSgens[counter+1]^(g_multiplicities[counter]) for counter in range(self._no_of_gens)])

    """
    Find the matrix of the diamond operator on H_1 for a matrix g which normalizes self.
    Here self can be regarded as a subgroup of PSU(2,1)(alpha) as a quotient; g need not normalize self as
    a subgroup of SU(2,1).
    This code does not check whether g normalizes self, and will produce meaningless results if it does not.
    """
    def H_1_diamond(self,g):
        row_list=[]
        for gen in self.H_1().gens():
            #print(gen,end='')
            gen_lift = self.H_1_lift(gen)
            gen_lift_conj = SUalpha2PSUalpha(g^-1*gen_lift*g) # project into self as a quotient.
            #pretty_print((gen_lift_conj))
            #print(gen_lift_conj in self)
            gen_conj = self.elt2H_1(gen_lift_conj)
            gen_conj = gen_conj.lift()
            # the following loop should not be needed, but sage is lifting group elements in a weird way.
            for counter in range(len(self._H_1_reduced_elementary_divisors)):
                generator_order = self._H_1_reduced_elementary_divisors[counter]
                if generator_order!=0:
                    gen_conj[counter] =gen_conj[counter] % generator_order
            #print(gen_conj)
            row_list.append(list(gen_conj))
        return matrix(ZZ,row_list)
            
# Common arithmetic subgroups of SU(2,1)(OK,alpha).
# each of these function returns an instance of ArithGp.

def Gamma1(beta):
    name = 'Gamma1('+str(beta)+')'
    latex_name = LatexExpr('\\Gamma_1(')+latex(beta)+LatexExpr(')')
    def key(g):
        R_finite.<onemod,omegamod> = OK.quotient_ring([beta])
        homg = MatOK(g).row(2).change_ring(R_finite)
        return str(homg)
    return ArithGp(name,key,latex_name=latex_name)

def Gamma0(beta):
    name = 'Gamma0('+str(beta)+')'
    latex_name = LatexExpr('\\Gamma_0(')+latex(beta)+LatexExpr(')')
    def key(g):
        R_finite.<onemod,omegamod> = OK.quotient_ring([beta])
        homg = MatOK(g).row(2).change_ring(R_finite)
        try:
            homg = homg* homg[0]^-1
        except ZeroDivisionError:
            try:
                homg = homg* homg[1]^-1
            except ZeroDivisionError:
                homg = homg* homg[2]^-1
        return str(homg)
    return ArithGp(name,key,latex_name=latex_name)

def GammaNC(beta):
    name = 'GammaNC('+str(beta)+')'
    latex_name = LatexExpr('\\Gamma_{\\mathrm{NC}}(')+latex(beta)+LatexExpr(')')
    def key(g):
        R_finite.<onemod,omegamod> = OK.quotient_ring([beta])
        homg = R_finite(phi(MatOK(g)))
        return str(homg)
    return ArithGp(name,key,latex_name=latex_name)

def Gamma(beta):
    name = 'Gamma('+str(beta)+')'
    latex_name = LatexExpr('\\Gamma(')+latex(beta)+LatexExpr(')')
    def key(g):
        R_finite.<onemod,omegamod> = OK.quotient_ring([beta])
        homg = MatOK(g).change_ring(R_finite)
        return str(homg)
    return ArithGp(name,key,latex_name=latex_name)


def Index3congruence(a,b,c,d):
    name = 'GammaIndex3'+str((a,b,c,d))
    latex_name = LatexExpr('\\Gamma_{\\mathrm{index\\ 3}}')+latex((a,b,c,d))
    def key(g):
        homg = F3((a*g[0,1]+b*g[0,2]+c*g[1,0]+d*g[2,0])/alpha)
        return str(homg)
    return ArithGp(name,key,latex_name=latex_name)


for n in range(1,100):
    for beta in K.elements_of_norm(n):
        G = Gamma0(3) & GammaNC(beta)
        print(G.index(), G.H_1().short_name())


