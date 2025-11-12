#!/usr/bin/env sage

print("Initializing Octonion Algebra over the Symbolic Ring (SR)...")
# Define the Octonion algebra over the Symbolic Ring
O = OctonionAlgebra(SR, -1, -1, -1)

# basis[0] is 1, basis[1:] are e1, ..., e7
basis = O.basis()
one = basis[0]
e = list(basis)[1:] # e[0] is e1, e[1] is e2, ..., e[6] is e7

# Define the specific basis vector e1
e1_vec = e[0]

print("Defining symbolic vectors X, Y, Z, H...")
# Define symbolic variables for the components
#try:
#    x_vars = list(var(','.join(['x%d' % i for i in range(1, 8)])))
#    y_vars = list(var(','.join(['y%d' % i for i in range(1, 8)])))
#    z_vars = list(var(','.join(['z%d' % i for i in range(1, 8)])))
#    h_vars = list(var(','.join(['h%d' % i for i in range(1, 8)])))
#except TypeError:
    # Handle case where vars might already be defined in a session
#    x_vars = [SR.symbol('x%d' % i) for i in range(1, 8)]
#    y_vars = [SR.symbol('y%d' % i) for i in range(1, 8)]
#    z_vars = [SR.symbol('z%d' % i) for i in range(1, 8)]
#    h_vars = [SR.symbol('h%d' % i) for i in range(1, 8)]

#if working over Im O, remember to project off real component

X = e[2]
Y = e[3]
Z = e[4]
H = e[5]

def ip(U, V):
    return (U * V.conjugate()).real_part()

def cross_product(U, V):
    return (U * V).imag_part()

def triple_cross(U, V, W):
    return (U * (V.conjugate() * W) - W * (V.conjugate() * U)) / 2

def B(U, V):
    return ip(U, H) * V + ip(V, H) * U - ip(U, V) * H

def J(U):
    return cross_product(e1_vec, U)

def C(U, V):
    return triple_cross(e1_vec, U, V)

print("Calculating terms (this may take a while)...")
result1 = J(B(X, Y)) - B(X, J(Y))
result2 = B(X, C(Y, Z)) - C(B(X, Y), Z) - C(Y, B(X, Z))

print("Result 1: ")
print(result1)
print("Result 2: ")
print(result2)