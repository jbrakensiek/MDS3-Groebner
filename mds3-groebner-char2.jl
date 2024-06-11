# Copyright (c) 2022 Joshua Brakensiek, Manik Dhar, Sivakanth Gopi
#
# This code is licensed under the MIT License
#
# To run this code, run "julia mds3-groebner-char2.jl"
#
# The Oscar package needs to be installed first:
# https://oscar.computeralgebra.de/install/

using Oscar;

# our construction works over characteristic 7
k = GF(2);
K, (gamma,a1,a2,a3,a4,a5,a6) = polynomial_ring(k, ["x","a1","a2","a3","a4","a5","a6"], internal_ordering=:degrevlex);

# defining the Reed-Solomon generators

a = [a1,a2,a3,a4,a5,a6];
b = [a[i]+gamma*a[i]*a[i]*a[i] for i in 1:6];

# building the polynomials p

M = matrix(K, [[1,b[1]+b[2],b[1]*b[2]],[1,b[3]+b[4],b[3]*b[4]],[1,b[5]+b[6],b[5]*b[6]]]);
q = det(M);
p = [coeff(q, [gamma], [i]) for i in 0:3]; # extracts the coefficients of gamma

for i in 0:3
    println("p", i," = ", p[i+1]) # julia is 1-indexed
end
println()

# building the ideal
Il = [p[1], p[2], p[3], p[4]] 
I = ideal(K, Il)

# compute the groebner basis
G,m = groebner_basis_with_transformation_matrix(I, ordering=degrevlex(K), complete_reduction=true);

# generate the polynomial we seek to prove is in the ideal
key = prod([a[i] - a[j] for i in 1:6 for j in 1:(i-1)])*prod([a[i] + a[j] + a[k] for i in 1:6 for j in 1:(i-1) for k in 1:(j-1)]);
println("Checking if ", factor(key), " is in the ideal")

# prove that key is in the ideal
Q, r = reduce_with_quotients(key, gens(G));
println("Remainder: ", r); # should be 0
   
# Optional: if you want to see the certificate for why key is in the ideal
#for z in Q*transpose(m)
#    println(factor(z))
#end
