# Copyright (c) 2022 Joshua Brakensiek, Manik Dhar, Sivakanth Gopi
#
# This code is licensed under the MIT License
#
# To run this code, run "julia claim_5_2.jl"
#
# The Oscar package needs to be installed first:
# https://oscar.computeralgebra.de/install/

using Oscar;

# suffices to verify identity in characteristic 0
k = QQ;
K, (gamma,x1,x2,x3,x4,x5,x6) = polynomial_ring(k, ["y","x_1","x_2","x_3","x_4","x_5","x_6"]; internal_ordering=:degrevlex);

# defining the Reed-Solomon generators

x = [x1,x2,x3,x4,x5,x6];
b = [x[i]+gamma*x[i]*x[i] for i in 1:6];

# building the polynomials p

M = matrix(K, [[1,b[1]+b[2],b[1]*b[2]],[1,b[3]+b[4],b[3]*b[4]],[1,b[5]+b[6],b[5]*b[6]]]);
q = det(M);

p0 = coeff(q, [gamma], [0])
p1 = coeff(q, [gamma], [1])
p2 = coeff(q, [gamma], [2])
p3 = coeff(q, [gamma], [3])

g = x2 * x3 + x2 * x4 - x2 * x5 - x2 * x6 - x3 * x4 + x5 * x6;

Q0 = (QQ(1)/QQ(2))*
    (- 2x1*x2^3*x3 - 2x1*x2^3*x4 + 2*x1*x2^3*x5 + 2x1*x2^3*x6 - x1*x2^2*x3^2 - x1*x2^2*x4^2 + x1*x2^2*x5^2
     + x1*x2^2*x6^2 + 2x1*x2*x3^2*x4 + 2x1*x2*x3*x4^2 - 2x1*x2*x5^2*x6 - 2x1*x2*x5*x6^2 - x1*x3^2*x4^2
     + x1*x5^2*x6^2 - 2x2^3*x3^2 - 2x2^3*x3*x4 - 2x2^3*x4^2 + 2x2^3*x5^2 + 2x2^3*x5*x6
     + 2x2^3*x6^2 - x2^2*x3^2*x4 - x2^2*x3*x4^2 + x2^2*x3*x4*x5 + x2^2*x3*x4*x6 - x2^2*x3*x5*x6
     - x2^2*x4*x5*x6 + x2^2*x5^2*x6 + x2^2*x5*x6^2 + 3x2*x3^2*x4^2 + x2*x3^2*x4*x5 + x2*x3^2*x4*x6
     - x2*x3^2*x5*x6 + x2*x3*x4^2*x5 + x2*x3*x4^2*x6 + x2*x3*x4*x5^2 + x2*x3*x4*x6^2
     - x2*x3*x5^2*x6 - x2*x3*x5*x6^2 - x2*x4^2*x5*x6 - x2*x4*x5^2*x6 - x2*x4*x5*x6^2
     - 3*x2*x5^2*x6^2 - x3^2*x4^2*x5 - x3^2*x4^2*x6 + x3^2*x4*x5*x6 + x3*x4^2*x5*x6
     - x3*x4*x5^2*x6 - x3*x4*x5*x6^2 + x3*x5^2*x6^2 + x4*x5^2*x6^2);
Q1 = 0;
Q2 = g * x2;
Q3 = -g / QQ(2);

S = [(3,6),(2,4),(2,6),(3,5),(4,5),(4,6),(2,5),(2,3)]

h = prod([x[i] - x[j] for (i,j) in S]);

r = h - p0*Q0 - p1*Q1 - p2*Q2 - p3*Q3;

println("h - p0*Q0 - p1*Q1 - p2*Q2 - p3*Q3: ", r);
