import sympy as s

s.init_printing()

rho, u, v, E, p, gamma = s.symbols('rho u v E p gamma')
Q1, Q2, Q3, Q4 = s.symbols('Q1 Q2 Q3 Q4')

u_eq = Q2 / Q1
v_eq = Q3 / Q1
E_eq = Q4 / Q1
p_eq = (gamma-1)*(Q4 - Q1*0.5*(u_eq**2 + v_eq**2))

Q = s.Matrix([Q1, Q2, Q3, Q4])

#F = s.Matrix([Q2, (Q2**2)/Q1 + p_eq, (Q2*Q3)/Q1, (Q4 + p)*u_eq])
F = s.Matrix([Q1*u_eq,
	Q1*(u_eq**2) + p_eq,
	Q1*u_eq*v_eq,
	(Q1*E_eq + p)*u_eq
	])

G = s.Matrix([Q1*v_eq,
	Q1*u_eq*v_eq,
	Q1*(v_eq**2) + p_eq,
	(Q1*E_eq + p)*v_eq
	])

A = F.jacobian(Q)
B = G.jacobian(Q)

A, B = s.simplify(A), s.simplify(B)

s.pprint(A)
print()
s.pprint(B)
