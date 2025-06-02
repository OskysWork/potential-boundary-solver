import sympy as s

s.init_printing()

rho, u, v, E, p, gamma = s.symbols('rho u v E p gamma')
x, y = s.symbols('x y')
ut, vt, Et = s.symbols('ut vt Et')
ug, vg, Eg = s.symbols('ug vg Eg')
Q1, Q2, Q3, Q4 = s.symbols('Q1 Q2 Q3 Q4')

ug = s.Function('ug')(x, y)
vg = s.Function('vg')(x, y)
Eg = s.Function('Eg')(x, y)

ut = s.Function('ut')(x, y)
vt = s.Function('vt')(x, y)
Et = s.Function('Et')(x, y)

u_eq = Q2 / Q1
v_eq = Q3 / Q1
E_eq = Q4 / Q1
p_eq = (gamma-1)*(Q4 - Q1*0.5*(u_eq**2 + v_eq**2))

Q = s.Matrix([Q1, Q2, Q3, Q4])

#Q_subs = {Q1: rho, Q2: rho*u, Q3: rho*v, Q4: rho*E}
Q_subs_comp = {Q1: rho, Q2: rho*(ut+ug), Q3: rho*(vt+vg), Q4: rho*(Et + Eg)}
Q_subs_t = {Q1: rho, Q2: rho*ut, Q3: rho*vt, Q4: rho*Et}
Q_subs_g = {Q1: rho, Q2: rho*ug, Q3: rho*vg, Q4: rho*Eg}

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

A_inv = A.inv()

A, B = s.simplify(A), s.simplify(B)

Q_dx = s.Matrix([
	s.Derivative(Q1, x),
	s.Derivative(Q2, x),
	s.Derivative(Q3, x),
	s.Derivative(Q4, x)
	])

Q_dy = s.Matrix([
	s.Derivative(Q1, y),
	s.Derivative(Q2, y),
	s.Derivative(Q3, y),
	s.Derivative(Q4, y)
	])

Q_dx_g = Q_dx.subs(Q_subs_g).doit()
Q_dy_g = Q_dy.subs(Q_subs_g).doit()

Q_dx_t = Q_dx.subs(Q_subs_t).doit()
Q_dy_t = Q_dy.subs(Q_subs_t).doit()

A_t = A.subs(Q_subs_t)
B_t = B.subs(Q_subs_t)

source_A = A_t*Q_dx_g
source_B = B_t*Q_dy_g
source = source_A + source_B

source_A = s.simplify(source_A)
source_B = s.simplify(source_B)
source = s.simplify(source)




A_comp = A.subs(Q_subs_comp)
B_comp = B.subs(Q_subs_comp)

A_g = A.subs(Q_subs_g)
B_g = B.subs(Q_subs_g)

A_tg = A_t + A_g
B_tg = B_t + B_g

A_comp, B_comp = s.simplify(A_comp), s.simplify(B_comp)
A_tg, B_tg = s.simplify(A_tg), s.simplify(B_tg)

A_inv_g = A_inv.subs(Q_subs_g)
F_comp = F.subs(Q_subs_comp)
#Qt = A_inv_g*F_comp
#Qt = s.simplify(Qt)

Qt_x = A_t * Q_dx_t
Qg_x = A_g * Q_dx_g





s.pprint(A_t)
print()
s.pprint(B_comp)
print("\nSource:\n")
s.pprint(source)
print("\nQt:\n")
s.pprint(Qt_x)
print("\nQg:\n")
s.pprint(Qg_x)
print("\nQ_dx:\n")
s.pprint(Q_dx_t)
