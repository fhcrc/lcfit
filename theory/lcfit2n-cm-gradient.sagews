︠5763dd17-42da-4004-b396-236710dba9b4︠
#
# This workbook derives the gradient of the normalized lcfit2 log-likelihood function, implemented in lcfit2n_gradient() in lcfit2.c.
#

# Declare variables.
var('c, m, t, t_0')
var('f_2', latex_name = 'f^{\prime\prime}(t_0)')

︡f490c6cf-0e17-4443-a806-7492f40c5340︡
︠2a967c2f-d93c-4ad8-b388-ac237fc40674s︠

# Define v.
v(t) = (c - m) * exp((-2 * (t - t_0) / (c - m)) * sqrt((-f_2 * c * m) / (c + m))); show(v)

︡70bd00f8-822e-4b50-be4c-2afcd92ec8a7︡
︠a29b8de2-d484-4d0e-bf12-db5b50f00fbds︠

# Define f in terms of c, m, and v.
f(t) = c * log(c + m + v) + m * log(c + m - v) - (c + m) * log(2 * (c + m)); show(f)

︡7d0767e4-923e-4c7b-bdcb-6220c5e1a260︡
︠1db93691-ee79-4638-8b90-3774627f6176s︠

# Define the normalized function f_n.
f_n(t) = f(t) - f(t_0); show(f_n)

︡45bf16d8-0943-4d7c-8e44-a3df00cc32da︡
︠44b57572-58eb-4d59-bffd-52d17ba87148s︠

# Compute the derivatives of f_n.
dc = f_n.diff(c)
dm = f_n.diff(m)

︡1c9e4789-c2fa-4bfc-9ea9-342a94a9a0e1︡
︠827aa573-bd27-462f-b71e-ee6674e77f3f︠

# Sympify.
import sympy
s_dc = sympy.sympify(dc)
s_dm = sympy.sympify(dm)

# Reset v so we can use it again.
reset('v')

# Declare substitution variables.
var('theta_tilde', latex_name = '\widetilde{\\theta}')
var('r, z, v')

# Store expressions for our substitutions.
theta_tilde_inv_expr = exp((-2 * (t - t_0) / (c - m)) * sqrt((-f_2 * c * m) / (c + m)))
r_expr = (2 / (c - m)) * sqrt((-f_2 * c * m) / (c + m))
z_expr = (-f_2 * c * m) / (c + m)
v_expr = (c - m) / theta_tilde

show(1/theta_tilde == theta_tilde_inv_expr)
show(r == r_expr)
show(z == z_expr)
show(v == v_expr)

︡4b84b51c-34ea-440b-8949-0cb7e521a589︡
︠6d6660b8-ec59-4394-927b-ac90601ebc3fs︠

#
# Simplify df/dc.
#

# Here's what we're starting with.
show(s_dc._sage_())

︡713118b8-28bb-4f07-ba42-c4f559538b3d︡
︠9fba9b97-fa50-40d7-b2f0-bcc8b497cc61s︠

# Substitute theta.
s_dc1 = s_dc.subs(theta_tilde_inv_expr, 1 / theta_tilde)
show(s_dc1._sage_())

︡496e9231-3b76-4afd-ab18-63d5aaa9acd8︡
︠22ab83be-a7c6-4882-be6a-106c8896bd32s︠

# Substitute r.
s_dc2 = s_dc1.subs(r_expr, r)
show(s_dc2._sage_())

︡0b8f4c50-c42a-4769-ace5-4c07749c9f30︡
︠7fb3f03e-554c-4dc8-baad-60a67d803ceb︠

# Substitute the radicand of r as z.
s_dc3 = s_dc2.subs(z_expr, z)
show(s_dc3._sage_())

︡ea51d659-ca92-40de-9379-bed9b848946b︡
︠edcb4529-de49-4a23-a38e-84aaa90a15e6s︠

# Substitute terms similar to the radicand z.
s_dc4 = s_dc3.subs((f_2 * m) / (c + m), -z / c)
show(s_dc4._sage_())

︡ee47debc-242e-438e-9d43-d622c7d77074︡
︠1d2e5e1e-3e94-40e5-8c1a-ddb96e6945d2s︠

# Substitute v.
s_dc5 = s_dc4.subs(v_expr, v)
show(s_dc5._sage_())

dc_simplified = s_dc5._sage_()

#
# df/dc =
#

︡2184ed02-42cb-4175-8e76-d113d1f7ffbf︡
︠f1695748-f853-42ba-bfce-2623b2f7a107s︠

#
# Simplify df/dm.
#

# Here's what we're starting with.
show(s_dm._sage_())

︡8e56edb8-6253-44cf-881c-e0635b8080cc︡
︠52b3da4a-5319-4c19-bf33-bc1e53e4b436s︠

# Substitute theta.
s_dm1 = s_dm.subs(theta_tilde_inv_expr, 1 / theta_tilde)
show(s_dm1._sage_())

︡13c43968-bcb3-418f-b899-155f92e12cb3︡
︠20d367b6-de70-4a07-b74d-55cee60ccf1as︠

# Substitute r.
s_dm2 = s_dm1.subs(r_expr, r)
show(s_dm2._sage_())

︡2a2fc7c0-8e2e-4738-acf8-e28777b0c4e7︡
︠bdcd1904-0d86-41ad-817c-907384d3943bs︠

# Substitute the radicand of r as z.
s_dm3 = s_dm2.subs(z_expr, z)
show(s_dm3._sage_())

︡7fa8ba89-7258-4107-988e-d361225a3e28︡
︠d0cde77f-2a0b-4349-90f6-6f51cb3f06e9s︠

# Substitute terms similar to the radicand z.
# Note that this substitution differs slightly from the one used for df/dc.
s_dm4 = s_dm3.subs((f_2 * c) / (c + m), -z / m)
show(s_dm4._sage_())

︡562f3677-34f6-490e-beed-f0c4ad1777a3︡
︠e062879a-8d21-4f07-a963-e080d65d2547s︠

# Substitute v.
s_dm5 = s_dm4.subs(v_expr, v)
show(s_dm5._sage_())

dm_simplified = s_dm5._sage_()

#
# df/dm =
#
︡1f5bf09f-6451-4e05-b995-4c62c0b03122︡
︠66b8da17-b91c-4eaf-b53e-5305c7f0ac0es︠

#
# Summary
#

show(dc_simplified)

︡cc0da9bd-cdfc-4bd1-a32a-3ad6cf544cab︡
︠ddef3772-6a82-4dc8-a37e-f9a1ac9ff996s︠

dc_simplified

︡d134cb8b-38e7-40b3-8602-ef1ff3257555︡
︠eda1741d-8b14-479e-abfd-035cda5b530cs︠

# For some reason sympy failed to substitute (c - m) / theta_tilde = v in the last log term.
show(dm_simplified)

︡03ad4383-f9ad-431f-95ca-5ce5eab8f7ab︡
︠dfda2512-54d9-4ffb-acdd-a976f5c562f4s︠

dm_simplified

︡b3c44b60-b2b2-4643-80c1-b5c2f0dd21e2︡









