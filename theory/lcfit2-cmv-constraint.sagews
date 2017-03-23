︠40a90f05-9575-4085-9fb9-d371ffe05c16︠
#
# This workbook derives the NLopt constraint function and its gradient for enforcing that c + m - v > 0.
# The constraint is implemented in lcfit2_cons_cmv_nlopt() in lcfit2_nlopt.c.
#

def print_c(fn):
    from sympy.utilities.codegen import codegen

    fn_s = fn._sympy_()
    [(_, fn_code), (_, _)] = codegen(("fn", fn_s), "C", "test")
    print fn_code

︡44b64a98-29d2-4f72-866d-2ee5c7a1fe5d︡
︠b8cade84-70a9-40cd-a5e1-64757d3667c3︠

var('c, m, t_0')
var('f_2', latex_name = 'f^{\prime\prime}(t_0)')

r = (2 / (c - m)) * sqrt((-f_2 * c * m) / (c + m))

# we need t_0 <= (1/r)*log((c+m)/(c-m))
# nlopt expects constraints of the form fc <= 0
fcons = t_0 - (1/r)*log((c+m)/(c-m))

︡ae201671-6928-4db4-ac12-06fea6b3f68b︡
︠f2586ada-7346-401d-a6f6-09233c2cce6a︠

print_c(fcons)

︡b0274d70-fdfb-4cd3-8315-5477b301e957︡
︠c347345e-a65e-4fc0-bab0-2fcb4dd46efds︠

dc = fcons.diff(c)
print_c(dc)


︡9a388780-361c-416b-b5ed-190adde8ee84︡
︠d14b9efb-66c8-4f2f-9ade-fd0f895b1a0c︠

dm = fcons.diff(m)
print_c(dm)

︡acf9f3df-6f22-471a-8d2e-839095653b2d︡









