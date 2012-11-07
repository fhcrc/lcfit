#include <stdio.h>
#include "quintic_C.h"

void print_arr(double *a, int len)
{
    int i;

    for(i=0; i < len; i++) {
        printf("%g\t", a[i]);
    }
    printf("\n");
}

void run_test(double dd[])
{
double sol[4], soli[4];
int Nsol;

/*
 *     dd(0:4)     (i)  vector containing the polynomial coefficients
 *     sol(1:4)    (o)  results, real part
 *     soli(1:4)   (o)  results, imaginary part
 *     Nsol        (o)  number of real solutions
 */
quartic(dd, sol, soli, &Nsol);

printf("%d real solutions and %d complex solutions\n", Nsol, 4-Nsol);
printf("real part\n");
print_arr(sol,4);
printf("imaginary part\n");
print_arr(soli,4);
}

main()
{

/*
http://www.1728.org/quartic2.htm

Example #1:
3X^4   + 6X^3   - 123X^2   - 126X + 1,080 = 0

Solutions:
X1 = 5
X2 = 3
X3 = -4
X4 = -6

From Maxima:
solve(3*X^4 + 6*X^3 - 123*X^2 - 126*X + 1080, X);
[X=3,X=-6,X=5,X=-4]

4 real solutions and 0 complex solutions
real part
5	3	-4	-6
imaginary part
0	0	0	0
*/

double dd1[] = {1080, -126, -123, 6, 3};
run_test(dd1);

/*
http://www.1728.org/quartic2.htm

Example #2:
-20X^4   + 5X^3   + 17X^2   - 29X + 87 = 0

X1 = 1.48758311033
X2 = 0.222210408124 + i*1.29967219908
X3 = 0.222210408124 - i*1.29967219908
X4 = -1.68200392658

From Maxima:
solve(-20*X^4 + 5*X^3 + 17*X^2 - 29*X + 87, X),numer;
[X=-1.68200394973571,X=1.487583147671779,X=0.22221040103197-1.299672233975684*%i,X=1.299672233975684*%i+
0.22221040103197]

quintic_C:
2 real solutions and 2 complex solutions
real part
1.48758	-1.682	0.22221	0.22221
imaginary part
0	0	1.29967	-1.29967
*/

double dd2[] = {87, -29, 17, 5, -20};
run_test(dd2);

return 0;
}
