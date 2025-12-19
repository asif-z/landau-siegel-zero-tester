//
// Created by kris on 12/19/25.
//

#include <flint/acb_dirichlet.h>
#include <flint/acb.h>
#include <flint/dirichlet.h>
#include <flint//arb.h>

int main(int argc, char* argv[])
{
    long q = 2;
    while (q<=1000000)
    {

        if (q%1000 == 0)
        {
            printf("======%ld======\n", q);
        }
        long prec = 30;

        arb_t rhs;
        arb_t log;
        arb_t den;
        arb_t eight;
        arb_init(rhs);
        arb_init(log);
        arb_init(den);
        arb_init(eight);
        arb_set_si(eight, 8);
        arb_log_ui(log, q, prec);
        arb_mul(den,eight,log,prec);
        arb_inv(rhs,den,prec);

        acb_t re;
        arb_t result;
        acb_t s;
        acb_init(re);
        acb_init(s);
        acb_one(s);
        arb_init(result);
        dirichlet_group_t G;
        dirichlet_group_init(G, q);
        dirichlet_char_t chi;
        dirichlet_char_init(chi, G);
        dirichlet_char_one(chi, G);
        do
        {
            /* use character chi */
            if (dirichlet_char_is_primitive(G, chi) == 1 && dirichlet_char_is_real(G, chi) == 1)
            {
                acb_dirichlet_l(re, s, G, chi, prec);
                acb_get_real(result, re);
                // arb_print(rhs);
                if (arb_ge(result,rhs)== 1)
                {
                    // printf("true @ %ld \n",q);
                } else
                {
                    printf("false @ %ld \n",q);
                }
            }
        }
        while (dirichlet_char_next(chi, G) >= 0);
        q++;
    }
}
