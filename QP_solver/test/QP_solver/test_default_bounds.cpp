#include <cstdlib>
#include <cassert>
#include <iostream>

#include <CGAL/QP_models.h>
#include <CGAL/QP_options.h>
#include <CGAL/QP_functions.h>

// choose exact integral type
#ifndef CGAL_USE_GMP
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#else
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz ET;
#endif

int fl;
int fu;
int l;
int u;
int r;

// inline void check(bool cond)
// {
//   if (!cond) {
//     std::cout << "fl = " << fl << ", "
//               << " l = " <<  l << ", "
//               << "fu = " << fu << ", "
//               << " u = " <<  u << ", "
//               << " r = " <<  r << std::endl;
//   }
// }

int main()
{
  // options
  CGAL::Quadratic_program_options options;
  options.set_auto_validation(true);

  for (fl = 0; fl < 2; ++fl)
    for (fu = 0; fu < 2; ++fu)
      for (l = -2; l <= 0; ++l)
        for (u = 0; u <= 2; ++u)
          for (r = -1;  r < 1; ++r) {
            // generate empty program with all possible
            // defaults, and test solver / bound status functions
            CGAL::Quadratic_program<int> qp = CGAL::Quadratic_program<int>
              (CGAL::Comparison_result(r), (fl!=0), l, (fu!=0), u);
            assert(qp.is_valid());

            // test solver
            CGAL::Quadratic_program_solution<ET> s =
              CGAL::solve_quadratic_program (qp, ET(), options);
            assert(CGAL::is_zero(s.objective_value()));
            assert(s.is_optimal()); assert(s.is_valid());

            // test bounds (program is empty, so everything should hold)
            assert(qp.is_nonnegative());
            assert(qp.is_nonpositive());
            assert(qp.is_free());

            // now manipulate program
            qp.set_c(0, 1);        //       min x_0
                                   //      l <= x_0 <= u
            // test solver
            s = CGAL::solve_quadratic_program (qp, ET(), options);
            assert(s.is_valid());
            if (fl) {
              assert(s.is_optimal());
              assert(s.objective_value() == l);
            } else
              assert(s.is_unbounded());

            // test bounds
            if (fl && fu) {
              assert(!qp.is_nonnegative());
              assert(!qp.is_nonpositive());
              assert(!qp.is_free());
            }
            if (fl && !fu) {
              if (l==0)
                assert(qp.is_nonnegative());
              else
                assert(!qp.is_nonnegative());
              assert(!qp.is_nonpositive());
              assert(!qp.is_free());
            }
            if (!fl && fu) {
              if (u==0)
                assert(qp.is_nonpositive());
              else
                assert(!qp.is_nonpositive());
              assert(!qp.is_nonnegative());
              assert(!qp.is_free());
            }
            if (!fl && !fu) {
              assert(qp.is_free());
              assert(!qp.is_nonnegative());
              assert(!qp.is_nonpositive());
            }

            // manipulate program
            qp.set_l(0, true, 0);  //      0 <= x_0 <= u

            // test solver
            s = CGAL::solve_quadratic_program (qp, ET(), options);
            assert(s.is_valid());
            assert(s.is_optimal());
            assert(CGAL::is_zero(s.objective_value()));

            // test bounds
            if (!fu)
              assert(qp.is_nonnegative());
            else
              assert(!qp.is_nonnegative());
            assert(!qp.is_free());
            assert(!qp.is_nonpositive());

            // manipulate program
            qp.set_l(0, false);    // -infty <= x_0 <= u

            // test solver
            s = CGAL::solve_quadratic_program (qp, ET(), options);
            assert(s.is_valid());
            assert(s.is_unbounded());

            // manipulate program
            qp.set_c(0, -1);      //        min -x_0

            // test solver
            s = CGAL::solve_quadratic_program (qp, ET(), options);
            assert(s.is_valid());
            if (fu) {
              assert(s.is_optimal());
              assert(s.objective_value() == -u);
            } else
              assert(s.is_unbounded());

            // test bounds
            if (fu) {
              if (u==0)
                assert(qp.is_nonpositive());
              else
                assert(!qp.is_nonpositive());
              assert(!qp.is_free());
            } else {
              assert(qp.is_free());
              assert(!qp.is_nonpositive());
            }
            assert(!qp.is_nonnegative());

            // manipulate program
            qp.set_c0(5);
            qp.set_u(0, true, 0); //  -infty <= x_0 <= 0

            // test solver
            s = CGAL::solve_quadratic_program (qp, ET(), options);
            assert(s.is_valid());
            assert(s.is_optimal());
            assert(s.objective_value() == 5);

            // test bounds
            assert(qp.is_nonpositive());
            assert(!qp.is_nonnegative());
            assert(!qp.is_free());

            qp.set_u(0, false);  //   -infty <= x_0 <= infty
            s = CGAL::solve_quadratic_program (qp, ET(), options);
            assert(s.is_valid());
            assert(s.is_unbounded());
            assert(qp.is_free());
            assert(!qp.is_nonnegative());
            assert(!qp.is_nonpositive());

            // test default_r: insert constraint 0x ~ 1; this is
            // feasible iff <=
            qp.set_a(0,0,0);
            qp.set_b(0,1);

            // test solver
            s = CGAL::solve_quadratic_program (qp, ET(), options);
            assert(s.is_valid());
            if (r == -1)
              assert(s.is_unbounded());
            else
              assert(s.is_infeasible());
          }
  return 0;
}
