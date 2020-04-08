
/*
* This example shows how to formulate and solve the following MIP problem:
*
*        Maximize
*   Objective: x1 + 2 x2 + 3 x3 + x4
*   Subject to
*           c1: - x1 +   x2 + x3 +  10 x4 <= 20
*           c2:   x1 - 3 x2 + x3          <= 30
*           c3:          x2      - 3.5 x4  =  0
*   Bounds
*           0 <= x1 <= 40
*           2 <= x4 <= 3
*   General
*           x4 is integer
*
*           Expected results: x1=40; x2=10.5; x3=19.5; x4=3;
*/


#include <iostream>

#ifdef CGAL_USE_SCIP

#include <CGAL/SCIP_mixed_integer_program_traits.h>
typedef CGAL::SCIP_mixed_integer_program_traits<double>                        MIP_Solver;

#elif defined(CGAL_USE_GLPK)

#include <CGAL/GLPK_mixed_integer_program_traits.h>
typedef CGAL::GLPK_mixed_integer_program_traits<double>                        MIP_Solver;

#endif


#if defined(CGAL_USE_GLPK) || defined(CGAL_USE_SCIP)

typedef typename MIP_Solver::Variable                        Variable;
typedef typename MIP_Solver::Linear_objective        Linear_objective;
typedef typename MIP_Solver::Linear_constraint        Linear_constraint;


int main()
{
        MIP_Solver solver;

        // Variable x1
        Variable* x1 = solver.create_variable(Variable::CONTINUOUS, 0, 40, "x1");

        // Variable x2
        // You can create first (using default parameters) and then assign values.
        Variable* x2 = solver.create_variable();
        x2->set_name("x2");        // This is optional (a default will be given)

        // Variable x3
        Variable* x3 = solver.create_variable();        // Uses all default parameters
        x3->set_name("x3");

        // Variable x4
        Variable* x4 = solver.create_variable(Variable::INTEGER, 2, 3, "x4");

        // Objective.
        // Be careful this is "MAXIMIZE"
        Linear_objective * obj = solver.create_objective(Linear_objective::MAXIMIZE);
        obj->add_coefficient(x1, 1.0);
        obj->add_coefficient(x2, 2.0);
        obj->add_coefficient(x3, 3.0);
        obj->add_coefficient(x4, 1.0);

        // Constraint c1: -x1 + x2 + x3 + 10 x4 <= 20
        Linear_constraint* c1 = solver.create_constraint(-Linear_constraint::infinity(), 20, "c1");
        c1->add_coefficient(x1, -1);
        c1->add_coefficient(x2, 1);
        c1->add_coefficient(x3, 1);
        c1->add_coefficient(x4, 10);

        // Constraint c2: x1 - 3 x2 + x3 <= 30
        Linear_constraint* c2 = solver.create_constraint(-Linear_constraint::infinity(), 30, "c2");
        c2->add_coefficient(x1, 1);
        c2->add_coefficient(x2, -3);
        c2->add_coefficient(x3, 1);

        // Constraint c3:  x2 - 3.5 x4 = 0
        Linear_constraint* c3 = solver.create_constraint(0, 0, "c3");
        c3->add_coefficient(x2, 1);
        c3->add_coefficient(x4, -3.5);

        // Solve
        if (solver.solve()) {
                std::cout << "result: " << std::endl;
                const std::vector<double>& results = solver.solution();
                for (std::size_t i = 0; i < results.size(); ++i) {
                        std::cout << "\tx" << i + 1 << ": " << results[i] << std::endl;
                }
                return EXIT_SUCCESS;
        }
        else {
                std::cerr << "solving problem failed" << std::endl;
                return EXIT_FAILURE;
        }
}


#else

int main(int , char**)
{
    std::cerr << "This test requires either GLPK or SCIP.\n";
    return EXIT_SUCCESS;
}

#endif  // defined(CGAL_USE_GLPK) || defined(CGAL_USE_SCIP)

