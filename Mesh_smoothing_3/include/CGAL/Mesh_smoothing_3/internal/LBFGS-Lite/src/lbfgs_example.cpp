#include "lbfgs.hpp"
#include <Eigen/Eigen>
#include <iostream>
#include <iomanip>

class MinimizationExample
{
public:
    int run(const int N)
    {
        double finalCost;
        Eigen::VectorXd x(N);

        /* Set the initial guess */
        for (int i = 0; i < N; i += 2)
        {
            x(i) = -1.2;
            x(i + 1) = 1.0;
        }

        /* Set the minimization parameters */
        lbfgs::lbfgs_parameter_t params;
        params.g_epsilon = 1.0e-8;
        params.past = 3;
        params.delta = 1.0e-8;

        /* Start minimization */
        int ret = lbfgs::lbfgs_optimize(x,
                                        finalCost,
                                        costFunction,
                                        nullptr,
                                        monitorProgress,
                                        this,
                                        params);

        /* Report the result. */
        std::cout << std::setprecision(4)
                  << "================================" << std::endl
                  << "L-BFGS Optimization Returned: " << ret << std::endl
                  << "Minimized Cost: " << finalCost << std::endl
                  << "Optimal Variables: " << std::endl
                  << x.transpose() << std::endl;

        return ret;
    }

private:
    static double costFunction(void *instance,
                               const Eigen::VectorXd &x,
                               Eigen::VectorXd &g)
    {
        const int n = x.size();
        double fx = 0.0;
        for (int i = 0; i < n; i += 2)
        {
            const double t1 = 1.0 - x(i);
            const double t2 = 10.0 * (x(i + 1) - x(i) * x(i));
            g(i + 1) = 20.0 * t2;
            g(i) = -2.0 * (x(i) * g(i + 1) + t1);
            fx += t1 * t1 + t2 * t2;
        }
        return fx;
    }

    static int monitorProgress(void *instance,
                               const Eigen::VectorXd &x,
                               const Eigen::VectorXd &g,
                               const double fx,
                               const double step,
                               const int k,
                               const int ls)
    {
        std::cout << std::setprecision(4)
                  << "================================" << std::endl
                  << "Iteration: " << k << std::endl
                  << "Function Value: " << fx << std::endl
                  << "Gradient Inf Norm: " << g.cwiseAbs().maxCoeff() << std::endl
                  << "Variables: " << std::endl
                  << x.transpose() << std::endl;
        return 0;
    }
};

int main(int argc, char **argv)
{
    MinimizationExample example;
    return example.run(200);
}
