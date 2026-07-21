#pragma once

#include <functional>

namespace Mesh_smoothing_3_internal {

#include "LBFGS-Lite/include/lbfgs.hpp" // todo: probably not good inside namespace. modify the code itself? Without namespace?


class Function_minimizer {
public:
    using Func_eval = std::function<double(Eigen::VectorXd const &x, Eigen::VectorXd &g)>;
    // return true to stop the optimization
    using Progress_call_back = std::function<bool(
        Eigen::VectorXd const &x,
        Eigen::VectorXd const &g,
        double f,
        double step,
        unsigned iter,
        unsigned nbEval
    )>;

    using Max_LS_step_bound = std::function<double(
        Eigen::VectorXd const &x,
        Eigen::VectorXd const &direction
    )>;

    using Direction_Precond = std::function<Eigen::VectorXd(
        Eigen::VectorXd const &x,
        Eigen::VectorXd const &direction
    )>;

    Function_minimizer(const Func_eval f) : _f(f) {}

    Func_eval const _f;
    Progress_call_back _call_back = [](Eigen::VectorXd const &,
                                       Eigen::VectorXd const &,
                                       const double, const double, const unsigned, const unsigned) { return false; };
    Max_LS_step_bound _max_step = nullptr;
    Direction_Precond _precond = nullptr;
    lbfgs::lbfgs_parameter_t _parameters;
    bool _has_warmstart = false;

    void set_max_iter(unsigned nb_iter) {_parameters.max_iterations = static_cast<int>(nb_iter); }

    void set_init_step(double step) { _parameters.init_step = step; }

    // 0 is success
    int get_status() const {return _status;}
    std::string get_message() const {
        if (_status == -1)
            return "Not run.";
        else
            return lbfgs::lbfgs_strerror(_status);
    }

    bool lbfgs_optimize(Eigen::VectorXd &x) {
        _status = lbfgs::lbfgs_optimize(x, _fval, static_function_call, static_bound_call, static_progress_call, static_cast<void*>(this), _parameters, static_direction_precond);
        return _status >= 0;
    }

    bool line_search(Eigen::VectorXd &x, Eigen::VectorXd const &d, double &step) {
        lbfgs::callback_data_t cd {
            static_cast<void*>(this),
            static_function_call,
            static_bound_call,
            static_progress_call
        };
        Eigen::VectorXd g = Eigen::VectorXd::Zero(x.size());
        double f = _f(x,g);
        Eigen::VectorXd gp = g;
        Eigen::VectorXd xp = x;
        _status = lbfgs::line_search_lewisoverton(x, f, g, step, d, xp, gp, _parameters.min_step, _parameters.max_step, cd, _parameters);
        return _status >= 0;
    }

    double minimum() const { return _fval; }

private:
    int _status = -1;
    double _fval = 0.;
    static double static_function_call(void *instance, Eigen::VectorXd const &x, Eigen::VectorXd &g) {
        Function_minimizer const * wrapper = static_cast<Function_minimizer const *>(instance);
        return wrapper->_f(x, g);
    }
    static int static_progress_call(void *instance,
                                const Eigen::VectorXd &x,
                                const Eigen::VectorXd &g,
                                const double fx,
                                const double step,
                                const int k,
                                const int ls)
    {
        Function_minimizer const * wrapper = static_cast<Function_minimizer const *>(instance);
        return wrapper->_call_back(x,g,fx,step, static_cast<unsigned>(k), static_cast<unsigned>(ls));
    }
    static double static_bound_call(void *instance,
        const Eigen::VectorXd &x,
        const Eigen::VectorXd &d)
    {
        Function_minimizer const * wrapper = static_cast<Function_minimizer const *>(instance);
        return wrapper->_max_step ? wrapper->_max_step(x,d) : wrapper->_parameters.max_step;
    }
    static Eigen::VectorXd static_direction_precond(void *instance,
        const Eigen::VectorXd &x,
        const Eigen::VectorXd &d)
    {
        Function_minimizer const * wrapper = static_cast<Function_minimizer const *>(instance);
        return wrapper->_precond ? wrapper->_precond(x,d) : d;
    }
};

}

