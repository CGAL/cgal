// Copyright (c) 2018  Liangliang Nan. All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Liangliang Nan

#ifndef CGAL_MIXED_INTEGER_PROGRAM_TRAITS_H
#define CGAL_MIXED_INTEGER_PROGRAM_TRAITS_H

#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

namespace CGAL {

        /// \cond SKIP_IN_MANUAL
        // forward declaration
        template <typename FT>
        class Mixed_integer_program_traits;

        /// The base class of solver element(e.g., Variable, Linear_constraint, and Linear_objective) in
        /// a mixed integer program

        template <typename FT>
        class Solver_entry
        {
        public:
                typedef        CGAL::Mixed_integer_program_traits<FT>        Solver;

        private:
                /// A solver element (e.g., variable, constraint, objective) cannot belong to multiple solvers.
                /// "solver" owns this entry.
                Solver_entry(
                        Solver* solver,
                        const std::string& name = "",
                        int idx = 0
                ) : solver_(solver), name_(name), index_(idx) {}

        public:
                const std::string& name() const { return name_; }
                void set_name(const std::string& n) { name_ = n; }

                int  index() const { return index_; }
                void set_index(int idx) { index_ = idx; }

                /// The solver that owns this entry
                const Solver* solver() const { return solver_; }
                Solver* solver() { return solver_; }

        private:
                Solver*                solver_; // The solver that owns this entry
                std::string        name_;
                int                        index_;

                template <typename T> friend class Variable;
                template <typename T> friend class Linear_expression;
                template <typename T> friend class Mixed_integer_program_traits;
        };


        /// The base class of solver elements that have bound constraints.
        template <typename FT>
        class Bound
        {
        private:
                Bound(FT lb = -infinity(), FT ub = +infinity());

        public:
                void set_bounds(FT lb, FT ub);
                void set_lower_bound(FT lb) { lower_bound_ = lb; }
                void set_upper_bound(FT ub) { upper_bound_ = ub; }

                void get_bounds(FT& lb, FT& ub) const;
                FT lower_bound() const { return lower_bound_; }
                FT upper_bound() const { return upper_bound_; }

                static FT infinity();

        private:
                FT                lower_bound_;
                FT                upper_bound_;

                static FT infinity_;

                template <typename T> friend class Variable;
                template <typename T> friend class Linear_constraint;
                template <typename T> friend class Mixed_integer_program_traits;
        };
        /// \endcond

        /// \ingroup PkgSolverInterfaceRef
        ///
        /// The variable of a mixed integer program.
        ///
        /// \cgalModels `MixedIntegerProgramVariable`
        template <typename FT>
        class Variable : public Solver_entry<FT>, public Bound<FT>
        {
                /// \cond SKIP_IN_MANUAL
        public:
                enum Variable_type { CONTINUOUS, INTEGER, BINARY };

                typedef        CGAL::Bound<FT>                             Bound;
                typedef        CGAL::Solver_entry<FT>                      Solver_entry;
                typedef        CGAL::Mixed_integer_program_traits<FT>      Solver;

        private:
                /// A variable cannot belong to several solvers.
                /// "solver" owns this variable.
                Variable(
                        Solver* solver,
                        Variable_type type = CONTINUOUS,
                        FT lb = -Bound::infinity(),
                        FT ub = +Bound::infinity(),
                        const std::string& name = "",
                        int idx = 0
                );

        public:
                Variable_type variable_type() const { return variable_type_; }
                void set_variable_type(Variable_type t);

                /// Returns the value of the variable in the current solution.
                /// \note (1) Valid only if the program was successfully solved.
                ///       (2) If the variable is integer and rounded == true, then the
                ///           value will be rounded to the nearest integer.
                FT solution_value(bool rounded = false) const;

                void set_solution_value(FT value) { solution_value_ = value; }

        private:
                Variable_type        variable_type_;
                FT                solution_value_;

                template <typename T> friend class Mixed_integer_program_traits;
                /// \endcond
        };


        /// \cond SKIP_IN_MANUAL
        /// The base class of Linear_constraint and Linear_objective.
        template <typename FT>
        class Linear_expression : public Solver_entry<FT>
        {
        public:
                typedef CGAL::Solver_entry<FT>                  Solver_entry;
                typedef CGAL::Variable<FT>                        Variable;
                typedef        CGAL::Mixed_integer_program_traits<FT>        Solver;

        private:
                /// An expression cannot belong to several solvers.
                /// "solver" owns this expression.
                Linear_expression(
                        Solver* solver,
                        const std::string& name = "",
                        int idx = 0
                );

        public:
                /// Adds a coefficient to a variable.
                void  add_coefficient(const Variable* var, FT coeff);

                const std::unordered_map<const Variable*, FT>& coefficients() const { return coefficients_; }
                void  set_coefficients(const std::unordered_map<const Variable*, FT>& coeffs) { coefficients_ = coeffs; }

                FT get_coefficient(const Variable* var) const;

                // The constant term
                void set_offset(FT value) { offset_ = value; }
                FT offset() const { return offset_; }

                /// Evaluates the value of this expression at the solution found.
                /// \note (1) valid only if the problem was successfully solved.
                ///       (2) if a variable is integer and rounded == true, then the
                ///           variable value will be rounded to the nearest integer.
                FT solution_value(bool rounded = false) const;

                virtual void clear() { coefficients_.clear(); offset_ = 0.0; }

        private:
                std::unordered_map<const Variable*, FT>        coefficients_;
                FT        offset_;

                template <typename T> friend class Linear_constraint;
                template <typename T> friend class Linear_objective;
                template <typename T> friend class Mixed_integer_program_traits;
        };
        /// \endcond

        /// \ingroup PkgSolverInterfaceRef
        ///
        /// The linear constraint of a mixed integer program.
        ///
        /// \cgalModels `MixedIntegerProgramLinearConstraint`
        template <typename FT>
        class Linear_constraint : public Linear_expression<FT>, public Bound<FT>
        {
                /// \cond SKIP_IN_MANUAL
        public:
                typedef        CGAL::Bound<FT>                         Bound;
                typedef        CGAL::Linear_expression<FT>             Linear_expression;
                typedef        CGAL::Mixed_integer_program_traits<FT>        Solver;

        private:
                /// A constraint cannot belong to several solvers.
                /// The "solver" owns this constraint.
                Linear_constraint(
                        Solver* solver,
                        FT lb = -Bound::infinity(),
                        FT ub = +Bound::infinity(),
                        const std::string& name = "",
                        int idx = 0
                );

                virtual ~Linear_constraint() {}

                template <typename T> friend class Mixed_integer_program_traits;
                /// \endcond
        };


        /// \ingroup PkgSolverInterfaceRef
        ///
        /// The linear objective of a mixed integer program.
        ///
        /// \cgalModels `MixedIntegerProgramLinearObjective`
        ///
        template <typename FT>
        class Linear_objective : public Linear_expression<FT>
        {
                /// \cond SKIP_IN_MANUAL
        public:
                typedef        CGAL::Mixed_integer_program_traits<FT>      Solver;
                typedef CGAL::Linear_expression<FT>                 Linear_expression;
                typedef typename Linear_expression::Solver_entry    Solver_entry;

                enum Sense { MINIMIZE, MAXIMIZE, UNDEFINED };

        private:
                /// An objective cannot belong to several solvers.
                /// "solver" owns this objective.
                Linear_objective(Solver* solver, Sense sense);
                virtual ~Linear_objective() {}

        public:
                void  set_sense(Sense sense) { sense_ = sense; }
                Sense sense() const { return sense_; }

                void clear();

        private:
                Sense sense_;

                template <typename T> friend class Mixed_integer_program_traits;
                /// \endcond
        };

        /// \ingroup PkgSolverInterfaceRef
        ///
        /// The class `CGAL::Mixed_integer_program_traits` provides an interface for
        /// formulating and solving (constrained or unconstrained) mixed integer
        /// programs. It can also be used for general linear programs.
        /// \note The solve() function is virtual and thus this class cannot be
        ///                  instantiated directly. Client code should use the inherited
        ///       classes, i.e., `CGAL::GLPK_mixed_integer_program_traits` or
        ///                  `CGAL::SCIP_mixed_integer_program_traits`. Alternatively, use
        ///       `CGAL::Mixed_integer_program_traits` as a base to derive a new model
        ///       (using e.g., <a href = "https://projects.coin-or.org/Cbc"> CBC </a>,
        ///       <a href = "http://www.gurobi.com/"> Gurobi </a> for better
        ///       performance).
        ///
        /// \cond SKIP_IN_MANUAL
        /// \tparam FT Number type
        /// \endcond
        ///
        /// \cgalModels `MixedIntegerProgramTraits`
        template <typename FT>
        class Mixed_integer_program_traits
        {
                /// \cond SKIP_IN_MANUAL
        public:
                typedef CGAL::Variable<FT>                                Variable;
                typedef CGAL::Linear_constraint<FT>                        Linear_constraint;
                typedef CGAL::Linear_objective<FT>                        Linear_objective;
                typedef typename Linear_objective::Sense                Sense;
                typedef typename Variable::Variable_type                Variable_type;

        public:
                Mixed_integer_program_traits();
                ~Mixed_integer_program_traits();

                /// Creates a single variable, add it to the solver, and returns its pointer.
                /// \note (1) If name is empty or not provided, a default name (e.g., x0, x1...)
                ///                  will be given.
                ///                  (2) Memory is managed by the solver and will be automatically released
                ///                  when the solver is destroyed.
                Variable* create_variable(
                        Variable_type type = Variable::CONTINUOUS,
                        FT lb = -Variable::infinity(),
                        FT ub = +Variable::infinity(),
                        const std::string& name = ""
                );

                /// Creates a set of variables and add them to the solver.
                /// \note (1) Variables will be given default names, e.g., x0, x1...
                ///                  (2) Memory is managed by the solver and will be automatically released
                ///                  when the solver is destroyed.
                std::vector<Variable*> create_variables(std::size_t n);

                /// Creates a single linear constraint, add it to the solver, and returns the pointer.
                /// \note (1) If 'name' is empty or not provided, a default name (e.g., c0, c1...) will be given.
                ///                  (2) Memory is managed by the solver and will be automatically released when the
                ///                  solver is destroyed.
                Linear_constraint* create_constraint(
                        FT lb = -Variable::infinity(),
                        FT ub = +Variable::infinity(),
                        const std::string& name = ""
                );

                /// Creates a set of linear constraints and add them to the solver.
                /// \note (1) Constraints with be given default names, e.g., c0, c1...
                ///                  (2) Memory is managed by the solver and will be automatically released
                ///                  when the solver is destroyed.
                std::vector<Linear_constraint*> create_constraints(std::size_t n);

                /// Creates the objective function and returns the pointer.
                /// \note Memory is managed by the solver and will be automatically released
                ///                  when the solver is destroyed.
                Linear_objective* create_objective(Sense sense = Linear_objective::MINIMIZE);

                std::size_t number_of_variables() const { return variables_.size(); }
                const std::vector<Variable*>& variables() const { return variables_; }
                std::vector<Variable*>& variables() { return variables_; }

                std::size_t number_of_constraints() const { return constraints_.size(); }
                const std::vector<Linear_constraint*>& constraints() const { return constraints_; }
                std::vector<Linear_constraint*>& constraints() { return constraints_; }

                const Linear_objective* objective() const;
                Linear_objective* objective();

                std::size_t number_of_continuous_variables() const;
                std::size_t number_of_integer_variables() const;
                std::size_t number_of_binary_variables() const;

                bool is_continuous() const;                                // Returns true if all variables are continuous
                bool is_mixed_integer_program() const;        // Returns true if mixed integer program
                bool is_integer_program() const;                // Returns true if integer program
                bool is_binary_program() const;                        // Returns true if binary program

                /// Clears all variables, constraints, and the objective.
                void clear();

                /// Solves the program. Returns false if fails.
                virtual bool solve() = 0;

                /// Returns the result.
                /// The result can also be retrieved using Variable::solution_value().
                /// \note (1) Result is valid only if the solver succeeded.
                ///       (2) Each entry in the result corresponds to the variable with the
                ///                         same index in the program.
                const std::vector<FT>& solution() const { return result_; }

                /// Returns the error message.
                /// \note This function should be called after call to solve().
                const std::string& error_message() const { return error_message_; }

        protected:
                Linear_objective*                                objective_;
                std::vector<Variable*>                        variables_;
                std::vector<Linear_constraint*>        constraints_;

                std::vector<FT>                result_;
                std::string                        error_message_;
                /// \endcond
        };

        //////////////////////////////////////////////////////////////////////////

        // implementation
#ifndef DOXYGEN_RUNNING
        template<typename FT>
        FT Bound<FT>::infinity_ = 1e20;                // values larger than this value are considered infinity

        template<typename FT>
        Bound<FT>::Bound(FT lb /* = -infinity() */, FT ub /* = +infinity() */)
                : lower_bound_(lb)
                , upper_bound_(ub)
        {
        }

        template<typename FT>
        FT Bound<FT>::infinity() {
                return infinity_;
        }

        template<typename FT>
        void Bound<FT>::set_bounds(FT lb, FT ub) {
                lower_bound_ = lb;
                upper_bound_ = ub;
        }

        template<typename FT>
        void Bound<FT>::get_bounds(FT& lb, FT& ub) const {
                lb = lower_bound_;
                ub = upper_bound_;
        }


        template<typename FT>
        Variable<FT>::Variable(
                Solver* solver,
                Variable_type type /* = CONTINUOUS */,
                FT lb /* = -infinity() */,
                FT ub /* = +infinity() */,
                const std::string& name /* = "" */,
                int idx /* = 0*/
        )
                : Solver_entry(solver, name, idx)
                , Bound(lb, ub)
                , variable_type_(type)
                , solution_value_(0.0)
        {
                if (type == BINARY)
                        Bound::set_bounds(0.0, 1.0);
        }


        template<typename FT>
        void Variable<FT>::set_variable_type(Variable_type type) {
                variable_type_ = type;
                if (type == BINARY)
                        Bound::set_bounds(0.0, 1.0);
        }


        template<typename FT>
        FT Variable<FT>::solution_value(bool rounded /* = false*/) const {
                if (rounded && variable_type_ != CONTINUOUS)
                        return std::round(solution_value_);
                else
                        return solution_value_;
        }

        //////////////////////////////////////////////////////////////////////////

        template<typename FT>
        Linear_expression<FT>::Linear_expression(Solver* solver, const std::string& name, int idx)
                : Solver_entry(solver, name, idx)
                , offset_(0.0)
        {
        }


        template<typename FT>
        void Linear_expression<FT>::add_coefficient(const Variable* var, FT coeff) {
                if (coefficients_.find(var) == coefficients_.end())
                        coefficients_[var] = coeff;
                else
                        coefficients_[var] += coeff;
        }


        template<typename FT>
        FT Linear_expression<FT>::get_coefficient(const Variable* var) const {
                typename std::unordered_map<const Variable*, FT>::const_iterator pos = coefficients_.find(var);
                if (pos != coefficients_.end())
                        return pos->second;
                else {
                        std::cerr << "linear expression does not own variable " << var->name() << " (" << var->index() << ")" << std::endl;
                        return 0.0;
                }
        }


        template<typename FT>
        FT Linear_expression<FT>::solution_value(bool rounded /* = false*/) const {
                FT solution = offset_;

                typename std::unordered_map<const Variable*, FT>::const_iterator it = coefficients_.begin();
                for (; it != coefficients_.end(); ++it) {
                        const Variable* var = it->first;
                        FT coeff = it->second;
                        solution += var->solution_value(rounded) * coeff;
                }
                return solution;
        }


        template<typename FT>
        void Linear_objective<FT>::clear() {
                Linear_expression::clear();
                Solver_entry::set_name("");
                Solver_entry::set_index(0);
        }


        template<typename FT>
        Linear_constraint<FT>::Linear_constraint(
                Solver* solver,
                FT lb /* = -infinity() */,
                FT ub /* = +infinity() */,
                const std::string& name/* = "" */,
                int idx /* = 0*/
        )
                : Linear_expression(solver, name, idx)
                , Bound(lb, ub)
        {
        }

        template<typename FT>
        Linear_objective<FT>::Linear_objective(Solver* solver, Sense sense)
                : Linear_expression(solver)
                , sense_(sense)
        {
        }

        template<typename FT>
        Mixed_integer_program_traits<FT>::Mixed_integer_program_traits() {
                // Intentionally set the objective to UNDEFINED, so it will allow me to warn
                // the user if he/she forgot to set the objective sense.
                objective_ = new Linear_objective(this, Linear_objective::UNDEFINED);
        }

        template<typename FT>
        Mixed_integer_program_traits<FT>::~Mixed_integer_program_traits() {
                clear();
                delete objective_;
        }

        template<typename FT>
        void Mixed_integer_program_traits<FT>::clear() {
                for (std::size_t i = 0; i < variables_.size(); ++i)
                        delete variables_[i];
                variables_.clear();

                for (std::size_t i = 0; i < constraints_.size(); ++i)
                        delete constraints_[i];
                constraints_.clear();

                objective_->clear();
        }

        namespace internal {
                /**
                * Converts an integer v to a string of specified 'width' by
                * filling with character 'fill'
                */
                template <class Int>
                inline std::string from_integer(Int v, int width, char fill) {
                        std::ostringstream string_stream;
                        string_stream << std::setfill(fill) << std::setw(width) << v;
                        return string_stream.str();
                }
        }

        template<typename FT>
        typename Mixed_integer_program_traits<FT>::Variable* Mixed_integer_program_traits<FT>::create_variable(
                Variable_type type /* = Variable::CONTINUOUS */,
                FT lb /* = -Variable::infinity() */,
                FT ub /* = +Variable::infinity() */,
                const std::string& name /* = "" */)
        {
                Variable* v = new Variable(this, type, lb, ub);

                std::size_t idx = variables_.size();
                v->set_index(static_cast<int>(idx));

                const std::string& fixed_name = name.empty() ? "x" + internal::from_integer(idx, 9, '0') : name;
                v->set_name(fixed_name);

                variables_.push_back(v);
                return v;
        }

        template<typename FT>
        std::vector<typename Mixed_integer_program_traits<FT>::Variable*> Mixed_integer_program_traits<FT>::create_variables(std::size_t n) {
                std::vector<Variable*> variables;
                for (std::size_t i = 0; i < n; ++i) {
                        Variable* v = create_variable();
                        variables.push_back(v);
                }
                return variables;
        }

        template<typename FT>
        typename Mixed_integer_program_traits<FT>::Linear_constraint* Mixed_integer_program_traits<FT>::create_constraint(
                FT lb /* = -Variable::infinity() */,
                FT ub /* = +Variable::infinity() */,
                const std::string& name /* = "" */)
        {
                Linear_constraint* c = new Linear_constraint(this, lb, ub);

                std::size_t idx = constraints_.size();
                c->set_index(static_cast<int>(idx));

                const std::string& fixed_name = name.empty() ? "c" + internal::from_integer(idx, 9, '0') : name;
                c->set_name(fixed_name);

                constraints_.push_back(c);
                return c;
        }

        template<typename FT>
        std::vector<typename Mixed_integer_program_traits<FT>::Linear_constraint*> Mixed_integer_program_traits<FT>::create_constraints(std::size_t n) {
                std::vector<Linear_constraint*> constraints;
                for (std::size_t i = 0; i < n; ++i) {
                        Linear_constraint* v = create_constraint();
                        constraints.push_back(v);
                }
                return constraints;
        }

        template<typename FT>
        typename Mixed_integer_program_traits<FT>::Linear_objective * Mixed_integer_program_traits<FT>::create_objective(Sense sense /* = Linear_objective ::MINIMIZE*/) {
                if (objective_)
                        delete objective_;

                objective_ = new Linear_objective(this, sense);
                return objective_;
        }


        template<typename FT>
        const typename Mixed_integer_program_traits<FT>::Linear_objective * Mixed_integer_program_traits<FT>::objective() const {
                if (!objective_)
                        std::cerr << "please call \'create_objective()\' to create an objective first" << std::endl;

                return objective_;
        }

        template<typename FT>
        typename Mixed_integer_program_traits<FT>::Linear_objective * Mixed_integer_program_traits<FT>::objective() {
                if (!objective_)
                        std::cerr << "please call \'create_objective()\' to create an objective first" << std::endl;

                return objective_;
        }

        template<typename FT>
        std::size_t Mixed_integer_program_traits<FT>::number_of_continuous_variables() const {
                std::size_t num_continuous_var = 0;
                for (std::size_t i = 0; i < variables_.size(); ++i) {
                        const Variable* v = variables_[i];
                        if (v->variable_type() == Variable::CONTINUOUS)
                                ++num_continuous_var;
                }
                return num_continuous_var;
        }

        template<typename FT>
        std::size_t Mixed_integer_program_traits<FT>::number_of_integer_variables() const {
                std::size_t num_iteger_var = 0;
                for (std::size_t i = 0; i < variables_.size(); ++i) {
                        const Variable* v = variables_[i];
                        if (v->variable_type() == Variable::INTEGER)
                                ++num_iteger_var;
                }
                return num_iteger_var;
        }

        template<typename FT>
        std::size_t Mixed_integer_program_traits<FT>::number_of_binary_variables() const {
                std::size_t num_binary_var = 0;
                for (std::size_t i = 0; i < variables_.size(); ++i) {
                        const Variable* v = variables_[i];
                        if (v->variable_type() == Variable::BINARY)
                                ++num_binary_var;
                }
                return num_binary_var;
        }

        // Returns true if all variables are continuous
        template<typename FT>
        bool Mixed_integer_program_traits<FT>::is_continuous() const {
                std::size_t num = number_of_continuous_variables();
                return (num > 0) && (num == variables_.size());
        }


        // Returns true if this is a mixed integer program
        template<typename FT>
        bool Mixed_integer_program_traits<FT>::is_mixed_integer_program() const {
                std::size_t num = number_of_continuous_variables();
                return (num > 0) && (num < variables_.size());
        }


        // Returns true if inter program
        template<typename FT>
        bool Mixed_integer_program_traits<FT>::is_integer_program() const {
                std::size_t num = number_of_integer_variables();
                return (num > 0) && (num == variables_.size());
        }


        // Returns true if binary program
        template<typename FT>
        bool Mixed_integer_program_traits<FT>::is_binary_program() const {
                std::size_t num = number_of_binary_variables();
                return (num > 0) && (num == variables_.size());
        }
#endif
} // namespace CGAL

#endif // CGAL_MIXED_INTEGER_PROGRAM_TRAITS_H
