// Copyright (c) 2018  Liangliang Nan. All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Liangliang Nan

#ifndef CGAL_GLPK_MIXED_INTEGER_PROGRAM_TRAITS_H
#define CGAL_GLPK_MIXED_INTEGER_PROGRAM_TRAITS_H

#include <CGAL/Mixed_integer_program_traits.h>

#if defined(CGAL_USE_GLPK) || defined(DOXYGEN_RUNNING)

#include <glpk.h>

#include <iostream>
#include <vector>

namespace CGAL {
namespace internal {

// Infers "bound type" (required by GLPK) from the bounds values
template <typename FT>
int bound_type(FT lb, FT ub)
{
  typedef CGAL::Variable<FT> Variable;

  if (lb <= -Variable::infinity() && ub >= Variable::infinity()) {
    return GLP_FR; // free (unbounded) variable
  } else if (lb > -Variable::infinity() && ub >= Variable::infinity()) {
    return GLP_LO; // variable with lower bound
  } else if (lb <= -Variable::infinity() && ub < Variable::infinity()) {
    return GLP_UP; // variable with upper bound
  } else {
    // lb > -Variable::infinity() && ub < Variable::infinity()
    if (lb == ub)
      return GLP_FX; // fixed variable
    else
      return GLP_DB; // double-bounded variable
  }
}

} // namespace internal

/// \ingroup PkgSolverInterfaceRef
///
/// This class provides an interface for formulating and solving
/// constrained or unconstrained mixed integer programs using
/// \ref thirdpartyGLPK, which must be available on the system.
///
/// \note For better performance, please consider using
///      `CGAL::SCIP_mixed_integer_program_traits`, or derive a new
///      model from `CGAL::Mixed_integer_program_traits`.
///
/// \cgalModels `MixedIntegerProgramTraits`
///
/// \sa `SCIP_mixed_integer_program_traits`
template <typename FT>
class GLPK_mixed_integer_program_traits
  : public Mixed_integer_program_traits<FT>
{
  /// \cond SKIP_IN_MANUAL
public:
  typedef CGAL::Mixed_integer_program_traits<FT>          Base_class;
  typedef typename Base_class::Variable                   Variable;
  typedef typename Base_class::Linear_constraint          Linear_constraint;
  typedef typename Base_class::Linear_objective           Linear_objective;
  typedef typename Linear_objective::Sense                Sense;
  typedef typename Variable::Variable_type                Variable_type;

public:
  /// Solves the program. Returns `false` if fails.
  virtual bool solve()
  {
    Base_class::error_message_.clear();

    glp_prob* lp = glp_create_prob();
    if (!lp) {
      Base_class::error_message_ = "failed creating GLPK program";
      return false;
    }

    std::size_t num_integer_variables = 0;

    // Creates variables

    // This "static_cast<>" suppresses many annoying warnings: "conversion from 'size_t' to 'int', possible loss of data"
    int num_variables = static_cast<int>(Base_class::variables_.size());

    glp_add_cols(lp, num_variables);
    for (int i = 0; i < num_variables; ++i) {
      const Variable* var = Base_class::variables_[i];
      glp_set_col_name(lp, i + 1, var->name().c_str());

      if (var->variable_type() == Variable::INTEGER) {
        glp_set_col_kind(lp, i + 1, GLP_IV);  // glpk uses 1-based arrays
        ++num_integer_variables;
      } else if (var->variable_type() == Variable::BINARY) {
        glp_set_col_kind(lp, i + 1, GLP_BV);  // glpk uses 1-based arrays
        ++num_integer_variables;
      } else {
        glp_set_col_kind(lp, i + 1, GLP_CV);  // Continuous variable
      }

      FT lb, ub;
      var->get_bounds(lb, ub);

      int type = internal::bound_type(lb, ub);
      glp_set_col_bnds(lp, i + 1, type, lb, ub);
    }

    // Adds constraints

    // This "static_cast<>" suppresses many annoying warnings: "conversion from 'size_t' to 'int', possible loss of data"
    int num_constraints = static_cast<int>(Base_class::constraints_.size());
    glp_add_rows(lp, num_constraints);

    for (int i = 0; i < num_constraints; ++i) {
      const Linear_constraint* c = Base_class::constraints_[i];
      const std::unordered_map<const Variable*, FT>& coeffs = c->coefficients();
      typename std::unordered_map<const Variable*, FT>::const_iterator cur = coeffs.begin();

      std::vector<int>  indices(coeffs.size() + 1, 0); // glpk uses 1-based arrays
      std::vector<FT> coefficients(coeffs.size() + 1, 0.0); // glpk uses 1-based arrays
      std::size_t idx = 1; // glpk uses 1-based arrays
      for (; cur != coeffs.end(); ++cur) {
        int var_idx = cur->first->index();
        FT coeff = cur->second;

        indices[idx] = var_idx + 1;   // glpk uses 1-based arrays
        coefficients[idx] = coeff;
        ++idx;
      }

      glp_set_mat_row(lp, i + 1, static_cast<int>(coeffs.size()), indices.data(), coefficients.data());

      FT lb, ub;
      c->get_bounds(lb, ub);

      int type = internal::bound_type(lb, ub);
      glp_set_row_bnds(lp, i + 1, type, lb, ub);

      glp_set_row_name(lp, i + 1, c->name().c_str());
    }

    // Sets objective

    // Determines the coefficient of each variable in the objective function
    const std::unordered_map<const Variable*, FT>& obj_coeffs = Base_class::objective_->coefficients();
    typename std::unordered_map<const Variable*, FT>::const_iterator cur = obj_coeffs.begin();
    for (; cur != obj_coeffs.end(); ++cur) {
      int var_idx = cur->first->index();
      FT coeff = cur->second;
      glp_set_obj_coef(lp, var_idx + 1, coeff); // glpk uses 1-based arrays
    }

    // Sets objective function sense
    bool minimize = (Base_class::objective_->sense() == Linear_objective::MINIMIZE);
    glp_set_obj_dir(lp, minimize ? GLP_MIN : GLP_MAX);
    int msg_level = GLP_MSG_ERR;
    int status = -1;
    if (num_integer_variables == 0) { // Continuous problem
      glp_smcp parm;
      glp_init_smcp(&parm);
      parm.msg_lev = msg_level;
      status = glp_simplex(lp, &parm);
    }
    else { // Solves as MIP problem
      glp_iocp parm;
      glp_init_iocp(&parm);
      parm.msg_lev = msg_level;
      parm.presolve = GLP_ON;
      // The routine glp_intopt is a driver to the MIP solver based on the branch-and-cut method,
      // which is a hybrid of branch-and-bound and cutting plane methods.
      status = glp_intopt(lp, &parm);
    }

    switch (status) {
      case 0: {
        if (num_integer_variables == 0) { // continuous problem
          Base_class::result_.resize(num_variables);
          for (int i = 0; i < num_variables; ++i) {
            FT x = glp_get_col_prim(lp, i + 1);  // glpk uses 1-based arrays
            Variable* v = Base_class::variables_[i];
            v->set_solution_value(x);
            Base_class::result_[i] = x;
          }
        }
        else { // MIP problem
          Base_class::result_.resize(num_variables);
          for (int i = 0; i < num_variables; ++i) {
            FT x = glp_mip_col_val(lp, i + 1); // glpk uses 1-based arrays

            Variable* v = Base_class::variables_[i];
            if (v->variable_type() != Variable::CONTINUOUS)
              x = static_cast<int>(std::round(x));

            v->set_solution_value(x);
            Base_class::result_[i] = x;
          }
        }
        break;
      }

      case GLP_EBOUND:
        Base_class::error_message_ =
            "Unable to start the search, because some FT-bounded variables have incorrect"
            "bounds or some integer variables have non - integer(fractional) bounds.";
        break;

      case GLP_EROOT:
        Base_class::error_message_ =
            "Unable to start the search, because optimal basis for initial LP relaxation is not"
            "provided. (This code may appear only if the presolver is disabled.)";
        break;

      case GLP_ENOPFS:
        Base_class::error_message_ =
            "Unable to start the search, because LP relaxation of the MIP problem instance has"
            "no primal feasible solution. (This code may appear only if the presolver is enabled.)";
        break;

      case GLP_ENODFS:
        Base_class::error_message_ =
            "Unable to start the search, because LP relaxation of the MIP problem instance has"
            "no dual feasible solution.In other word, this code means that if the LP relaxation"
            "has at least one primal feasible solution, its optimal solution is unbounded, so if the"
            "MIP problem has at least one integer feasible solution, its(integer) optimal solution"
            "is also unbounded. (This code may appear only if the presolver is enabled.)";
        break;

      case GLP_EFAIL:
        Base_class::error_message_ =
            "The search was prematurely terminated due to the solver failure.";
        break;

      case GLP_EMIPGAP:
        Base_class::error_message_ =
            "The search was prematurely terminated, because the relative mip gap tolerance has been reached.";
        break;

      case GLP_ETMLIM:
        Base_class::error_message_ =
            "The search was prematurely terminated, because the time limit has been exceeded.";
        break;

      case GLP_ESTOP:
        Base_class::error_message_ =
            "The search was prematurely terminated by application. (This code may appear only"
            "if the advanced solver interface is used.)";
        break;

      default:
        Base_class::error_message_ =
            "optimization was stopped with status code " + std::to_string(status);
        break;
    }

    glp_delete_prob(lp);

    return (status == 0);
  }
  /// \endcond
};

} // namespace CGAL

#endif // CGAL_USE_GLPK or DOXYGEN_RUNNING

#endif // CGAL_GLPK_MIXED_INTEGER_PROGRAM_TRAITS_H
