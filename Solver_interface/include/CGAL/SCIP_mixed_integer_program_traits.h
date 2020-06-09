// Copyright (c) 2018  Liangliang Nan. All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Liangliang Nan

#ifndef CGAL_SCIP_MIXED_INTEGER_PROGRAM_TRAITS_H
#define CGAL_SCIP_MIXED_INTEGER_PROGRAM_TRAITS_H

#include <CGAL/Mixed_integer_program_traits.h>

#if defined(CGAL_USE_SCIP) || defined(DOXYGEN_RUNNING)

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace CGAL {

/// \ingroup PkgSolverInterfaceRef
///
/// This class provides an interface for formulating and solving
/// constrained or unconstrained mixed integer programs using
/// \ref thirdpartySCIP (which must be available on the system).
///
/// \cgalModels `MixedIntegerProgramTraits`
///
/// \sa `GLPK_mixed_integer_program_traits`
template <typename FT>
class SCIP_mixed_integer_program_traits
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

    Scip* scip = 0;
    SCIP_CALL(SCIPcreate(&scip));
    SCIP_CALL(SCIPincludeDefaultPlugins(scip));

    // Disables scip output to stdout
    SCIPmessagehdlrSetQuiet(SCIPgetMessagehdlr(scip), TRUE);

    // Uses wall clock time because getting CPU user seconds
    // involves calling times() which is very expensive
    SCIP_CALL(SCIPsetIntParam(scip, "timing/clocktype", SCIP_CLOCKTYPE_WALL));

    // Creates empty problem
    SCIP_CALL(SCIPcreateProbBasic(scip, "Solver_interface"));

    // Creates variables
    std::vector<SCIP_VAR*> scip_variables;
    for (std::size_t i = 0; i < Base_class::variables_.size(); ++i) {
      const Variable* var = Base_class::variables_[i];
      SCIP_VAR* v = 0;

      double lb, ub;
      var->get_bounds(lb, ub);

      switch (var->variable_type())
      {
        case Variable::CONTINUOUS:
          SCIP_CALL(SCIPcreateVar(scip, &v, var->name().c_str(), lb, ub, 0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, 0, 0, 0, 0, 0));
          break;
        case Variable::INTEGER:
          SCIP_CALL(SCIPcreateVar(scip, &v, var->name().c_str(), lb, ub, 0.0, SCIP_VARTYPE_INTEGER, TRUE, FALSE, 0, 0, 0, 0, 0));
          break;
        case Variable::BINARY:
          SCIP_CALL(SCIPcreateVar(scip, &v, var->name().c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, 0, 0, 0, 0, 0));
          break;
      }
      // Adds the SCIP_VAR object to the scip problem
      SCIP_CALL(SCIPaddVar(scip, v));

      // Stores the SCIP_VAR pointer for later access
      scip_variables.push_back(v);
    }

    // Adds constraints

    std::vector<SCIP_CONS*> scip_constraints;
    for (std::size_t i = 0; i < Base_class::constraints_.size(); ++i) {
      const Linear_constraint* c = Base_class::constraints_[i];
      const std::unordered_map<const Variable*, double>& coeffs = c->coefficients();
      typename std::unordered_map<const Variable*, double>::const_iterator cur = coeffs.begin();

      std::vector<SCIP_VAR*>  cstr_variables(coeffs.size());
      std::vector<double>    cstr_values(coeffs.size());
      std::size_t idx = 0;
      for (; cur != coeffs.end(); ++cur) {
        std::size_t var_idx = cur->first->index();
        double coeff = cur->second;
        cstr_variables[idx] = scip_variables[var_idx];
        cstr_values[idx] = coeff;
        ++idx;
      }

      // Creates SCIP_CONS object
      SCIP_CONS* cons = 0;
      const std::string& name = c->name();

      double lb, ub;
      c->get_bounds(lb, ub);

      SCIP_CALL(SCIPcreateConsLinear(scip, &cons, name.c_str(), int(coeffs.size()), cstr_variables.data(), cstr_values.data(), lb, ub, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

      // Adds the constraint to scip
      SCIP_CALL(SCIPaddCons(scip, cons));

      // Stores the constraint for later on
      scip_constraints.push_back(cons);
    }

    // Sets objective

    // Determines the coefficient of each variable in the objective function
    const std::unordered_map<const Variable*, double>& obj_coeffs = Base_class::objective_->coefficients();
    typename std::unordered_map<const Variable*, double>::const_iterator cur = obj_coeffs.begin();
    for (; cur != obj_coeffs.end(); ++cur) {
      const Variable* var = cur->first;
      double coeff = cur->second;
      SCIP_CALL(SCIPchgVarObj(scip, scip_variables[var->index()], coeff));
    }

    // Sets the objective sense
    bool minimize = (Base_class::objective_->sense() == Linear_objective::MINIMIZE);
    SCIP_CALL(SCIPsetObjsense(scip, minimize ? SCIP_OBJSENSE_MINIMIZE : SCIP_OBJSENSE_MAXIMIZE));

    // Turns presolve on (it's the SCIP default).
    bool presolve = true;
    if (presolve)
      SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrounds", -1)); // maximal number of presolving rounds (-1: unlimited, 0: off)
    else
      SCIP_CALL(SCIPsetIntParam(scip, "presolving/maxrounds", 0));  // disable presolve

    bool status = false;
    // This tells scip to start the solution process
    if (SCIPsolve(scip) == SCIP_OKAY) {
      // Gets the best found solution from scip
      SCIP_SOL* sol = SCIPgetBestSol(scip);
      if (sol) {
        // If optimal or feasible solution is found.
        Base_class::result_.resize(Base_class::variables_.size());
        for (std::size_t i = 0; i < Base_class::variables_.size(); ++i) {
          FT x = SCIPgetSolVal(scip, sol, scip_variables[i]);

          Variable* v = Base_class::variables_[i];
          if (v->variable_type() != Variable::CONTINUOUS)
            x = static_cast<int>(std::round(x));

          v->set_solution_value(x);
          Base_class::result_[i] = x;
        }
        status = true;
      }
    }

    // Reports the status: optimal, infeasible, etc.
    SCIP_STATUS scip_status = SCIPgetStatus(scip);
    switch (scip_status) {
      case SCIP_STATUS_OPTIMAL:
        // Provides info only if fails.
        break;
      case SCIP_STATUS_GAPLIMIT:
        // To be consistent with the other solvers.
        // Provides info only if fails.
        break;
      case SCIP_STATUS_INFEASIBLE:
        Base_class::error_message_ = "model was infeasible";
        break;
      case SCIP_STATUS_UNBOUNDED:
        Base_class::error_message_ = "model was unbounded";
        break;
      case SCIP_STATUS_INFORUNBD:
        Base_class::error_message_ = "model was either infeasible or unbounded";
        break;
      case SCIP_STATUS_TIMELIMIT:
        Base_class::error_message_ = "aborted due to time limit";
        break;
      default:
        Base_class::error_message_ = "aborted with status: " + std::to_string(scip_status);
        break;
    }

    SCIP_CALL(SCIPresetParams(scip));

    // Since the SCIPcreateVar captures all variables, we have to release them now
    for (std::size_t i = 0; i < scip_variables.size(); ++i)
      SCIP_CALL(SCIPreleaseVar(scip, &scip_variables[i]));
    scip_variables.clear();

    // The same for the constraints
    for (std::size_t i = 0; i < scip_constraints.size(); ++i)
      SCIP_CALL(SCIPreleaseCons(scip, &scip_constraints[i]));
    scip_constraints.clear();

    // After releasing all vars and cons we can free the scip problem.
    // Remember this has always to be the last call to scip
    SCIP_CALL(SCIPfree(&scip));

    return status;
  }
  /// \endcond
};

} // namespace CGAL

#endif // CGAL_USE_SCIP or DOXYGEN_RUNNING

#endif // CGAL_SCIP_MIXED_INTEGER_PROGRAM_TRAITS_H
