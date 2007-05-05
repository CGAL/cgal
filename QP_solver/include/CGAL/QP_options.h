// Copyright (c) 1997-2007  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://gaertner@scm.gforge.inria.fr/svn/cgal/trunk/QP_solver/include/CGAL/QP_options.h $
// $Id: QP_options.h 38125 2007-04-14 15:39:55Z gaertner $
// 
//
// Author(s)     : Bernd Gaertner <gaertner@inf.ethz.ch>

#ifndef CGAL_QP_OPTIONS_H
#define CGAL_QP_OPTIONS_H
// this file defines a class for passing options to the linear and
// quadratic programming solvers

#include<string>
#include<CGAL/QP_solver/basic.h>

CGAL_BEGIN_NAMESPACE

enum Quadratic_program_pricing_strategy 
{ 
  QP_CHOOSE_DEFAULT,
  QP_DANTZIG, 
  QP_FILTERED_DANTZIG, 
  QP_PARTIAL_DANTZIG, 
  QP_PARTIAL_FILTERED_DANTZIG,
  QP_BLAND
};

class Quadratic_program_options 
{
public:
  // default constructor
  // -------------------
  Quadratic_program_options ()
    : verbosity_ (0), pricing_strategy_ (QP_CHOOSE_DEFAULT),
      validation_flag_ (false)
  {}

  // set/get verbosity
  // -----------------
  int get_verbosity () const 
  {
    return verbosity_;
  }
  void set_verbosity (int verbosity) 
  {
    CGAL_qpe_assertion ( 0 <= verbosity && verbosity <= 5);
    verbosity_ = verbosity;
  }

  // set/get pricing strategy
  // ------------------------
  Quadratic_program_pricing_strategy get_pricing_strategy() const
  {
    return pricing_strategy_;
  }

  void set_pricing_strategy
  (Quadratic_program_pricing_strategy pricing_strategy)
  {
    pricing_strategy_ = pricing_strategy;
  }

  // set/get validation flag
  // -----------------------
  bool get_auto_validation() const
  {
    return validation_flag_;
  }

  void set_auto_validation (bool validate)
  {
    validation_flag_ = validate;
  }

private:
  // verbosity
  // ---------
  //    0: silent
  //    1: short iteration summary (recommened for the user)
  // >= 2: output of internal details (not recommend for the user) 
  int verbosity_;   

  // pricing_strategy
  // ----------------
  Quadratic_program_pricing_strategy pricing_strategy_;

  // validation
  // ----------
  bool validation_flag_;
};

// output
// ------
std::ostream& operator<< (std::ostream& o, 
			  const Quadratic_program_options& options)
{
  o << "   Verbosity:        " << options.get_verbosity() << "\n";
  o << "   Pricing strategy: ";
  switch (options.get_pricing_strategy()) {
  case QP_CHOOSE_DEFAULT:
    o << "QP_CHOOSE_DEFAULT";
    break;
  case QP_DANTZIG:
    o << "QP_DANTZIG";
    break;
  case QP_FILTERED_DANTZIG:
    o << "QP_FILTERED_DANTZIG";
    break;
  case QP_PARTIAL_DANTZIG: 
    o << "QP_PARTIAL_FILTERED_DANTZIG";
    break;
  case QP_PARTIAL_FILTERED_DANTZIG:
    o << "QP_PARTIAL_FILTERED_DANTZIG";
    break;
  case QP_BLAND:
    o << "QP_BLAND";
  }
  o << "\n";
  o << "   Auto-validation:  " << options.get_auto_validation() << std::endl;
  return o;
}

CGAL_END_NAMESPACE

#endif // CGAL_QP_OPTIONS_H
