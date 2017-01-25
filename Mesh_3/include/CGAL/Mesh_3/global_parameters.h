// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_MESH_3_GLOBAL_PARAMETERS_H
#define CGAL_MESH_3_GLOBAL_PARAMETERS_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/config.h>
#include <CGAL/Mesh_3/config.h>

#ifdef BOOST_PARAMETER_MAX_ARITY
#  if (BOOST_PARAMETER_MAX_ARITY < 12)
#    error "BOOST_PARAMETER_MAX_ARITY must be at least 12 for CGAL::Mesh_3"
#  endif
#else
#  define  BOOST_PARAMETER_MAX_ARITY 12
#endif

#include <boost/parameter.hpp>


namespace CGAL {

namespace parameters {

template <typename T>
struct Base
{
  Base(T t) : t_(t) {}
  T operator()() const { return t_; }
private:
  T t_;
};
  
#define CGAL_MESH_BOOLEAN_PARAMETER(Class, function_true, function_false)     \
  struct Class : public Base<bool> { Class(bool b) : Base<bool>(b){} };       \
  inline Class function_true() { return Class(true); }                        \
  inline Class function_false() { return Class(false); }

#define CGAL_MESH_DOUBLE_PARAMETER(Class, function, precondition)             \
  struct Class : public Base<double>                                          \
  { Class(double d) : Base<double>(d) { precondition(d); } };                 \
  inline Class function(double d) { return Class(d); }

// see <CGAL/config.h>
CGAL_PRAGMA_DIAG_PUSH
// see <CGAL/Mesh_3/config.h>
CGAL_MESH_3_IGNORE_BOOST_PARAMETER_NAME_WARNINGS

BOOST_PARAMETER_NAME( c3t3 )
BOOST_PARAMETER_NAME( domain )
BOOST_PARAMETER_NAME( criteria )
  
BOOST_PARAMETER_NAME( (time_limit, tag) time_limit_ )
BOOST_PARAMETER_NAME( (sliver_bound, tag) sliver_bound_)
BOOST_PARAMETER_NAME( (sliver_criterion, tag) sliver_criterion_)
BOOST_PARAMETER_NAME( (perturbation_vector, tag) perturbation_vector_) 
BOOST_PARAMETER_NAME( (freeze_bound, tag) freeze_bound_)
BOOST_PARAMETER_NAME( (do_freeze, tag) do_freeze_)
BOOST_PARAMETER_NAME( (max_iteration_number, tag) max_iteration_number_ )
BOOST_PARAMETER_NAME( (convergence, tag) convergence_)

BOOST_PARAMETER_NAME( (dump_after_init_prefix, tag ) dump_after_init_prefix_)
BOOST_PARAMETER_NAME( (dump_after_refine_surface_prefix, tag ) dump_after_refine_surface_prefix_)
BOOST_PARAMETER_NAME( (dump_after_refine_prefix, tag ) dump_after_refine_prefix_)
BOOST_PARAMETER_NAME( (dump_after_glob_opt_prefix, tag ) dump_after_glob_opt_prefix_)
BOOST_PARAMETER_NAME( (dump_after_perturb_prefix, tag ) dump_after_perturb_prefix_)
BOOST_PARAMETER_NAME( (dump_after_exude_prefix, tag ) dump_after_exude_prefix_)
BOOST_PARAMETER_NAME( (number_of_initial_points, tag) number_of_initial_points_)

CGAL_PRAGMA_DIAG_POP
} // end namespace parameters





} //namespace CGAL

#endif // CGAL_MESH_3_GLOBAL_PARAMETERS_H
