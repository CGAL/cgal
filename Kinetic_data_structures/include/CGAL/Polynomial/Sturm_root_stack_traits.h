// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_STURM_ROOT_ENUMERATOR_TRAITS_H
#define CGAL_POLYNOMIAL_STURM_ROOT_ENUMERATOR_TRAITS_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/Root_stack_default_traits.h>
#include <CGAL/Polynomial/internal/Rational/Euclidean_Sturm_sequence.h>
#include <CGAL/Polynomial/internal/Rational/Monic_Sturm_sequence.h>
#include <CGAL/Polynomial/internal/Rational/Primitive_part_Sturm_sequence.h>
#include <CGAL/Polynomial/internal/Rational/Reduced_Sturm_sequence.h>
#include <CGAL/Polynomial/internal/Rational/Subresultant_Sturm_sequence.h>
#include <CGAL/Polynomial/internal/Rational/Sturm_root_counter.h>

namespace CGAL { namespace POLYNOMIAL {

template<class Polynomial>
class Sturm_root_stack_traits
  : public Root_stack_default_traits<Polynomial>
{
private:
  typedef Root_stack_default_traits<Polynomial>  Base;
  typedef Sturm_root_stack_traits<Polynomial>    Self;

public:
  typedef internal::Sturm_root_counter<Self> Root_count;
  typedef typename Base::Function            Function;


  typedef internal::Monic_Sturm_sequence<Self> Sturm_sequence;
 
  Sturm_sequence Sturm_sequence_object(const Function &f,
				       const Function &g) const
  {
    return Sturm_sequence(f, g, *this);
  }

  typedef internal::Standard_sequence<Sturm_sequence> Standard_sequence;
  Standard_sequence standard_sequence_object(const Function &f) const
  {
    return Standard_sequence(f, *this);
  }

  typedef internal::Sign_Sturm_sequence<Sturm_sequence> Sign_Sturm_sequence;
  Sign_Sturm_sequence sign_Sturm_sequence_object(const Function &f, const
						 Function &g) const
  {
    return Sign_Sturm_sequence(f, g, *this);
  }

  typedef internal::Root_bound_evaluator<Function> Root_bound;
  Root_bound root_bound_object(bool b = true) const
  {
    return Root_bound(b);
  }
};

} } //namespace CGAL::POLYNOMIAL
#endif                                            // CGAL_POLYNOMIAL_STURM_ROOT_ENUMERATOR_TRAITS_H
