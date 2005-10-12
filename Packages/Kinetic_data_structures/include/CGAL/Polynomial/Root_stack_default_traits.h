// Copyright (c) 2005  Stanford University (USA).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_ROOT_ENUMERATOR_DEFAULT_TRAITS_H
#define CGAL_POLYNOMIAL_ROOT_ENUMERATOR_DEFAULT_TRAITS_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Root_stack_traits_base.h>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE

template <class Polynomial>
class Root_stack_default_traits: public internal::Root_stack_traits_base<Polynomial>{
private:
  typedef internal::Root_stack_traits_base<Polynomial> Base;

public:
};

CGAL_POLYNOMIAL_END_NAMESPACE

#endif
