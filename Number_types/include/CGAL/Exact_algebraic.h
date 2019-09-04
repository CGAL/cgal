// Copyright (c) 2019  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Andreas Fabri

#include <CGAL/number_type_basic.h>

#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_integer.h>
#  include <CGAL/leda_rational.h>
#  include <CGAL/leda_real.h>
#endif
#ifdef CGAL_USE_CORE
#  include <CGAL/CORE_Expr.h>
#endif


namespace CGAL {

/*!
\ingroup nt_cgal

`Exact_algebraic` is an exact algebraic number type, constructible from `double`.

It is a typedef of another number type. Its exact definition depends on
the availability the third-party libraries %CORE, and %LEDA. %CGAL must
be configured with at least one of those libraries.

\cgalModels `FieldWithSqrt` 
\cgalModels `RealEmbeddable` 
\cgalModels `Fraction` 
\cgalModels `FromDoubleConstructible` 

*/
#if DOXYGEN_RUNNING

typedef unspecified_type Exact_algebraic;

#else // not DOXYGEN_RUNNING

#ifdef CGAL_USE_CORE
  typedef CORE::Expr Exact_algebraic;
#endif
  
#ifdef CGAL_USE_LEDA
typedef leda_real Exact_algebraic;
#endif

#endif
  
} /* end namespace CGAL */
