// Copyright (c) 2014  
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
//
//
// Author(s)     : Laurent Rineau

#include <CGAL/config.h>
#if CGAL_USE_GMP
#  include <CGAL/Gmpz.h>
#elif CGAL_USE_LEDA
#  include <CGAL/leda_integer.h>
#elif CGAL_USE_CORE
#  include <CGAL/CORE_BigInt.h>
#else
#  error CGAL is configured with none of GMP, LEDA and CORE. <CGAL/Exact_integer.h> cannot be used.
#endif

namespace CGAL {

/*!
\ingroup nt_cgal

`Exact_integer` is an exact integer number type.

It is a typedef of another number type. Its exact definition depends on
the availability the third-party libraries %GMP, %CORE, and %LEDA. %CGAL must
be configured with at least one of those libraries.

\cgalModels `EuclideanRing` 
\cgalModels `RealEmbeddable` 

*/
#if DOXYGEN_RUNNING

typedef unspecified_type Exact_integer;

#else // not DOXYGEN_RUNNING

#if CGAL_USE_GMP

typedef Gmpz Exact_integer;

#elif CGAL_USE_LEDA

typedef leda_integer Exact_integer;

#elif CGAL_USE_CORE

typedef CORE::BigInt Exact_integer;

#endif // CGAL_USE_CORE
#endif // not DOXYGEN_RUNNING

} /* end namespace CGAL */
