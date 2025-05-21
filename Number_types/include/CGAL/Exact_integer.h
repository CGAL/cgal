// Copyright (c) 2014
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_EXACT_INTEGER_H
#define CGAL_EXACT_INTEGER_H

#include <CGAL/Number_types/internal/Exact_type_selector.h>

namespace CGAL {

/*!
\ingroup nt_cgal

`Exact_integer` is an exact integer number type.

It is a typedef of another number type. Its exact definition depends on
the availability the third-party libraries \gmp, \core, and \leda. \cgal must
be configured with at least one of those libraries.

\cgalModels{EuclideanRing,RealEmbeddable}

*/
#if DOXYGEN_RUNNING

typedef unspecified_type Exact_integer;

#else // not DOXYGEN_RUNNING

using Exact_integer = internal::Exact_NT_backend<internal::Default_exact_nt_backend>::Integer;

#endif // not DOXYGEN_RUNNING

} /* end namespace CGAL */

#endif // CGAL_EXACT_INTEGER_H
