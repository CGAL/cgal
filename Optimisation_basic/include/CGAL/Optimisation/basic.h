// Copyright (c) 1997-2001
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
// Author(s)     : Sven Schoenherr <sven@inf.ethz.ch>

#ifndef CGAL_OPTIMISATION_BASIC_H
#define CGAL_OPTIMISATION_BASIC_H

// includes
#include <CGAL/basic.h>
#include <CGAL/Optimisation/assertions.h>
#include <CGAL/Optimisation/debug.h>
#include <CGAL/IO/Verbose_ostream.h>

namespace CGAL {

// Function declarations
// =====================

// is_valid failure function
// -------------------------
inline
bool
_optimisation_is_valid_fail( CGAL::Verbose_ostream& verr,
                             const char*            message)
{
    verr << "FAILED." << std::endl;
    verr << "  --> " << message << std::endl;
    verr << "  object is NOT valid!" << std::endl;
    return false;
}
} //namespace CGAL

#endif // CGAL_OPTIMISATION_BASIC_H
