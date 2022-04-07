// Copyright (c) 1999
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
// Author(s)     : Geert-Jan Giezeman and Sven Schoenherr

#include <CGAL/config.h>

#ifndef CGAL_ASSERTIONS_BEHAVIOUR_H
#define CGAL_ASSERTIONS_BEHAVIOUR_H

// workaround against the definition of EXIT in <opencv2/core/internal.hpp>
#ifdef EXIT
#  undef EXIT
#endif


namespace CGAL {

enum Failure_behaviour { ABORT, EXIT, EXIT_WITH_SUCCESS, CONTINUE,
                         THROW_EXCEPTION };

// failure handler declarations
// ==========================
// failure handler
// ---------------
typedef
    void
    (*Failure_function)(
        const char*, const char*, const char*, int, const char*);

CGAL_EXPORT
Failure_function
set_error_handler( Failure_function handler);

CGAL_EXPORT
Failure_function
set_warning_handler( Failure_function handler);

// failure behaviour handler
// -------------------------
CGAL_EXPORT
Failure_behaviour
set_error_behaviour(Failure_behaviour eb);

CGAL_EXPORT
Failure_behaviour
set_warning_behaviour(Failure_behaviour eb);

} //namespace CGAL

#endif // CGAL_ASSERTIONS_BEHAVIOUR_H
