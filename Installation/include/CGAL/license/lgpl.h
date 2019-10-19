// Copyright (c) 2016  GeometryFactory Sarl (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Andreas Fabri, Laurent Rineau

#ifndef CGAL_LICENSE_LGPL_H
#define CGAL_LICENSE_LGPL_H

#include <CGAL/config.h>
#include <CGAL/license.h>



#ifdef CGAL_LGPL_COMMERCIAL_LICENSE

#  if CGAL_LGPL_COMMERCIAL_LICENSE < CGAL_RELEASE_DATE

#    if defined(CGAL_LICENSE_WARNING)
       CGAL_pragma_warning("Your commercial license for CGAL does not cover "
                           "this release of CGAL LGPL packages.")
#    endif

#    ifdef CGAL_LICENSE_ERROR
#      error "Your commercial license for CGAL does not cover this release \
of CGAL LGPL packages. You get this error, as you defined CGAL_LICENSE_ERROR."
#    endif // CGAL_LICENSE_ERROR

#  endif // CGAL_LGPL_COMMERCIAL_LICENSE < CGAL_RELEASE_DATE

#else // not CGAL_LGPL_COMMERCIAL_LICENSE

#  if defined(CGAL_LGPL_LICENSE_WARNING)
     CGAL_pragma_warning("You use the CGAL LGPL packages under the terms "
                         "of the LGPLv3+.")
#  endif // CGAL_LGPL_LICENSE_WARNING

#  ifdef CGAL_LGPL_LICENSE_ERROR
#    error "You use the CGAL AABB Tree package under the terms of the GPLv3+. \
You get this error, as you defined CGAL_LGPL_LICENSE_ERROR."
#  endif // CGAL_LGPL_LICENSE_ERROR

#endif // not CGAL_LGPL_COMMERCIAL_LICENSE

#endif // CGAL_LICENSE_CHECK_LGPL_H
