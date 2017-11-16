// Copyright (c) 2016  GeometryFactory SARL (France).
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s) : Andreas Fabri
//
// Warning: this file is generated, see include/CGAL/licence/README.md


#ifndef CGAL_LICENSE_RIDGES_3_H
#define CGAL_LICENSE_RIDGES_3_H

#include <CGAL/config.h>
#include <CGAL/license.h>




#ifdef CGAL_RIDGES_3_COMMERCIAL_LICENSE

#  if CGAL_RIDGES_3_COMMERCIAL_LICENSE < CGAL_RELEASE_DATE

#    if defined(CGAL_LICENSE_WARNING)

       CGAL_pragma_warning("Your commercial license for CGAL does not cover "
                           "this release of the Approximation of Ridges and Umbilics on Triangulated Surface Meshes package.")
#    endif

#    ifdef CGAL_LICENSE_ERROR
#      error "Your commercial license for CGAL does not cover this release \
of the Approximation of Ridges and Umbilics on Triangulated Surface Meshes package. \
You get this error, as you defined CGAL_LICENSE_ERROR."
#    endif // CGAL_LICENSE_ERROR

#  endif // CGAL_RIDGES_3_COMMERCIAL_LICENSE < CGAL_RELEASE_DATE

#else // no CGAL_RIDGES_3_COMMERCIAL_LICENSE

#  if defined(CGAL_LICENSE_WARNING)
     CGAL_pragma_warning("\nThe macro CGAL_RIDGES_3_COMMERCIAL_LICENSE is not defined."
                          "\nYou use the CGAL Approximation of Ridges and Umbilics on Triangulated Surface Meshes package under "
                          "the terms of the GPLv3+.")
#  endif // CGAL_LICENSE_WARNING

#  ifdef CGAL_LICENSE_ERROR
#    error "The macro CGAL_RIDGES_3_COMMERCIAL_LICENSE is not defined.\
 You use the CGAL Approximation of Ridges and Umbilics on Triangulated Surface Meshes package under the terms of \
the GPLv3+. You get this error, as you defined CGAL_LICENSE_ERROR."
#  endif // CGAL_LICENSE_ERROR

#endif // no CGAL_RIDGES_3_COMMERCIAL_LICENSE

#endif // CGAL_LICENSE_CHECK_RIDGES_3_H
