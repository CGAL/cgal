// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.

// See the file LICENSE.LGPL distributed with CGAL.
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_GMPFR_MAKE_UNIQUE_H
#define CGAL_RS_GMPFR_MAKE_UNIQUE_H

#include <CGAL/Gmpfr.h>

// Make sure a number does not share references. If it does, copy it.
#ifdef CGAL_GMPFR_NO_REFCOUNT
#  define CGAL_RS_GMPFR_MAKE_UNIQUE(_number,_tempvar) {};
#else // CGAL_GMPFR_NO_REFCOUNT
#  define CGAL_RS_GMPFR_MAKE_UNIQUE(_number,_tempvar) \
        if(!_number.is_unique()){ \
                Gmpfr _tempvar(0,_number.get_precision()); \
                mpfr_set(_tempvar.fr(),_number.fr(),GMP_RNDN); \
                _number=_tempvar; \
                CGAL_assertion_code(_tempvar=Gmpfr();) \
                CGAL_assertion(_number.is_unique()); \
        }
#endif // CGAL_GMPFR_NO_REFCOUNT

#endif // CGAL_RS_GMPFR_MAKE_UNIQUE_H
