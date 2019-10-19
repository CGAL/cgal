// Copyright (c) 2006-2013 INRIA Nancy-Grand Est (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
