// Copyright (c) 2006-2010 Inria Lorraine (France). All rights reserved.
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_BASIC_H
#define CGAL_RS_BASIC_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>

#ifndef CGAL_RS_VERB
#ifdef CGAL_RS_DEBUG
#define CGAL_RS_VERB    2
#else
#define CGAL_RS_VERB    0
#endif
#endif

// the default precision of RS to calculate a root (precision is 2^n)
#ifndef CGAL_RS_DEF_PREC
#define CGAL_RS_DEF_PREC        0
#endif

// the minimum, used when calculating a sign
#ifndef CGAL_RS_MIN_PREC
#define CGAL_RS_MIN_PREC        0
#endif

#if ( defined CGAL_HAS_THREADS && !defined CGAL_RS_NO_TLS )
#  ifdef _MSC_VER
#    ifdef _WINDLL
#      error "Can't build CGAL_RS as thread safe."
#      define CGALRS_THREAD_ATTR
#    else
#      define CGALRS_THREAD_ATTR __declspec(thread)
#    endif
#  else
#    define CGALRS_THREAD_ATTR __thread
#  endif
#else
#  define CGALRS_THREAD_ATTR
#endif

namespace RS{

enum rs_sign{
    RS_NEGATIVE=        CGAL::NEGATIVE,
    RS_ZERO=            CGAL::ZERO,
    RS_POSITIVE=        CGAL::POSITIVE,
    RS_UNKNOWN
};

// this function must only be called when s is not RS_UNKNOWN
inline CGAL::Sign convert_rs_sign(rs_sign s){
        CGAL_precondition(s!=RS_UNKNOWN);
        switch(s){
                case RS_NEGATIVE:return CGAL::NEGATIVE;break;
                case RS_POSITIVE:return CGAL::POSITIVE;break;
                default:return CGAL::ZERO;
        }
}

} // namespace RS

#endif  // CGAL_RS_BASIC_H
