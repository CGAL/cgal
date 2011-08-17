// Copyright (c) 2010-2011  GeometryFactory Sarl (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
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
//
//
// Author(s)     : Laurent Rineau

extern "C" { 
#include <fenv.h>
}

namespace CGAL {

// brute-force replacement for C99 (which does not require an inline-function)
inline int
feclearexcept(int exceptions) {
    // TODO: clear only given exceptions
    asm volatile("fnclex");
    return 0;
}

inline int
fetestexcept(int exceptions) {
    int status;
    asm volatile("fnstsw %0" : "=m" (status));
    return status & exceptions;
}

} // end namespace CGAL
