// Copyright (c) 2010-2011  GeometryFactory Sarl (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Laurent Rineau

extern "C" { 
#include <fenv.h>
}

namespace CGAL {

// replacement for C99

inline int
feclearexcept(int exceptions) {
    int mxcsr;
    asm volatile("stmxcsr %0" : "=m" (mxcsr) );
    mxcsr &= ~exceptions;
    asm volatile("ldmxcsr %0" : : "m" (mxcsr) );
    return 0;
}

inline int
fetestexcept(int exceptions) {
    int status;
    asm volatile("stmxcsr %0" : "=m" (status) );
    return status & exceptions;
}

} // end namespace CGAL
