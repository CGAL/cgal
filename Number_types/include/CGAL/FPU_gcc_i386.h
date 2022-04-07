// Copyright (c) 2010-2011  GeometryFactory Sarl (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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
