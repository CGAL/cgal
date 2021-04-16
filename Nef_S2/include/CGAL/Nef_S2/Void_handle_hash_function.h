// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel  <seel@mpi-sb.mpg.de>

#ifndef CGAL_VOID_HANDLE_HASH_FUNCTION_H
#define CGAL_VOID_HANDLE_HASH_FUNCTION_H

#include <CGAL/license/Nef_S2.h>

namespace CGAL {

struct Void_handle_hash_function {
    std::size_t operator() (void* h) const {
        return std::size_t(h);
    }
};

} //namespace CGAL
#endif //CGAL_VOID_HANDLE_HASH_FUNCTION_H
