// Copyright (c) 2019 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : André Nusser <anusser@mpi-inf.mpg.de>
//                 Marvin Künnemann <marvin@mpi-inf.mpg.de>
//                 Karl Bringmann <kbringma@mpi-inf.mpg.de>
//
// =============================================================================

#pragma once
#include <CGAL/license/Polyline_distance.h>#include <CGAL/license/Polyline_distance.h>
#include <iostream>

// helper macro

#define CENTER_NOP \
    do {           \
    } while (0)

// printing macros

#ifdef NDEBUG
#define DEBUG(x) CENTER_NOP
#else
#define DEBUG(x)                     \
    do {                             \
        std::cout << x << std::endl; \
    } while (0)
#endif

#if defined(NVERBOSE) && defined(NDEBUG)
#define PRINT(x) CENTER_NOP
#else
#define PRINT(x)                     \
    do {                             \
        std::cout << x << std::endl; \
    } while (0)
#endif

#define ERROR(x)                                  \
    do {                                          \
        std::cerr << "Error: " << x << std::endl; \
        std::exit(EXIT_FAILURE);                  \
    } while (0)

#include <cassert>
