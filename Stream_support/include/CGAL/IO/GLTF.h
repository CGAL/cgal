// Copyright (c) 2025
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Dema Nasrawt and Aidan Pawlak

#ifndef CGAL_IO_GLTF_H
#define CGAL_IO_GLTF_H

// read_gltf.h must come first: it defines TINYGLTF_IMPLEMENTATION and
// includes tiny_gltf.h, so write_gltf.h can include it without re-defining
// the implementation bodies.
// read_gltf.h must come first: it defines TINYGLTF_IMPLEMENTATION and
// includes tiny_gltf.h (then undefines the macros), so write_gltf.h's
// second inclusion of tiny_gltf.h skips the implementation block.
#include <CGAL/IO/GLTF/read_gltf.h>
#include <CGAL/IO/GLTF/write_gltf.h>

#endif
