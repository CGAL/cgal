// Copyright (c) 2026 The CGAL Project
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Benjamin Phipps
//
// Compile tinygltf + stb implementation in exactly one translation unit.

#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <CGAL/IO/GLTF/tiny_gltf.h>
