// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_VECTOR_DETVEC_SMALL_H
#define CGAL_VECTOR_DETVEC_SMALL_H
#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/Dimension.h>
#include <CGAL/determinant_of_vectors.h>

#define CGAL_ALLOWED_INCLUSION 1

#define CGAL_CLASS Add_determinant_of_vectors_small_dim
#define CGAL_TAG Has_determinant_of_vectors_tag
#define CGAL_FUNC determinant_of_vectors
#define CGAL_SIGN_FUNC sign_of_determinant_of_vectors
#define CGAL_SHIFT 0

#include <CGAL/NewKernel_d/Vector/determinant_of_vectors_small_dim_internal.h>

#undef CGAL_CLASS
#undef CGAL_TAG
#undef CGAL_FUNC
#undef CGAL_SIGN_FUNC
#undef CGAL_SHIFT

#define CGAL_CLASS Add_determinant_of_vectors_omit_last_small_dim
#define CGAL_TAG Has_determinant_of_vectors_omit_last_tag
#define CGAL_FUNC determinant_of_vectors_omit_last
#define CGAL_SIGN_FUNC sign_of_determinant_of_vectors_omit_last
#define CGAL_SHIFT 1

#include <CGAL/NewKernel_d/Vector/determinant_of_vectors_small_dim_internal.h>

#undef CGAL_CLASS
#undef CGAL_TAG
#undef CGAL_FUNC
#undef CGAL_SIGN_FUNC
#undef CGAL_SHIFT

#undef CGAL_ALLOWED_INCLUSION

#endif
