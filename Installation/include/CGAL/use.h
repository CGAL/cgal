// Copyright (c) 2007   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_USE_H
#define CGAL_USE_H

namespace CGAL { namespace internal {

template <typename ...T> inline
void use(T&&...) {}

template<typename> void use_type() {}

} }

/// CGAL_USE() is a macro which aims at removing "variable is unused" warnings.
#define CGAL_USE(...) ::CGAL::internal::use(__VA_ARGS__)

/// CGAL_USE_TYPE() is a macro which aims at removing "typedef locally
/// defined but not used" warnings.
#define CGAL_USE_TYPE(T) ::CGAL::internal::use_type<T>()

#endif // CGAL_USE_H
