// Copyright (c) 2007   INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_USE_H
#define CGAL_USE_H

namespace CGAL { namespace internal {

template < typename T > inline
void use(const T&) {}

template<typename> void use_type() {}

} }

/// CGAL_USE() is a macro which aims at removing "variable is unused" warnings.
#define CGAL_USE(x) ::CGAL::internal::use(x)

/// CGAL_USE_TYPE() is a macro which aims at removing "typedef locally
/// defined but not used" warnings.
#define CGAL_USE_TYPE(T) ::CGAL::internal::use_type<T>()

#endif // CGAL_USE_H
