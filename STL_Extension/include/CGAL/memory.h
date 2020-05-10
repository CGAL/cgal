// Copyright (c) 1999, 2002, 2003
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
// Author(s)     : Michael Seel, Sylvain Pion

#ifndef CGAL_MEMORY_H
#define CGAL_MEMORY_H

#include <memory>

// CGAL_ALLOCATOR(t) defines the default allocator used by CGAL.
// CGAL_MEMORY(t) overloads the new and delete operators for a given class.

// When LEDA is there, the user could define these macros as
// leda_allocator< T >  and  LEDA_MEMORY(T) for example.

// For debugging with GCC, the following allocator can be useful :
// std::__allocator<T, std::__debug_alloc<std::__malloc_alloc_template<0> > >

#ifndef CGAL_ALLOCATOR
#  define CGAL_ALLOCATOR(T) std::allocator< T >
#endif

#ifndef CGAL_MEMORY
#  define CGAL_MEMORY(T)
#endif

#endif // CGAL_MEMORY_H
