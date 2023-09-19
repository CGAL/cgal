// Copyright (c) 2005
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
// Author(s)     : Sylvain Pion

//| This flag is set, if a compiler cannot distinguish the signature
//| of overloaded function templates, which have one template parameter
//| to be passed explicitly when being called.
//|
//| This bug appears for example on g++ 3.3 and 3.4 (but not on more recent
//| g++ version). This bug appears also on Sun CC 5.90.

template < typename T >
struct A {};

template < typename T, typename U >
T enum_cast(const U&) { return T(); }

template < typename T, typename U >
T enum_cast(const A<U>&) { return T(); }

int main()
{
  A<double> a;
  int i = enum_cast<int>(a);
  (void) i;
  return 0;
}
