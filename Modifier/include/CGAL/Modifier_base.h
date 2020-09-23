// Copyright (c) 1997
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
// Author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>

#ifndef CGAL_MODIFIER_BASE_H
#define CGAL_MODIFIER_BASE_H 1

namespace CGAL {

template <class R>
class Modifier_base {
    // Abstract Base Class for protected internal access.
    // It defines the common interface for all modifiers.
public:
    typedef R Representation;
    virtual void operator()( R& rep) = 0;
        // Postcondition: `rep' is a valid representation.
    virtual ~Modifier_base() {}
};

} //namespace CGAL

#endif // CGAL_MODIFIER_BASE_H //
// EOF //
