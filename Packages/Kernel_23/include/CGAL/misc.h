// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Andreas Fabri
//                 Stefan Schirra
 

#ifndef CGAL_MISC_H
#define CGAL_MISC_H

CGAL_BEGIN_NAMESPACE

// A helper class:

template <class Target, class Source>
struct converter
{
    static inline Target do_it(const Source& s)
    { return static_cast<Target>(s); }
};

template <class Target, class Source>
inline
Target
convert_to (const Source& s)
{ return converter<Target, Source>::do_it(s); }

CGAL_END_NAMESPACE

#endif // CGAL_MISC_H
