
// Copyright (c) 2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Geert-Jan Giezeman


#ifndef CGAL_TRIANGLE_2_TRIANGLE_2_INTERSECTION_H
#define CGAL_TRIANGLE_2_TRIANGLE_2_INTERSECTION_H

#include <CGAL/Object.h>
#include <CGAL/Triangle_2_Triangle_2_do_intersect.h>

CGAL_BEGIN_NAMESPACE

template <class R>
Object
intersection(const Triangle_2<R> &tr1, const Triangle_2<R>&tr2);

// template <class R>
// bool
// do_intersect(const Triangle_2<R> &tr1, const Triangle_2<R>&tr2);

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Triangle_2_Triangle_2_intersection.C>
#endif
#endif
