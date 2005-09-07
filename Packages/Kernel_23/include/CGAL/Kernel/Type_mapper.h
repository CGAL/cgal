// Copyright (c) 2005  Utrecht University (The Netherlands),
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_KERNEL_TYPE_MAPPER_H
#define CGAL_KERNEL_TYPE_MAPPER_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

// This is a tool to obtain the K2::Point_2 from K1 and K1::Point_2.
// Similarly for other kernel types.

template < typename T, typename K1, typename K2 >
struct Type_mapper
{
  typedef T type; // By default, assume same type (e.g. Object).
};

template < typename K1, typename K2 >
struct Type_mapper < typename K1::Point_2, K1, K2 >
{
  typedef typename K2::Point_2 type;
};

template < typename K1, typename K2 >
struct Type_mapper < typename K1::Vector_2, K1, K2 >
{
  typedef typename K2::Vector_2 type;
};


template < typename K1, typename K2 >
struct Type_mapper < typename K1::Direction_2, K1, K2 >
{
  typedef typename K2::Direction_2 type;
};

template < typename K1, typename K2 >
struct Type_mapper < typename K1::Segment_2, K1, K2 >
{
  typedef typename K2::Segment_2 type;
};

template < typename K1, typename K2 >
struct Type_mapper < typename K1::Ray_2, K1, K2 >
{
  typedef typename K2::Ray_2 type;
};


template < typename K1, typename K2 >
struct Type_mapper < typename K1::Line_2, K1, K2 >
{
  typedef typename K2::Line_2 type;
};


template < typename K1, typename K2 >
struct Type_mapper < typename K1::Triangle_2, K1, K2 >
{
  typedef typename K2::Triangle_2 type;
};


template < typename K1, typename K2 >
struct Type_mapper < typename K1::Iso_rectangle_2, K1, K2 >
{
  typedef typename K2::Iso_rectangle_2 type;
};

template < typename K1, typename K2 >
struct Type_mapper < typename K1::Circle_2, K1, K2 >
{
  typedef typename K2::Circle_2 type;
};

// TODO : more specializations...

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_TYPE_MAPPER_H
