// Copyright (c) 2007  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// $URL$
// $Id$
// 
// Author(s)     : Sylvain Pion

// ---------------------------------------------------------------------
// This program is used by install_cgal.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set, if a compiler cannot properly perform a partial
//| specialization relying on deducable arguments.
//| This bug appears for example on Sunpro 5.9 beta for x86 linux.
//| (more precisely: Sun C++ 5.9 Linux_i386 Build40_1 2007/02/09).

template < typename T, typename K1, typename K2 >
struct Type_mapper
{
  typedef T type;
};

template < typename K1, typename K2 >
struct Type_mapper < typename K1::Point_2, K1, K2 >
{ typedef typename K2::Point_2 type; };

template < typename K1, typename K2 >
struct Type_mapper < typename K1::Point_3, K1, K2 >
{ typedef typename K2::Point_3 type; };

template < typename K >
struct Point_2 {};

template < typename K >
struct Point_3 {};

template < typename NT >
struct Cartesian
{
  typedef Cartesian<NT>  Self;
  typedef ::Point_2<Self>  Point_2;
  typedef ::Point_3<Self>  Point_3;
};

int main()
{
  typedef Cartesian<int> CI;
  typedef Cartesian<double> CD;
  typedef CI::Point_2 P2I;
  typedef Type_mapper<P2I, CI, CD>::type P2D;
  return 0;
}
