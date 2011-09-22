// Copyright (c) 1999,2000  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
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
// $URL$
// $Id$
// 
//
// Author(s)     : various

//| If a compiler (or assembler or linker) has problems with long names
//| CGAL_CFG_LONGNAME_BUG is set.

#ifdef _MSC_VER
#  pragma warning( error : 4503)
#endif

#define LONG_NAME \
Wwwwwwwwwooooooooo_vvvvvvveeeeeerrrryyyy_llllooooonnnnnnggggg\
_nnnnnaaaammmmmeeeeWwwwwwwwwooooooooo_vvvvvvveeeeeerrrryyyy_l\
lllooooonnnnnnggggg_nnnnnaaaammmmmeeeeWwwwwwwwwooooooooo_vvvv\
vvveeeeeerrrryyyy_llllooooonnnnnnggggg_nnnnnaaaammmmmeeeeWwww\
wwwwwooooooooo_vvvvvvveeeeeerrrryyyy_llllooooonnnnnnggggg_nnn\
nnaaaammmmmeeeeWwwwwwwwwooooooooo_vvvvvvveeeeeerrrryyyy_llllo\
oooonnnnnnggggg_nnnnnaaaammmmmeeeeWwwwwwwwwooooooooo_vvvvvvve\
eeeeerrrryyyy_llllooooonnnnnnggggg_nnnnnaaaammmmmeeee

template < class A >
struct LONG_NAME
{
  LONG_NAME (int i) : a(i) {}
  A a;
};

int main ()
{
  LONG_NAME< LONG_NAME< LONG_NAME< LONG_NAME< LONG_NAME< LONG_NAME<
  LONG_NAME< LONG_NAME< LONG_NAME< LONG_NAME< int > > > > > > > > > > a (1);
  (void) a;
  return 0;
}
