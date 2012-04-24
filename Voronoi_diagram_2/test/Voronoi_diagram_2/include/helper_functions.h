// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef VDA_HELPER_FUNCTIONS_H
#define VDA_HELPER_FUNCTIONS_H 1

template<class T> void kill_warning(const T&) {}

template<class T>
void test_is_convertible_to(const T&) {}

template<class Iterator>
void test_iterator(Iterator first, Iterator beyond)
{
  for (Iterator it = first; it != beyond; ++it) {}
  for (Iterator it = first; it != beyond; it++) {}
  for (Iterator it = beyond; it != first; --it) {}
  for (Iterator it = beyond; it != first; it--) {}

  bool b = first == beyond;
  kill_warning(b);
}

template<class Circulator>
void test_circulator(Circulator start)
{
  Circulator c = start;
  do { c++; } while ( c != start );
  do { ++c; } while ( c != start );
  do { c--; } while ( c != start );
  do { --c; } while ( c != start );
  bool b = c == start;
  kill_warning(b);
}

#endif // VDA_HELPER_FUNCTIONS_H
