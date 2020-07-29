// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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
