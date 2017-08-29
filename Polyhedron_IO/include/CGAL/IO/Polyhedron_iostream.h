// Copyright (c) 1997  ETH Zurich (Switzerland).
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
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_POLYHEDRON_IOSTREAM_H
#define CGAL_IO_POLYHEDRON_IOSTREAM_H 1

#include <CGAL/license/Polyhedron.h>


#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/print_OFF.h>
#include <CGAL/IO/scan_OFF.h>
#include <iostream>

namespace CGAL {


template < class Traits,
           class Items,
           template < class T, class I, class A>
           class HDS, class Alloc>
bool
write_off( std::ostream& out, const Polyhedron_3<Traits,Items,HDS,Alloc>& P){
  // writes P to `out' in PRETTY, ASCII or BINARY format
    // as the stream indicates.
    File_header_OFF header( is_binary( out), ! is_pretty( out), false);
    CGAL::print_polyhedron_with_header_OFF( out, P, header);
    return out.good();
}


template < class Traits,
           class Items,
           template < class T, class I, class A>
           class HDS, class Alloc>
bool 
read_off(std::istream& in,
         Polyhedron_3<Traits,Items,HDS,Alloc>& P) {
    // reads a polyhedron from `in' and appends it to P.
    CGAL::scan_OFF( in, P);
    return in.good();
}


template < class Traits,
           class Items,
           template < class T, class I, class A>
           class HDS, class Alloc>
std::ostream&
operator<<( std::ostream& out, const Polyhedron_3<Traits,Items,HDS,Alloc>& P)
{
  write_off(out,P);
  return out;
}


template < class Traits,
           class Items,
           template < class T, class I, class A>
           class HDS, class Alloc>
std::istream& operator>>(std::istream& in,
                         Polyhedron_3<Traits,Items,HDS,Alloc>& P)
{
  read_off(in,P);
  return in;
}

} //namespace CGAL
#endif // CGAL_IO_POLYHEDRON_IOSTREAM_H //
// EOF //
