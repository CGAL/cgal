// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

#ifndef CGAL_IO_POLYHEDRON_IOSTREAM_H
#define CGAL_IO_POLYHEDRON_IOSTREAM_H 1

#include <CGAL/license/Polyhedron.h>


#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/print_OFF.h>
#include <CGAL/IO/scan_OFF.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <iostream>

namespace CGAL {


template < class Traits,
           class Items,
           template < class T, class I, class A>
           class HDS, class Alloc, class NamedParameters>
bool
write_off( std::ostream& out, const Polyhedron_3<Traits,Items,HDS,Alloc>& P, const NamedParameters& np){
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // writes P to `out' in PRETTY, ASCII or BINARY format
    // as the stream indicates.
    File_header_OFF header( is_binary( out), ! is_pretty( out), false);
    typename CGAL::GetVertexPointMap<Polyhedron_3<Traits,Items,HDS,Alloc>, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, P));
    CGAL::print_polyhedron_with_header_OFF( out, P, header, vpm);
    return out.good();
}

template < class Traits,
           class Items,
           template < class T, class I, class A>
           class HDS, class Alloc>
bool
write_off( std::ostream& out, const Polyhedron_3<Traits,Items,HDS,Alloc>& P){
  return write_off(out, P, parameters::all_default());
}

template < class Traits,
           class Items,
           template < class T, class I, class A>
           class HDS, class Alloc,
           class NamedParameters>
bool
read_off(std::istream& in,
         Polyhedron_3<Traits,Items,HDS,Alloc>& P,
         NamedParameters np) {
    // reads a polyhedron from `in' and appends it to P.
    typedef typename CGAL::GetVertexPointMap<Polyhedron_3<Traits,Items,HDS,Alloc>, NamedParameters>::type Vpm;
    using parameters::choose_parameter;
    using parameters::get_parameter;

    Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_property_map(CGAL::vertex_point, P));
    CGAL::scan_OFF( in, P);
    if(!parameters::is_default_parameter(get_parameter(np, internal_np::vertex_point)))
    {
      typedef typename boost::graph_traits<Polyhedron_3<Traits,Items,HDS,Alloc> >::vertex_descriptor Vertex;
      typename property_map_selector<Polyhedron_3<Traits,Items,HDS,Alloc>, boost::vertex_point_t>::type
          def_vpm = get_property_map(CGAL::vertex_point, P);
      for(Vertex v : vertices(P))
      {
        put(vpm, v, get(def_vpm, v));
      }
    }
    return in.good();
}

template < class Traits,
           class Items,
           template < class T, class I, class A>
           class HDS, class Alloc>
bool
read_off(std::istream& in,
         Polyhedron_3<Traits,Items,HDS,Alloc>& P)
{
  return read_off(in, P, parameters::all_default());
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
