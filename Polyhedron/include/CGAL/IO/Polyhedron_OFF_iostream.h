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

#ifndef CGAL_IO_POLYHEDRON_OFF_IOSTREAM_H
#define CGAL_IO_POLYHEDRON_OFF_IOSTREAM_H

#include <CGAL/license/Polyhedron.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/IO/print_OFF.h>
#include <CGAL/IO/scan_OFF.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/Has_conversion.h>
#include <CGAL/property_map.h>

#include <fstream>
#include <iostream>

namespace CGAL {

namespace IO {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

template <class Traits,
          class Items,
          template < class T, class I, class A> class HDS,
          class Alloc, class CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(std::istream& in,
              Polyhedron_3<Traits, Items, HDS, Alloc>& P,
              const CGAL_BGL_NP_CLASS& np)
{
  typedef typename boost::graph_traits<Polyhedron_3<Traits, Items, HDS, Alloc> >::vertex_descriptor Vertex;

  // reads a polyhedron from `in' and appends it to P.
  typedef typename property_map_selector<Polyhedron_3<Traits, Items, HDS, Alloc>,
                                                      boost::vertex_point_t>::type      Def_VPM;
  typedef typename boost::property_traits<Def_VPM>::value_type                          Def_point;
  typedef typename Kernel_traits<Def_point>::Kernel                                     Def_kernel;

  typedef typename CGAL::GetVertexPointMap<Polyhedron_3<Traits, Items, HDS, Alloc>,
                                                        CGAL_BGL_NP_CLASS>::type        VPM;
  typedef typename boost::property_traits<VPM>::value_type                              Point;
  typedef typename Kernel_traits<Point>::Kernel                                         Kernel;

  typedef typename CGAL::internal::Converter_selector<Def_kernel, Kernel>::type         Converter;

  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  const bool verbose = choose_parameter(get_parameter(np, internal_np::verbose), true);

  if(!(is_default_parameter(get_parameter(np, internal_np::vertex_color_map))) ||
     !(is_default_parameter(get_parameter(np, internal_np::face_color_map))) ||
     !(is_default_parameter(get_parameter(np, internal_np::vertex_normal_map))) ||
     !(is_default_parameter(get_parameter(np, internal_np::vertex_texture_map))))
  {
    return CGAL::IO::internal::read_OFF_BGL(in, P, np);
  }

  CGAL::scan_OFF(in, P, verbose);

  if(!parameters::is_default_parameter(get_parameter(np, internal_np::vertex_point)))
  {
    Def_VPM def_vpm = get_property_map(CGAL::vertex_point, P);
    VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_property_map(CGAL::vertex_point, P));

    Converter to_vpm_type;
    for(Vertex v : vertices(P))
      put(vpm, v, to_vpm_type(get(def_vpm, v)));
  }

  return !in.fail();
}

template <class Traits,
          class Items,
          template < class T, class I, class A> class HDS,
          class Alloc>
bool read_OFF(std::istream& in, Polyhedron_3<Traits, Items, HDS, Alloc>& P)
{
  return read_OFF(in, P, parameters::all_default());
}

template <class Traits,
          class Items,
          template < class T, class I, class A> class HDS,
          class Alloc, class CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool read_OFF(const std::string& fname,
              Polyhedron_3<Traits, Items, HDS, Alloc>& P,
              const CGAL_BGL_NP_CLASS& np)
{
  std::ifstream in(fname);
  return read_OFF(in, P, np);
}

template <class Traits,
          class Items,
          template < class T, class I, class A> class HDS,
          class Alloc>
bool read_OFF(const std::string& fname, Polyhedron_3<Traits, Items, HDS, Alloc>& P)
{
  std::ifstream in(fname);
  return read_OFF(in, P, parameters::all_default());
}

} // namespace IO

template <class Traits,
          class Items,
          template < class T, class I, class A> class HDS,
          class Alloc>
std::istream& operator>>(std::istream& in, Polyhedron_3<Traits, Items, HDS, Alloc>& P)
{
  IO::read_OFF(in, P);
  return in;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

namespace IO {

template < class Traits,
           class Items,
           template < class T, class I, class A> class HDS,
           class Alloc, class CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(std::ostream& out,
               const Polyhedron_3<Traits, Items, HDS, Alloc>& P,
               const CGAL_BGL_NP_CLASS& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  if(!(is_default_parameter(get_parameter(np, internal_np::vertex_color_map))) ||
     !(is_default_parameter(get_parameter(np, internal_np::face_color_map))) ||
     !(is_default_parameter(get_parameter(np, internal_np::vertex_normal_map))) ||
     !(is_default_parameter(get_parameter(np, internal_np::vertex_texture_map))))
  {
    return CGAL::IO::internal::write_OFF_BGL(out, P, np);
  }

  // writes P to `out' in PRETTY, ASCII or BINARY format as the stream indicates.
  File_header_OFF header(is_binary(out), ! is_pretty(out), false);
  typename CGAL::GetVertexPointMap<Polyhedron_3<Traits, Items, HDS, Alloc>, CGAL_BGL_NP_CLASS>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, P));

  return CGAL::print_polyhedron_with_header_OFF(out, P, header, vpm);
}

template <class Traits,
          class Items,
          template < class T, class I, class A> class HDS,
          class Alloc>
bool write_OFF(std::ostream& out, const Polyhedron_3<Traits, Items, HDS, Alloc>& P)
{
  return write_OFF(out, P, parameters::all_default());
}

} // namespace IO

template <class Traits,
          class Items,
          template < class T, class I, class A>
          class HDS, class Alloc>
std::ostream& operator<<(std::ostream& out, const Polyhedron_3<Traits, Items, HDS, Alloc>& P)
{
  IO::write_OFF(out, P);
  return out;
}

} // namespace CGAL

#endif // CGAL_IO_POLYHEDRON_OFF_IOSTREAM_H
