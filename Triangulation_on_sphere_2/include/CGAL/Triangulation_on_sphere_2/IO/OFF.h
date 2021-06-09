// Copyright (c) 1997-2021 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_TOS2_IO_OFF_H
#define CGAL_TOS2_IO_OFF_H

#include <CGAL/license/Triangulation_on_sphere_2.h>

#include <CGAL/triangulation_assertions.h>

#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/IO/helpers.h>

#include <fstream>
#include <string>
#include <unordered_map>

#ifdef DOXYGEN_RUNNING
#define CGAL_BGL_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_BGL_NP_CLASS NamedParameters
#define CGAL_DEPRECATED
#endif

namespace CGAL {

template <typename Gt, typename Tds>
class Triangulation_on_sphere_2;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Write

namespace IO {

/*!
  \ingroup PkgPointSet3IOOFF

  \brief writes the content of a triangulation on the sphere into an output stream in the \ref IOStreamOFF.

  \tparam Gt the geometric traits type of the triangulation
  \tparam Tds the triangulation data structure type of the triangulation
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param os the output stream
  \param dt the triangulation
  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e., how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{the precision of the stream `os`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \return `true` if the writing was successful, `false` otherwise.
 */
template <typename Gt, typename Tds, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(std::ostream& os,
               const CGAL::Triangulation_on_sphere_2<Gt, Tds>& dt,
               const CGAL_BGL_NP_CLASS& np)
{
  typedef Triangulation_on_sphere_2<Gt,Tds>             Tr;
  typedef typename Tr::Vertex_handle                    Vertex_handle;
  typedef typename Tr::Vertices_iterator                Vertex_iterator;
  typedef typename Tr::All_faces_iterator               Face_iterator;
  typedef typename Tr::size_type                        size_type;

  set_stream_precision_from_NP(os, np);

  const bool use_colors = parameters::choose_parameter(parameters::get_parameter(np, internal_np::output_color), true);

  typename Gt::Construct_point_3 cp3 = dt.geom_traits().construct_point_3_object();

  const size_type n = dt.number_of_vertices();

  std::stringstream output;
  set_stream_precision_from_NP(output, np);

  // write the vertices
  std::unordered_map<Vertex_handle, size_type> index_of_vertex;
  size_type i = 0;
  for(Vertex_iterator it = dt.vertices_begin(); it != dt.vertices_end(); ++it, ++i)
  {
    Vertex_handle vh = it;
    output << cp3(vh->point()) << " 0 0 0 \n"; // '0 0 0' for colors
    index_of_vertex[vh] = i;
  }

  CGAL_triangulation_assertion(i == n);

  size_type number_of_triangles = 0;
  for(Face_iterator fit = dt.all_faces_begin() ; fit != dt.all_faces_end() ; ++fit)
  {
    output << "3 "
           << index_of_vertex[fit->vertex(0)] << " "
           << index_of_vertex[fit->vertex(1)] << " "
           << index_of_vertex[fit->vertex(2)];
    if(use_colors)
    {
      if(fit->is_ghost())
        output << " 229 117 0\n";
      else
        output << " 0 117 229\n";
    }
    else
    {
      output << "\n";
    }

    ++number_of_triangles;
  }

  if(use_colors)
    os << "COFF \n";
  else
    os << "OFF \n";

  os << n << " "
     << number_of_triangles << " 0\n"
     << output.str();

  return !os.fail();
}

/// \cond SKIP_IN_MANUAL

template <typename Gt, typename Tds>
bool write_OFF(std::ostream& os, const CGAL::Triangulation_on_sphere_2<Gt, Tds>& dt)
{
  return write_OFF(os, dt, parameters::all_default());
}

/// \endcond

/*!
  \ingroup PkgPointSet3IOOFF

  \brief writes the content of a triangulation on the sphere into an output file in the \ref IOStreamOFF.

  \tparam Gt the geometric traits type of the triangulation
  \tparam Tds the triangulation data structure type of the triangulation
  \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  \param fname the path to the output file
  \param point_set the point set
  \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  \cgalNamedParamsBegin
    \cgalParamNBegin{stream_precision}
      \cgalParamDescription{a parameter used to set the precision (i.e., how many digits are generated) of the output stream}
      \cgalParamType{int}
      \cgalParamDefault{`6`}
    \cgalParamNEnd
  \cgalNamedParamsEnd

  \return `true` if the writing was successful, `false` otherwise.
*/
template <typename Gt, typename Tds, typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool write_OFF(const std::string& fname,
               const CGAL::Triangulation_on_sphere_2<Gt, Tds>& dt,
               const CGAL_BGL_NP_CLASS& np)
{
  std::ofstream os(fname); // stream precision will be set in the ostream overload
  return write_OFF(os, dt, np);
}

/// \cond SKIP_IN_MANUAL

template <typename Gt, typename Tds>
bool write_OFF(const std::string& fname, const CGAL::Triangulation_on_sphere_2<Gt, Tds>& dt)
{
  std::ofstream os(fname);
  return write_OFF(os, dt, parameters::all_default());
}

/// \endcond

} } // namespace CGAL::IO

#endif // CGAL_TOS2_IO_OFF_H
