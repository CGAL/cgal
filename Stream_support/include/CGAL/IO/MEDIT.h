// Copyright (c) 2015-2020  Geometry Factory
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Mael Rouxel-Labbé, Andreas Fabri

#ifndef CGAL_IO_MEDIT_H
#define CGAL_IO_MEDIT_H

#include <CGAL/assertions.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Container_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <array>
#include <map>
#include <boost/unordered_map.hpp>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

namespace IO {
namespace internal {

template<typename SurfacePatchIndex>
struct Facet
{
  int v0, v1, v2;
  SurfacePatchIndex surface_patch_index;
};

template <class FWPI>
struct get_surface_path_index_type
{
  using type = int;
};

template <typename SurfacePatchIndex>
struct get_surface_path_index_type<Facet<SurfacePatchIndex>>
{
  using type = SurfacePatchIndex;
};

template <typename I>
struct get_surface_path_index_type<std::array<I,4>>
{
  using type = I;
};

template <typename A, typename I>
struct get_surface_path_index_type<std::tuple<A,A,A,I>>
{
  using type = I;
};

template<typename CurveIndex>
struct Edge
{
  int v0, v1;
  CurveIndex curve_index;
};

template<typename CornerIndex>
struct Vertex_with_corner_index
{
  int v;
  CornerIndex corner_index;

  bool operator==(const Vertex_with_corner_index& rhs) const
  {
    return (v == rhs.v) && (corner_index == rhs.corner_index);
  }
};

template <typename T>
int get_zero(const T& t)
{
  return std::get<0>(t);
}

template <typename T>
int get_one(const T& t)
{
  return std::get<1>(t);
}

template <typename T>
int get_two(const T& t)
{
  return std::get<2>(t);
}

template<typename CornerIndex>
int get_zero(const Vertex_with_corner_index<CornerIndex>& t)
{
  return t.v;
}

template<typename CornerIndex>
CornerIndex get_one(const Vertex_with_corner_index<CornerIndex>& t)
{
  return t.corner_index;
}

template<typename CurveIndex>
int get_zero(const Edge<CurveIndex>& t)
{
  return t.v0;
}

template<typename CurveIndex>
int get_one(const Edge<CurveIndex>& t)
{
  return t.v1;
}

template<typename CurveIndex>
CurveIndex get_two(const Edge<CurveIndex>& t)
{
  return t.curve_index;
}

template<class PointWithIndexRange,
         class CellWithIndexRange,
         class FacetWithIndexRange, // either Facet or a tuple/array
         class EdgeWithIndexRange, // either Edge or a tuple/array
         class RidgeRange,
         class CornerRange> // either Vertex_with_corner_index or a tuple/pair/array
bool read_MEDIT(std::istream& is,
                PointWithIndexRange& points,
                CellWithIndexRange& cells,
                FacetWithIndexRange& facets,
                bool read_facets,
                EdgeWithIndexRange& edges,
                bool read_edges,
                RidgeRange& ridges,
                bool read_ridges,
                CornerRange& corners,
                bool read_corners,
                bool verbose,
                bool& is_CGAL_mesh)
{
  using Point_3 = typename PointWithIndexRange::value_type::first_type;
  using FT      = typename Kernel_traits<Point_3>::Kernel::FT;
  using Facet   = std::array<int, 3>;
  using Cell_with_ref = typename std::iterator_traits<typename CellWithIndexRange::const_iterator>::value_type;

  if(!is)
    return false;

  int dim;
  int nv, nf, ntet, ref, nvertices, nedges, nridges;
  int offset = static_cast<int>(points.size());
  std::string word;

  is >> word >> dim; // MeshVersionFormatted 1
  is >> word >> dim; // Dimension 3

  CGAL_assertion(dim == 3);

  std::string line;
  while(std::getline(is, line) && line != "End")
  {
    // remove trailing whitespace, in particular a possible '\r' from Windows
    // end-of-line encoding
    while(!line.empty() && std::isspace(line.back())) {
      line.pop_back();
    }
    if(line.empty())
      continue;

    // remove whitespaces at the beginning of the line
    for (std::size_t i=0; i<line.size(); ++i)
    {
      if (!std::isspace(line[i]))
      {
        if (i!=0)
          line = line.substr(i);
        break;
      }
    }

    if (line.at(0) == '#' &&
        line.find("CGAL::Mesh_complex_3_in_triangulation_3") != std::string::npos)
    {
      is_CGAL_mesh = true; // with CGAL meshes, domain 0 should be kept
      continue;
    }

    // skip non-CGAL comments
    if (line.at(0)=='#') continue;

    if(line.find("Vertices") != std::string::npos)
    {
      is >> nv;
      if(verbose)
        std::cerr << "Reading "<< nv << " vertices" << std::endl;
      for(int i=0; i<nv; ++i)
      {
        FT x,y,z;
        if(!(is >> x >> y >> z >> ref))
        {
          if(verbose)
            std::cerr << "Issue while reading vertices" << std::endl;
          return false;
        }
        points.emplace_back(std::make_pair(Point_3(x,y,z),ref));
      }
    }

    if(line.find("Triangles") != std::string::npos)
    {
      if(read_facets){
        is >> nf;
        CGAL::internal::reserve(facets, nf);

        if(verbose)
          std::cerr << "Reading "<< nf << " triangles" << std::endl;

        for(int i=0; i<nf; ++i)
        {
          int n[3];
          typename get_surface_path_index_type<typename std::iterator_traits<typename FacetWithIndexRange::iterator>::value_type>::type surface_patch_id;
          if(!(is >> n[0] >> n[1] >> n[2] >> surface_patch_id))
          {
            if(verbose)
              std::cerr << "Issue while reading triangles" << std::endl;
            return false;
          }
          Facet facet;
          facet[0] = offset + n[0] - 1;
          facet[1] = offset + n[1] - 1;
          facet[2] = offset + n[2] - 1;

          if(verbose)
            std::cout << "Looking at face #" << i << ": " << n[0] << " " << n[1] << " " << n[2] << std::endl;

          CGAL_warning_code(
          for(int j=0; j<3; ++j)
            for(int k=0; k<3; ++k)
              if(j != k)
                CGAL_warning(n[j] != n[k]);
          )

          // find the circular permutation that puts the smallest index in the first place.
          int n0 = (std::min)({facet[0],facet[1], facet[2]});
          while(facet[0] != n0)
          {
            std::rotate(std::begin(facet), std::next(std::begin(facet)), std::end(facet));
          }
          facets.push_back({facet[0], facet[1], facet[2], surface_patch_id});
        }
      }else{
        is >> nf;
        std::string buffer;
        for(int i=0; i<nf; ++i)
          std::getline(is, buffer);
      }
    }
    if(line.find("Tetrahedra") != std::string::npos)
    {
      is >> ntet;

      if(verbose)
        std::cerr << "Reading "<< ntet << " cells" << std::endl;

      for(int i=0; i<ntet; ++i)
      {
        int n[4];
        int reference;

        if(!(is >> n[0] >> n[1] >> n[2] >> n[3] >> reference))
        {
          if(verbose)
            std::cerr << "Issue while reading cells" << std::endl;
          return false;
        }

        if(verbose)
          std::cout << "Looking at tet #" << i << ": " << n[0] << " " << n[1] << " " << n[2] << " " << n[3] << std::endl;

        CGAL_warning_code(
        for(int j=0; j<4; ++j)
          for(int k=0; k<4; ++k)
            if(j != k)
              CGAL_warning(n[j] != n[k]);
        )

        Cell_with_ref t;
        t[0] = offset + n[0] - 1;
        t[1] = offset + n[1] - 1;
        t[2] = offset + n[2] - 1;
        t[3] = offset + n[3] - 1;
        t[4] = reference;

        cells.push_back(t);
      }
    }

    if(line.find("Corners") != std::string::npos)
    {
        is >> nvertices;
        int ci;
        for(int i=0; i<nvertices; ++i){
          is >> ci;
          if(read_corners){
           corners.push_back(offset + i);
          }
        }
    }

    if(line.find("Edges") != std::string::npos)
    {
      if(read_edges){
        is >> nedges;
        CGAL::internal::reserve(edges, nedges);

        if(verbose)
          std::cerr << "Reading "<< nedges << " edges" << std::endl;

      if(verbose && nedges == 0)
        std::cerr << "Warning: Edges section is empty" << std::endl;

      for(int i = 0; i < nedges; ++i)
      {
        int n[2];
        /* AF: was Curve_index*/ int  curve_index = 0;
        if(!(is >> n[0] >> n[1] >> curve_index))
        {
          if(verbose)
            std::cerr << "Issue while reading edges" << std::endl;
          return false;
        }
        edges.push_back({offset + n[0] - 1, offset + n[1] - 1, curve_index});
        CGAL_assertion(edges.size() == static_cast<std::size_t>(i + 1));
      }
    }
    else
    {
      is >> nedges;
      std::string buffer;
      for(int i=0; i<nedges; ++i)
        std::getline(is, buffer);
    }

    if(line.find("Ridges") != std::string::npos)
    {
      is >> nridges;
        int ci;
        for(int i=0; i<nridges; ++i){
          is >> ci;
          if(read_ridges){
           ridges.push_back(offset + i);
          }
        }
    }
  }

  }

  if (verbose)
  {
    std::cout << points.size() - std::size_t(offset) << " points" << std::endl;
    std::cout << cells.size() << " cells" << std::endl;
    std::cout << facets.size() << " facets" << std::endl;
    std::cout << edges.size() << " edges" << std::endl;
    std::cout << ridges.size() << " ridges" << std::endl;
    std::cout << corners.size() << " corners" << std::endl;
  }

  if(cells.empty())
    return false;

  return true;
}

} // namespace internal



/*!
 * \ingroup PkgStreamSupportIoFuncsMEDIT
 *
 * \brief reads the content of `is` into `points` and the named parameters using the \ref IOStreamMedit.
 *
 * The subdomain index corresponds to what \medit calls cell reference,
 * facet patch index corresponds to triangle reference,
 * edge curve index to edge reference,
 * and vertex corner index to vertex reference.
 *
 * \note Currently, only tetrahedral cells are supported.
 * \note The cell soup is not cleared, and the data from the stream are appended.
 *
 * \tparam PointRange a model of the concept `BackInsertionSequence` whose value type is the point type
 *
 * \param is the input stream
 * \param points points of the soup of cells
 *
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *
 *
 *   \cgalParamNBegin{tetrahedra}
 *     \cgalParamDescription{a non-const reference wrapper of a container of quintuples of integers that will be filled by this function.
 *                           Each element represents a tetrahedron with vertices corresponding to the first four integers, and the last integer being the subdomin index of the tetrahedron.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of `BackInsertionSequence` with a value type constructible using a braced initializer list of five integers.}
 *     \cgalParamDefault{facets are ignored}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{facets}
 *     \cgalParamDescription{a non-const reference wrapper of a container of quadruples of integers that will be filled by this function.
 *                           Each element represents a facet with vertices corresponding to the first three integers, and the last integer being the surface patch index of the facet.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of `BackInsertionSequence` with a value type constructible using a braced initializer list of four integers.}
 *     \cgalParamDefault{facets are ignored}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{edges}
 *     \cgalParamDescription{a non-const reference wrapper of a container of triples of integers that will be filled by this function.
 *                           Each element represents an edge with vertices corresponding  to the first two integers, and the last integer being the curve index of the edge.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of `BackInsertionSequence`  with a value type constructible using a braced initializer list of three integers.}
 *     \cgalParamDefault{edges are ignored}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{ridges}
 *     \cgalParamDescription{a non-const reference wrapper of a container of integers that will be filled by this function.
 *                           Each element in the container indicates the subdomain index of the corresponding cell at the same position.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of `BackInsertionSequence` able to store `int`.}
 *     \cgalParamDefault{ridges are ignored}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{corners}
 *     \cgalParamDescription{a non-const reference wrapper of a container of pairs of integers that will be filled by this function.
 *                           Each element represents a corner at the vertex corresponding to the first integer, and the last integer being the corner index.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of `BackInsertionSequence`  with a value type constructible using a braced initializer list of two integers.}
 *     \cgalParamDefault{vertices are ignored}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{verbose}
 *     \cgalParamDescription{if `true`, prints information about the reading process.}
 *     \cgalParamType{Boolean}
 *     \cgalParamDefault{`false`}
 *   \cgalParamNEnd
 *
 * \cgalNamedParamsEnd
 *
 * \returns `true` if the reading was successful, `false` otherwise.
 *
 *  \see \ref IOStreamMedit
 */

template<class PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_MEDIT(std::istream& is,
                PointRange& points,
                const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::get_parameter_reference;

  const bool verbose = choose_parameter(get_parameter(np, internal_np::verbose), false);

  // cells
  std::vector<std::array<int, 5>> default_tetrahedra;
  using Tetrahedra = typename internal_np::Lookup_named_param_def<internal_np::tetrahedra_t, CGAL_NP_CLASS, std::vector<std::array<int,5>>>::reference;
  Tetrahedra tetrahedra = choose_parameter(get_parameter_reference(np, internal_np::tetrahedra), default_tetrahedra);

  // triangles
  std::vector<std::array<int, 4>> default_triangles;
  using Triangles = typename internal_np::Lookup_named_param_def<internal_np::triangles_t, CGAL_NP_CLASS, std::vector<std::array<int,4>>>::reference;
  Triangles triangles = choose_parameter(get_parameter_reference(np, internal_np::triangles), default_triangles);

  // edges
  std::vector<std::array<int, 3>> default_edges;
  using Edges = typename internal_np::Lookup_named_param_def<internal_np::edges_t, CGAL_NP_CLASS, std::vector<std::array<int,3>>>::reference;
  Edges edges = choose_parameter(get_parameter_reference(np, internal_np::edges), default_edges);

  // ridges
  std::vector<int> default_ridges;
  using Ridges = typename internal_np::Lookup_named_param_def<internal_np::ridges_t, CGAL_NP_CLASS, std::vector<int>>::reference;
  Ridges ridges = choose_parameter(get_parameter_reference(np, internal_np::ridges), default_ridges);

  // corners
  std::vector<std::array<int, 2>> default_corners;
  using Corners = typename internal_np::Lookup_named_param_def<internal_np::corners_t, CGAL_NP_CLASS, std::vector<std::array<int,2>>>::reference;
  Corners corners = choose_parameter(get_parameter_reference(np, internal_np::corners), default_corners);


  constexpr bool is_triangles_map = !parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::triangles_t>::value;

  constexpr bool is_edges_map = !parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::edges_t>::value;

  constexpr bool is_ridges_map = !parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::ridges_t>::value;

  constexpr bool is_corners_map = !parameters::is_default_parameter<CGAL_NP_CLASS, internal_np::corners_t>::value;


  bool is_CGAL_mesh;

  return internal::read_MEDIT(is, points, tetrahedra,
                              triangles, is_triangles_map,
                              edges, is_edges_map,
                              ridges, is_ridges_map,
                              corners, is_corners_map,
                              verbose, is_CGAL_mesh);
}


/*!
 * \ingroup PkgStreamSupportIoFuncsMEDIT
 *
 * \brief writes a soup of indexed cells using the \ref IOStreamMedit.
 *
 * The subdomain index corresponds to what \medit calls cell reference,
 * facet patch index corresponds to triangle reference,
 * edge curve index to edge reference,
 * and vertex corner index to vertex reference.

 * \note Currently only tetrahedral cells are supported.
 *
 * \tparam PointRange a model of the concept `ConstRange` whose value type is the point type
 * \tparam CellRange a model of the concept `ConstRange`
 *                   whose `value_type` is a model of the concept `RandomAccessContainer`
 *                   whose `value_type` is `std::size_t`.
 *
 * \param os the output stream
 * \param points points of the soup of cells
 * \param cells a range where each element is a cell described by indices into `points`
 *
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{ridges}
 *     \cgalParamDescription{a const reference wrapper of a container of integers of the same size as `cells`.
 *                           Each element in the container indicates the subdomain index of the corresponding cell at the same position.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of the concept `SequenceContainer` of integer.}
 *     \cgalParamDefault{all cells will be given the subdomain index `1`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{facets}
 *     \cgalParamDescription{a const reference wrapper of a container of quadruples of integers that will be written by this function.
 *                           Each element corresponds to a facet with vertices corresponding to the three first integers, and the last integer being the surface patch index of the facet.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of `SequenceContainer` with a value type where the elements of the quadruples can be accessed with `std::get<int>()`}
 *     \cgalParamDefault{facets are not written}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{edges}
 *     \cgalParamDescription{a const reference wrapper of a container of triples of integers that will be written by this function.
 *                           Each element corresponds to an edge with vertices corresponding  to the two first integers, and the last integer being the curve index of the edge.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of `SequenceContainer` with a value type where the elements of the triples can be accessed with `std::get<int>()`}
 *     \cgalParamDefault{edges are not written}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{corners}
 *     \cgalParamDescription{a const reference wrapper of a container of pairs of integers that will be written by this function.
 *                           Each element corresponds to a corner at the vertex corresponding to the first integer, the second integer being the corner index.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of `SequenceContainer` with a value type where the elements of the pairs can be accessed with `std::get<int>()`}
 *     \cgalParamDefault{corners are not written}
 *   \cgalParamNEnd
 *
 * \cgalNamedParamsEnd
 *
 * \returns `true` if the writing was successful, `false` otherwise.
 *
 *  \see \ref IOStreamMedit
 */
template<class PointRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_MEDIT(std::ostream& os,
                 const PointRange& points,
                 const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                 , std::enable_if_t<!is_named_function_parameter<PointRange>>* = nullptr
#endif
                 )

{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::get_parameter_reference;
  using parameters::is_default_parameter;

  using Point_3 = typename PointRange::value_type;

  if(!os)
    return false;


  std::vector<int> default_ridges;
  std::vector<int> default_corners;
  std::vector<std::array<int,3>> default_edges;
  std::vector<std::array<int,4>> default_facets;
  std::vector<std::array<int,5>> default_cells;

  using Cells = typename internal_np::Lookup_named_param_def<internal_np::tetrahedra_t, CGAL_NP_CLASS, std::vector<std::array<int,5>>>::reference;
  Cells cells = choose_parameter(get_parameter_reference(np, internal_np::tetrahedra), default_cells);

  using Facets = typename internal_np::Lookup_named_param_def<internal_np::triangles_t, CGAL_NP_CLASS, std::vector<std::array<int,4>>>::reference;
  Facets facets = choose_parameter(get_parameter_reference(np, internal_np::triangles), default_facets);

  using Edges = typename internal_np::Lookup_named_param_def<internal_np::edges_t, CGAL_NP_CLASS, std::vector<std::array<int,3>>>::reference;
  Edges edges = choose_parameter(get_parameter_reference(np, internal_np::edges), default_edges);

  using Ridges = typename internal_np::Lookup_named_param_def<internal_np::ridges_t, CGAL_NP_CLASS, std::vector<int>>::reference;
  Ridges ridges = choose_parameter(get_parameter_reference(np, internal_np::ridges), default_ridges);

  using corners = typename internal_np::Lookup_named_param_def<internal_np::corners_t, CGAL_NP_CLASS, std::vector<internal::Vertex_with_corner_index<int>>>::reference;
  corners vertices = choose_parameter(get_parameter_reference(np, internal_np::corners), default_corners);


  os << "MeshVersionFormatted 1\nDimension 3\n";

  os << "Vertices\n";
  os << points.size() << "\n";
  for (const auto& p : points){

    os << p.first << " " << p.second << "\n";
  }

  os << "Tetrahedra\n";
  os << cells.size() << "\n";
  for (std::size_t k=0; k<cells.size(); ++k)
  {
    os << cells[k][0]+1 << " "
       << cells[k][1]+1 << " "
       << cells[k][2]+1 << " "
       << cells[k][3]+1 << " "
       << cells[k][4]   << "\n";
  }

  if(!edges.empty()){
    os << "Edges\n" << edges.size() << "\n";
    for(const auto& e : edges)
      os << internal::get_zero(e) +1 << " " << internal::get_one(e) + 1 << " " << internal::get_two(e)  << "\n";

    os << "Ridges\n" << edges.size() << "\n";
    for(std::size_t i = 1; i <= edges.size(); ++i)
      os << i << "\n";
  }

  if(!facets.empty()){
    os << "Triangles\n" << facets.size() << "\n";
    for(const auto& f : facets)
      os << std::get<0>(f) +1 << " " << std::get<1>(f) + 1 << " " << std::get<2>(f) + 1 << " " <<  std::get<3>(f)  << "\n";
  }

  if(!vertices.empty()){
    os << "Corners\n" << vertices.size() << "\n";
    for(const auto& c : vertices)
      os << c + 1 << "\n";
  }

  if(!ridges.empty()){
    os << "Ridges\n" << ridges.size() << "\n";
    for(const auto& r : ridges)
      os << r + 1 << "\n";
  }

  os <<"End\n";
  return true;
}

} // namespace IO

} // namespace CGAL

#endif // CGAL_IO_MEDIT_H
