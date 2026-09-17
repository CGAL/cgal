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
#include <CGAL/Named_function_parameters.h>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <array>
#include <boost/unordered_map.hpp>

namespace CGAL {

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// Read

namespace IO {
namespace internal {

template<typename SurfacePatchIndex>
struct Facet_with_index
{
  int v0, v1, v2;
  SurfacePatchIndex surface_patch_index;
};

template<typename CurveIndex>
struct Edge_with_index
{
  int v0, v1;
  CurveIndex curve_index;
};

template<typename CornerIndex>
struct Corner_with_index
{
  int v;

  bool operator==(const Corner_with_index& rhs) const
  {
    return (v  == rhs.v);
  }
};


template<class PointRange,
         class TetrahedronRange,
         class FacetWithIndex, // either Facet_with_index or a tuple/array
         class EdgeWithIndexRange, // either Edge_with_index or a tuple/array
         class CornerWithIndexRange> // either Corner_with_index or a tuple/pair/array
bool read_MEDIT(std::istream& is,
                PointRange& points,
                TetrahedronRange& tetrahedra,
                std::vector<int>& subdomains,
                std::vector<FacetWithIndex>& facets_with_indices,
                bool read_facets_with_indices,
                EdgeWithIndexRange& edges_with_indices,
                CornerWithIndexRange& corners_with_indices,
                bool verbose,
                bool& is_CGAL_mesh)
{
  using Point_3 = typename PointRange::value_type;
  using FT = typename Kernel_traits<Point_3>::Kernel::FT;
  using Facet        = std::array<int, 3>;
  using Tet_with_ref = typename std::iterator_traits<typename TetrahedronRange::const_iterator>::value_type;

  if(!is)
    return false;

  int dim;
  int nv, nf, ntet, ref, ncorners, nedges;
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
        points.emplace_back(x,y,z);
      }
    }

    if(line.find("Triangles") != std::string::npos)
    {
      if(read_facets_with_indices){
        is >> nf;
        facets_with_indices.reserve(nf);

        if(verbose)
          std::cerr << "Reading "<< nf << " triangles" << std::endl;

        for(int i=0; i<nf; ++i)
        {
          int n[3];
          int surface_patch_id;
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
          facets_with_indices.push_back({facet[0], facet[1], facet[2], surface_patch_id});
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
        std::cerr << "Reading "<< ntet << " tetrahedra" << std::endl;

      for(int i=0; i<ntet; ++i)
      {
        int n[4];
        int reference;

        if(!(is >> n[0] >> n[1] >> n[2] >> n[3] >> reference))
        {
          if(verbose)
            std::cerr << "Issue while reading tetrahedra" << std::endl;
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

        Tet_with_ref t;
        t[0] = offset + n[0] - 1;
        t[1] = offset + n[1] - 1;
        t[2] = offset + n[2] - 1;
        t[3] = offset + n[3] - 1;

        tetrahedra.push_back(t);
        subdomains.push_back(reference);
      }
    }

    if(line.find("Corners") != std::string::npos)
    {
      is >> ncorners;
      if(verbose && ncorners == 0)
        std::cerr << "Warning: Corners section is ignored" << std::endl;

      for(int i = 0; i < ncorners; ++i)
      {
        int n;
        if(!(is >> n))
        {
          if(verbose)
            std::cerr << "Issue while reading corners" << std::endl;
          return false;
        }
        // typename CornerWithIndex::value_type cwi = {offset + n};
        corners_with_indices.push_back( {offset + n - 1 } );
      }
    }

    if(line.find("Edges") != std::string::npos)
    {
      is >> nedges;
      if(verbose && nedges == 0)
        std::cerr << "Warning: Edges section is ignored" << std::endl;

      for(int i = 0; i < nedges; ++i)
      {
        int n[2], curve_index;
        if(!(is >> n[0] >> n[1] >> curve_index))
        {
          if(verbose)
            std::cerr << "Issue while reading edges" << std::endl;
          return false;
        }
        edges_with_indices.push_back({offset + n[0] - 1, offset + n[1] - 1, curve_index});
        CGAL_assertion(edges_with_indices.size() == static_cast<std::size_t>(i + 1));
      }
    }

  }

  if (verbose)
  {
    std::cout << points.size() - std::size_t(offset) << " points" << std::endl;
    std::cout << tetrahedra.size() << " cells" << std::endl;
    std::cout << facets_with_indices.size() << " border facets" << std::endl;
    std::cout << edges_with_indices.size() << " edges" << std::endl;
    std::cout << corners_with_indices.size() << " corners" << std::endl;
  }

  if(tetrahedra.empty())
    return false;

  CGAL_assertion(tetrahedra.size() == subdomains.size());

  return true;
}

} // namespace internal



/*!
 * \ingroup PkgStreamSupportIoFuncsMEDIT
 *
 * \brief reads the content of `is` into `points` and `tetrahedra`.
 *
 * See \cgalCite{frey:inria-00069921} for a comprehensive description of the medit (`.mesh`) file format.
 *
 * \attention The tetrahedron soup is not cleared, and the data from the stream are appended.
 *
 * \tparam PointRange a model of the concept `BackInsertionSequence` whose value type is the point type
 * \tparam TetrahedronRange a model of the concept `BackInsertionSequence` whose `value_type` is `std::array<int,4>`
 *
 * \param is the input stream
 * \param points points of the soup of cells
 * \param tetrahedra each element in it describes a tetrahedron
 *        using the indices of the points in `points`
 *
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *
 *   \cgalParamNBegin{subdomains}
 *     \cgalParamDescription{a non-const reference wrapper to a container of integers that will be filled by this function.
 *                           Each element in the container indicates the subdomain index of the corresponding tetrahedron at the same position.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of `BackInsertionSequence` able to store `int`.}
 *     \cgalParamDefault{subdomains are ignored}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{facets_with_indices}
 *     \cgalParamDescription{a non-const reference wrapper to a container of quadruple of integers that will be filled by this function.
 *                           Each element corresponds to a facet with vertices corresponding to the three first integers, and the last integer being the surface patch index of the facet.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of `BackInsertionSequence` able to store  quadruple of `int`.}
 *     \cgalParamDefault{facets are ignored}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{edges_with_indices}
 *     \cgalParamDescription{a non-const reference wrapper to a container of triple of integers that will be filled by this function.
 *                           Each element corresponds to an edge with vertices corresponding  to the two first integers, and the last integer being the curve index of the edge.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of `BackInsertionSequence` able to store  triple of `int`.}
 *     \cgalParamDefault{edges are ignored}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{corners_with_indices}
 *     \cgalParamDescription{a non-const reference wrapper to a container of pair of integers that will be filled by this function.
 *                           Each element corresponds to a corner at the vertex corresponding to the first integer, the last one being the corner index.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of `BackInsertionSequence` able to store  pair of `int`.}
 *     \cgalParamDefault{corners are ignored}
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

template<class PointRange, class TetrahedronRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool read_MEDIT(std::istream& is,
                PointRange& points,
                TetrahedronRange& tetrahedra,
                const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::get_parameter_reference;

  const bool verbose = choose_parameter(get_parameter(np, internal_np::verbose), false);

  // subdomains
  std::vector<int> default_subdomains;
  using Subdomains = typename internal_np::Lookup_named_param_def<internal_np::subdomains_t, CGAL_NP_CLASS, std::vector<int>>::reference;
  Subdomains subdomains = choose_parameter(get_parameter_reference(np, internal_np::subdomains), default_subdomains);

  // facets
  std::vector<std::array<int, 4>> default_facets_with_indices;
  using Facets_with_indices = typename internal_np::Lookup_named_param_def<internal_np::facets_with_indices_t, CGAL_NP_CLASS, std::vector<std::array<int,4>>>::reference;
  Facets_with_indices facets_with_indices = choose_parameter(get_parameter_reference(np, internal_np::facets_with_indices), default_facets_with_indices);

  // edges
  std::vector<std::array<int, 3>> default_edges_with_indices;
  using Edges_with_indices = typename internal_np::Lookup_named_param_def<internal_np::edges_with_indices_t, CGAL_NP_CLASS, std::vector<std::array<int,3>>>::reference;
  Edges_with_indices edges_with_indices = choose_parameter(get_parameter_reference(np, internal_np::edges_with_indices), default_edges_with_indices);

  // corners
  std::vector<std::array<int, 2>> default_corners_with_indices;
  using Corners_with_indices = typename internal_np::Lookup_named_param_def<internal_np::corners_with_indices_t, CGAL_NP_CLASS, std::vector<std::array<int,2>>>::reference;
  Corners_with_indices corners_with_indices = choose_parameter(get_parameter_reference(np, internal_np::corners_with_indices), default_corners_with_indices);

  constexpr bool read_facets_with_indices = false;

  bool is_CGAL_mesh;

  return internal::read_MEDIT(is, points, tetrahedra, subdomains,
                              facets_with_indices, read_facets_with_indices,
                              edges_with_indices, corners_with_indices,
                              verbose, is_CGAL_mesh);
}


/*!
 * \ingroup PkgStreamSupportIoFuncsMEDIT
 *
 * \brief writes the `points` and `tetrahedra`.
 *
 * See \cgalCite{frey:inria-00069921} for a comprehensive description of the medit (`.mesh`) file format.
 *
 *
 * \tparam PointRange a model of the concept `ConstRange` whose value type is the point type
 * \tparam TetrahedronRange a model of the concept `ConstRange`
 *                   whose `value_type` is a model of the concept `RandomAccessContainer`
 *                   whose `value_type` is `std::size_t`.
 *
 * \param os the output stream
 * \param points points of the soup of cells
 * \param tetrahedra each element in it describes a cell
 *        using the indices of the points in `points`
 *
 * \param np optional \ref bgl_namedparameters "Named Parameters" described below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{subdomains}
 *     \cgalParamDescription{a reference wrapper to a container of integers of the same size as `tetrahedra`.
 *                           Each element in the container indicates the subdomain index of the corresponding tetrahedron at the same position.}
 *     \cgalParamType{a `std::reference_wrapper` to a model of the concept `RandomAccessContainer` of integer.}
 *     \cgalParamDefault{all tetrahedra will have the subdomain id `1`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \returns `true` if the writing was successful, `false` otherwise.
 *
 *  \see \ref IOStreamMedit
 */
template<class PointRange, class TetrahedronRange, typename CGAL_NP_TEMPLATE_PARAMETERS>
bool write_MEDIT(std::ostream& os,
                 const PointRange& points,
                 const TetrahedronRange& tetrahedra,
                 const CGAL_NP_CLASS& np = parameters::default_values()
#ifndef DOXYGEN_RUNNING
                 , std::enable_if_t<!is_named_function_parameter<TetrahedronRange>>* = nullptr
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

  // subdomains
  std::vector<int> default_subdomains;
  std::vector<internal::Corner_with_index<int>> default_corners;
  std::vector<std::array<int,3>> default_edges;

  using Subdomains = typename internal_np::Lookup_named_param_def<internal_np::subdomains_t, CGAL_NP_CLASS, std::vector<int>>::reference;
  Subdomains subdomains = choose_parameter(get_parameter_reference(np, internal_np::subdomains), default_subdomains);

  using Corners_with_indices = typename internal_np::Lookup_named_param_def<internal_np::corners_with_indices_t, CGAL_NP_CLASS, std::vector<internal::Corner_with_index<int>>>::reference;
  Corners_with_indices corners = choose_parameter(get_parameter_reference(np, internal_np::corners_with_indices), default_corners);

  using Edges_with_indices = typename internal_np::Lookup_named_param_def<internal_np::edges_with_indices_t, CGAL_NP_CLASS, std::vector<std::array<int,3>>>::reference;
  Edges_with_indices edges = choose_parameter(get_parameter_reference(np, internal_np::edges_with_indices), default_edges);


  if constexpr (is_default_parameter<CGAL_NP_CLASS, internal_np::subdomains_t>::value)
    default_subdomains.resize(tetrahedra.size(), 1);

  os << "MeshVersionFormatted 1\nDimension 3\nVertices\n";
  os << points.size() << "\n";
  for (const Point_3& p : points)
    os << p << " 0\n";
  os << "Triangles\n0\nTetrahedra\n";
  os << tetrahedra.size() << "\n";
  for (std::size_t k=0; k<tetrahedra.size(); ++k)
    os << tetrahedra[k][0]+1 << " "
       << tetrahedra[k][1]+1 << " "
       << tetrahedra[k][2]+1 << " "
       << tetrahedra[k][3]+1 << " " << subdomains[k] << "\n";

  if(!corners.empty()){
    os << "Corners\n" << corners.size() << "\n";
    for(const auto& c : corners)
      os << c.v + 1 << "\n";
  }

  if(!edges.empty()){
    os << "Edges\n" << edges.size() << "\n";
    for(const auto& e : edges)
      os << std::get<0>(e) +1 << " " << std::get<1>(e) + 1 << " " << std::get<2>(e)  << "\n";
  }
  os <<"End\n";
  return true;
}

} // namespace IO

} // namespace CGAL

#endif // CGAL_IO_MEDIT_H
