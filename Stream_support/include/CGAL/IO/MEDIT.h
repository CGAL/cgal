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

template<class PointRange, class CellRange, class SurfacePatchIndex_>
bool read_MEDIT(std::istream& is,
                PointRange& points,
                CellRange& finite_cells,
                std::vector<int>& subdomains,
                boost::unordered_map<std::array<int,3>,SurfacePatchIndex_ >& border_facets,
                bool read_border_facets,
                bool verbose,
                bool& is_CGAL_mesh)
{
  using Point_3 = typename PointRange::value_type;
  using FT = typename Kernel_traits<Point_3>::Kernel::FT;
  using Surface_patch_index = SurfacePatchIndex_;
  using Facet        = std::array<int, 3>;
  using Tet_with_ref = std::array<int, 4>;

  if(!is)
    return false;

  int dim;
  int nv, nf, ntet, ref;
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
      if(read_border_facets){
        bool has_negative_surface_patch_ids = false;
        Surface_patch_index max_surface_patch_id{0};
        is >> nf;

        if(verbose)
          std::cerr << "Reading "<< nf << " triangles" << std::endl;

        for(int i=0; i<nf; ++i)
        {
          int n[3];
          Surface_patch_index surface_patch_id;
          if(!(is >> n[0] >> n[1] >> n[2] >> surface_patch_id))
          {
            if(verbose)
              std::cerr << "Issue while reading triangles" << std::endl;
            return false;
          }
          has_negative_surface_patch_ids |= (surface_patch_id < 0);
          max_surface_patch_id = (std::max)(max_surface_patch_id, surface_patch_id);
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
          do
          {
            std::rotate(std::begin(facet), std::next(std::begin(facet)), std::end(facet));
          }
          while(facet[0] != n0);

          border_facets.emplace(facet, surface_patch_id);
        }
        if(has_negative_surface_patch_ids)
        {
          if(verbose)
            std::cerr << "Warning: negative surface patch ids" << std::endl;
          for(auto& facet_and_patch_id  : border_facets) {
            if(facet_and_patch_id.second < 0)
              facet_and_patch_id.second = max_surface_patch_id - facet_and_patch_id.second;
          }
        }
      }else{
        is >> nf;
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
        }
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

        finite_cells.push_back(t);
        subdomains.push_back(reference);
      }
    }
  }

  if (verbose)
  {
    std::cout << points.size() - std::size_t(offset) << " points" << std::endl;
    std::cout << border_facets.size() << " border facets" << std::endl;
    std::cout << finite_cells.size() << " cells" << std::endl;
  }

  if(finite_cells.empty())
    return false;

  CGAL_assertion(finite_cells.size() == subdomains.size());

  return true;
}

} // namespace internal



/*!
 * \ingroup PkgStreamSupportIoFuncsMEDIT
 *
 * \brief reads the content of `is` into `points` and `finite_cells`.
 *
 * See \cgalCite{frey:inria-00069921} for a comprehensive description of the medit (`.mesh`) file format.
 *
 * \attention The cell soup is not cleared, and the data from the stream are appended.
 *
 * \tparam PointRange a model of the concept `BackInsertionSequence` whose value type is the point type
 * \tparam CellRange a model of the concept `BackInsertionSequence` whose `value_type` is `std::array<int,4>`
 *
 * \param is the input stream
 * \param points points of the soup of cells
 * \param finite_cells each element in it describes a cell
 *        using the indices of the points in `points`
 * \param subdomains each element in it describes the subdomain index of the corresponding cell in `finite_cells`
 * \param verbose if `true`, prints information about the reading process
 *
 * \returns `true` if the reading was successful, `false` otherwise.
 *
 *  \see \ref IOStreamMedit
 */

template<class PointRange, class CellRange>
bool read_MEDIT(std::istream& is,
                PointRange& points,
                CellRange& finite_cells,
                std::vector<int>& subdomains,
                bool verbose = false)

{
  boost::unordered_map<std::array<int,3>,int > border_facets;
  constexpr bool read_border_facets = false;
  bool is_CGAL_mesh;
  return internal::read_MEDIT(is, points, finite_cells, subdomains, border_facets,
                              read_border_facets, verbose, is_CGAL_mesh);
}

} // namespace IO

} // namespace CGAL

#endif // CGAL_IO_MEDIT_H
