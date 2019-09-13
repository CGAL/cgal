// Copyright (c) 2019 GeometryFactory (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_INTERNAL_ADD_IMAGINARY_LAYER_H
#define CGAL_INTERNAL_ADD_IMAGINARY_LAYER_H

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Bbox_3.h>

#include <CGAL/Tetrahedral_remeshing/internal/triangulation_3_helpers.h>

#include <boost/unordered_map.hpp>
#include <vector>
#include <limits>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
  template<typename VertexIterator, typename Tr>
  void set_labels_of_incident_cells(VertexIterator begin,
                                    VertexIterator end,
                                    const Tr& tr,
                                    const int& label)
  {
    typedef typename Tr::Cell_handle Cell_handle;
    for (VertexIterator vit = begin; vit != end; ++vit)
    {
      std::vector<Cell_handle> cells;
      tr.finite_incident_cells(*vit, std::back_inserter(cells));

      for (std::size_t i = 0; i < cells.size(); ++i)
        cells[i]->set_subdomain_index(label);
    }
  }

  template<typename VertexIterator>
  void set_dimension(VertexIterator begin,
                     VertexIterator end,
                     const short& dimension)
  {
    for (VertexIterator vit = begin; vit != end; ++vit)
      (*vit)->set_dimension(dimension);
  }

  template<typename PointIterator, typename Tr, typename OutputIterator>
  OutputIterator insert_points(PointIterator begin,
                               PointIterator end,
                               Tr& tr,
                               OutputIterator oit)
  {
    typedef typename Tr::Point Point;

    CGAL_assertion(tr.is_valid());
    int i = 1;
    for (PointIterator pit = begin; pit != end; ++pit, ++i)
    {
      *oit++ = tr.insert(Point(*pit));
    }
    CGAL_assertion(tr.is_valid());

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "(" << i << " points inserted successfully)" << std::endl;
#endif

    return oit;
  }

  template<typename VertexNormalsMap, typename OutputIterator>
  OutputIterator compute_offset_points(const VertexNormalsMap& normals,
                                       const double& offset,
                                       OutputIterator oit)
  {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    std::ofstream ofs("imaginary_points.off");
    ofs << "OFF" << std::endl;
    ofs << normals.size() << " 0 0" << std::endl;
#endif

    for (typename VertexNormalsMap::const_iterator nit = normals.begin();
         nit != normals.end(); ++nit)
    {
      *oit++ = point((*nit).first->point()) + offset * (*nit).second;

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      ofs << ((*nit).first->point() + offset * (*nit).second) << std::endl;
#endif
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    ofs.close();
#endif
    return oit;
  }

  template<typename T3, typename VertexNormalsMap>
  void compute_normals_on_convex_hull(const T3& tr,
                                      VertexNormalsMap& normals)
  {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    std::ofstream ofs("imaginary_normals.xyz");
#endif
    namespace PMP = CGAL::Polygon_mesh_processing;

    typedef typename T3::Geom_traits Gt;
    typedef typename Gt::Vector_3    Vector_3;
    const Gt& gt = tr.geom_traits();

    typedef typename T3::Vertex_handle          Vertex_handle;
    typedef typename T3::Facet                  Facet;
    typedef typename T3::Finite_facets_iterator Finite_facets_iterator;

    boost::unordered_map<Facet, Vector_3> fnormals;
    for (Finite_facets_iterator fit = tr.finite_facets_begin();
         fit != tr.finite_facets_end();
         ++fit)
    {
      Facet f = *fit;
      if (tr.is_infinite(f.first))
      {
        f = tr.mirror_facet(f);
        fnormals.insert(std::make_pair(f, normal(f, gt)));
      }
      else if (tr.is_infinite(f.first->neighbor(f.second)))
        fnormals.insert(std::make_pair(f, normal(f, gt)));
    }

    std::vector<Vertex_handle> vertices;
    tr.finite_adjacent_vertices(tr.infinite_vertex(),
                                std::back_inserter(vertices));
    for (std::size_t i = 0; i < vertices.size(); ++i)
    {
      Vertex_handle vi = vertices[i];
      std::vector<Facet> inc_facets;
      tr.finite_incident_facets(vi, std::back_inserter(inc_facets));

      Vector_3 ni = gt.construct_vector_3_object()(CGAL::NULL_VECTOR);
      for (std::size_t j = 0; j < inc_facets.size(); ++j)
      {
        typename boost::unordered_map<Facet, Vector_3>::iterator
          fnit = fnormals.find(inc_facets[j]);
        if (fnit != fnormals.end())
          ni = gt.construct_sum_of_vectors_3_object()(ni, fnit->second);
        else
        {
          //check for mirror_facet
          fnit = fnormals.find(tr.mirror_facet(inc_facets[j]));
          if (fnit != fnormals.end())
            ni = gt.construct_sum_of_vectors_3_object()(ni, fnit->second);
        }
      }

      if (!typename Gt::Equal_3()(ni, CGAL::NULL_VECTOR))
      {
        ni = gt.construct_opposite_vector_3_object()(ni);
        PMP::internal::normalize(ni, gt);
        normals.insert(std::make_pair(vi, ni));
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        ofs << vi->point() << " " << ni << std::endl;
#endif
      }
    }
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    ofs.close();
#endif
  }

  template<typename Tr>
  double compute_bbox_max_size(const Tr& tr)
  {
    typename Tr::Finite_vertices_iterator vit = tr.finite_vertices_begin();
    const typename Tr::Point& p = vit->point();
    CGAL::Bbox_3 bbox = p.bbox();
    ++vit;
    for ( ; vit != tr.finite_vertices_end(); ++vit)
    {
      bbox = bbox + vit->point().bbox();
    }

    return (std::max)((std::max)(bbox.xmax() - bbox.xmin(),
                                 bbox.ymax() - bbox.ymin()),
                                 bbox.zmax() - bbox.zmin());

  }

  template<typename T3, typename Index/*Subdomain_index*/>
  void add_layer_of_imaginary_tets(T3& tr, const Index& imaginary_index)
  {
    typedef typename T3::Geom_traits Gt;
    typedef typename Gt::Point_3     Point_3;
    typedef typename Gt::Vector_3    Vector_3;

    typedef typename T3::Vertex_handle Vertex_handle;

    //compute normals
    boost::unordered_map<Vertex_handle, Vector_3> normals;
    compute_normals_on_convex_hull(tr, normals);

    //compute bbox max size
    const double& offset = 0.04 * compute_bbox_max_size(tr);

    //compute points to be inserted
    std::vector<Point_3> offset_points;
    compute_offset_points(normals,
      offset,
      std::back_inserter(offset_points));

    //insert vertices on offset
    //note we only need to insert them in the T3, because they
    //are all outside convex hull. The rest of the T3 will not be modified
    std::vector<Vertex_handle> offset_vertices;
    insert_points(offset_points.begin(), offset_points.end(),
      tr, std::back_inserter(offset_vertices));

    CGAL_assertion(tr.is_valid());

    //set labels
    set_labels_of_incident_cells(offset_vertices.begin(),
                                 offset_vertices.end(),
                                 tr,
                                 imaginary_index);

    set_dimension(offset_vertices.begin(), offset_vertices.end(), 3);
  }

}//end namespace internal
}//end namespace Tetrahedral_remeshing
}//end namesapce CGAL

#endif
