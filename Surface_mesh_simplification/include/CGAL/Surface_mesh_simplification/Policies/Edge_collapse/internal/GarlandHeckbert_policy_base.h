// Copyright (c) 2019  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Baskin Burak Senbaslar,
//                 Mael Rouxel-Labbé,
//                 Julian Komaromy

#ifndef CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_INTERNAL_GARLANDHECKBERT_POLICIES_BASE_H
#define CGAL_SURFACE_MESH_SIMPLIFICATION_POLICIES_INTERNAL_GARLANDHECKBERT_POLICIES_BASE_H

#include <CGAL/license/Surface_mesh_simplification.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/internal/Common.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/internal/GarlandHeckbert_core.h>

#include <iostream>
#include <CGAL/tags.h>

#include <Eigen/Dense>

#include <boost/optional/optional.hpp>

namespace CGAL {
namespace Surface_mesh_simplification {
namespace internal {
  
template<typename Mat_4>
Mat_4 combine_matrices(const Mat_4& a, const Mat_4& b) {
  return a + b;    
}

template<typename VCM, 
         typename GeomTraits, 
         typename QuadricImpl>
class GarlandHeckbert_cost_base
{
  public:

  // Tells the main function of 'Edge_collapse' that these
  // policies must call "initialize" and "update" functions.
  typedef CGAL::Tag_true Update_tag;  

  typedef VCM Vertex_cost_map;
  typedef typename GeomTraits::FT FT;
  
  typedef Eigen::Matrix<FT, 4, 4, Eigen::DontAlign> Mat_4;
  typedef Eigen::Matrix<FT, 4, 1> Col_4;
  
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Vector_3 Vector_3;

  GarlandHeckbert_cost_base() : m_cost_matrices() { } 
  
  GarlandHeckbert_cost_base(FT dm) 
    : m_cost_matrices(), discontinuity_multiplier(dm) { } 
  
  GarlandHeckbert_cost_base(Vertex_cost_map vcm, FT dm) 
    : m_cost_matrices(vcm), discontinuity_multiplier(dm) { } 

  template<typename TM>
  static bool is_discontinuity_edge(const typename boost::graph_traits<TM>::halfedge_descriptor& h, 
                                    const TM& tmesh)
  {
    return is_border_edge(h, tmesh);
  }

  // initialize all quadrics 
  template<typename TM,
           typename VPM> 
  void initialize(const TM& tmesh, const VPM& vpm, const GeomTraits& gt) const
  {
    // the graph library types used here
    typedef boost::graph_traits<TM> GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
    typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;
    typedef typename GraphTraits::face_descriptor face_descriptor;
    typedef typename boost::property_traits<VPM>::reference Point_reference;
    
    Mat_4 zero_mat = Mat_4::Zero();
    
    for(vertex_descriptor v : vertices(tmesh))
    {
      put(m_cost_matrices, v, zero_mat);
    }

    for(face_descriptor f : faces(tmesh))
    {
      if (f == GraphTraits::null_face())
      {
        continue;
      }

      const halfedge_descriptor h = halfedge(f, tmesh);

      // get the three vertices of the triangular face
      const Point_reference p = get(vpm, source(h, tmesh));
      const Point_reference q = get(vpm, target(h, tmesh));
      const Point_reference r = get(vpm, target(next(h, tmesh), tmesh));

      const Vector_3 face_normal = unit_normal(p, q, r);
      
      // construtct the (4 x 4) matrix representing the plane quadric
      const Mat_4 quadric = construct_quadric(face_normal, r, gt);

      for(halfedge_descriptor shd : halfedges_around_face(h, tmesh))
      {
        const vertex_descriptor vs = source(shd, tmesh);
        const vertex_descriptor vt = target(shd, tmesh);

        put(m_cost_matrices, vt, combine_matrices(get(m_cost_matrices, vt), quadric));

        if(!is_discontinuity_edge(shd, tmesh))
        {
          continue;
        }

        //TODO behaviour at discontinuous edges is largely untested at this stage
        //TODO helper function to add to target
        const Vector_3 edge_vector = Vector_3(get(vpm, vs), get(vpm, vt));
        const Vector_3 discontinuity_normal = cross_product(edge_vector, face_normal); 
        
        // normalize  
        const Vector_3 normalized = discontinuity_normal 
          / sqrt(discontinuity_normal.squared_length());
        const Mat_4 discontinuous_quadric 
          = discontinuity_multiplier * construct_quadric(normalized, get(vpm, vs), gt);
        
        put(m_cost_matrices, vs, combine_matrices(get(m_cost_matrices, vs), discontinuous_quadric));
        put(m_cost_matrices, vt, combine_matrices(get(m_cost_matrices, vt), discontinuous_quadric));
      }
    }
  }

  template<typename Profile>
  boost::optional<typename Profile::FT> operator()(const Profile& profile,
      const boost::optional<typename Profile::Point>& placement) const 
  {
    typedef boost::optional<typename Profile::FT> Optional_FT;

    if(!placement) {
      // return empty value
      return boost::optional<typename Profile::FT>();
    }

    CGAL_precondition(!get(m_cost_matrices, profile.v0()).isZero(0));
    CGAL_precondition(!get(m_cost_matrices, profile.v1()).isZero(0));

    const Mat_4 combined_matrix = combine_matrices(get(m_cost_matrices, profile.v0()),
                                                   get(m_cost_matrices, profile.v1()));
    const Col_4 pt = point_to_homogenous_column(*placement);
    const Optional_FT cost = (pt.transpose() * combined_matrix * pt)(0, 0);

    return cost;
  }
  
  template<typename Profile, typename VertexDescriptor>
  void update_after_collapse(const Profile& profile,
      const VertexDescriptor new_v) const
  {
    put(m_cost_matrices, new_v,
        combine_matrices(get(m_cost_matrices, profile.v0()),
                         get(m_cost_matrices, profile.v1())));
  }
  
  protected:
  void init_vcm(const Vertex_cost_map& vcm) {
    m_cost_matrices = vcm;
  }
  
  private: 
  Vertex_cost_map m_cost_matrices;
  
  FT discontinuity_multiplier;
  
  Col_4 point_to_homogenous_column(const Point_3& point) const {
    return Col_4(point.x(), point.y(), point.z(), FT(1));
  } 
  
  // use CRTP to call the quadric implementation 
  Mat_4 construct_quadric(const Vector_3& normal, const Point_3& point, const GeomTraits& gt) const {
    return static_cast<const QuadricImpl*>(this)->construct_quadric_from_normal(normal, point, gt);
  }
};

//TODO should probably take a TirangleMesh as Parameter and define the VCM type in base classes
template<typename VCM, 
         typename GeomTraits, 
         typename QuadricImpl>
class GarlandHeckbert_placement_base
{
  public:
  // define type required by the Get_cost concept
  typedef VCM Vertex_cost_map;
  
  // matrix and column vector types
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef Eigen::Matrix<FT, 4, 4, Eigen::DontAlign> Mat_4;
  typedef Eigen::Matrix<FT, 4, 1> Col_4;

  GarlandHeckbert_placement_base() { }
  GarlandHeckbert_placement_base(Vertex_cost_map cost_matrices) : m_cost_matrices(cost_matrices) { }

  template<typename Profile>
  boost::optional<typename Profile::Point> operator()(const Profile& profile) const
  {
    CGAL_precondition(!get(m_cost_matrices, profile.v0()).isZero(0));
    CGAL_precondition(!get(m_cost_matrices, profile.v1()).isZero(0));

    // the combined matrix has already been computed in the evaluation of the cost...
    const Mat_4 combinedMatrix = combine_matrices(
        get(m_cost_matrices, profile.v0()),
        get(m_cost_matrices, profile.v1()));

    const Col_4 p0 = point_to_homogenous_column(profile.p0());
    const Col_4 p1 = point_to_homogenous_column(profile.p1());

    const Col_4 opt = construct_optimum(combinedMatrix, p0, p1);
    
    boost::optional<typename Profile::Point> pt = typename Profile::Point(
        opt(0) / opt(3),
        opt(1) / opt(3),
        opt(2) / opt(3));

    return pt;
  }
  
  protected:
  void init_vcm(const Vertex_cost_map& vcm) {
    m_cost_matrices = vcm;
  }

  private:
  Vertex_cost_map m_cost_matrices;
    
  // use CRTP to call the quadric implementation 
  Col_4 construct_optimum(const Mat_4& mat, const Col_4& p0, const Col_4& p1) const {
    return static_cast<const QuadricImpl*>(this)->construct_optimal_point(mat, p0, p1);
  }
  
  Col_4 point_to_homogenous_column(const Point_3& point) const {
    return Col_4(point.x(), point.y(), point.z(), FT(1));
  } 
};

} // namespace internal
} // namespace Surface_mesh_simplification
} // namespace CGAL
#endif
