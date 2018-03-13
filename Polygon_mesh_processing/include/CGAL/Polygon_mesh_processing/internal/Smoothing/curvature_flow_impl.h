// Copyright (c) 2017 GeometryFactory (France).
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
// Author(s)     : Konstantinos Katrioplas (konst.katrioplas@gmail.com)

#ifndef CGAL_POLYGON_MESH_PROCESSING_CURVATURE_FLOW_IMPL_H
#define CGAL_POLYGON_MESH_PROCESSING_CURVATURE_FLOW_IMPL_H

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/Weights.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <CGAL/Eigen_solver_traits.h>
#include <CGAL/Eigen_matrix.h>
#include <CGAL/barycenter.h>
#include <Eigen/Sparse>

#include <fstream>
#include <unordered_set>
#include <unordered_map>


#include <CGAL/Polygon_mesh_processing/internal/Smoothing/constraints_map.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<typename PolygonMesh,
         typename VertexPointMap,
         typename VertexConstraintMap,
         typename GeomTraits>
class Shape_smoother{

private:

  typedef typename GeomTraits::FT NT;
  typedef typename GeomTraits::Point_3 Point;
  typedef typename GeomTraits::Triangle_3 Triangle;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;

  typedef typename boost::property_map<PolygonMesh, boost::vertex_index_t>::type IndexMap;

  typedef typename Eigen::SparseMatrix<double> Eigen_matrix;
  typedef typename Eigen::VectorXd Eigen_vector;

  typedef typename Eigen::BiCGSTAB<Eigen_matrix, Eigen::IncompleteLUT<double> > Eigen_solver;
  // typedef typename Eigen::SimplicialLDLT<Eigen_matrix> Eigen_solver;

public:
  Shape_smoother(PolygonMesh& mesh,
                 VertexPointMap& vpmap,
                 VertexConstraintMap& vcmap) :
      mesh_(mesh),
      vpmap_(vpmap),
      vcmap_(vcmap),
      weight_calculator_(mesh, vpmap),
      inc_areas_calculator_(mesh),
      nb_vert_(static_cast<int>(vertices(mesh).size())) {}

  template<typename FaceRange>
  void init_smoothing(const FaceRange& face_range)
  {
    check_face_range(face_range);
  }

  void setup_system(Eigen_matrix& A, const Eigen_matrix& L, Eigen_matrix& D,
                    Eigen_vector& bx, Eigen_vector& by, Eigen_vector& bz,
                    const double& time)
  {
    calc_mass_matrix(D);

    compute_coeff_matrix(A, L, D, time);

    compute_rhs(bx, by, bz, D);

  }

  void solve_system(const Eigen_matrix& A,
                    Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz,
                    const Eigen_vector& bx, const Eigen_vector& by, const Eigen_vector& bz)
  {

    Eigen_solver solver;
    solver.compute(A);


    Xx = solver.solve(bx);
    Xy = solver.solve(by);
    Xz = solver.solve(bz);



    if(solver.info() != Eigen::Success)
    {
      std::cerr << "Not Solved!" << std::endl;
      return;
    }


    /*
    std::cout << "Xx old:\n";
    for(int i = 0; i < Xx.rows(); ++i)
    {
      std::cout << Xx[i] << std::endl;
    }
    std::cout <<"\n";
    */
  }

  void calc_stiff_matrix(Eigen_matrix& mat)
  {
    typedef Eigen::Triplet<double> Triplet;
    std::vector<Triplet> tripletList;
    tripletList.reserve(8 * nb_vert_);
    // todo: calculate exactly how many non zero entries there will be.

    for(face_descriptor f : frange_)
    {
      for(halfedge_descriptor hi : halfedges_around_face(halfedge(f, mesh_), mesh_))
      {
        vertex_descriptor v_source = source(hi, mesh_);
        vertex_descriptor v_target = target(hi, mesh_);

        if(!is_constrained(v_source) && !is_constrained(v_target))
        {
          auto i_source = vimap_[v_source];
          auto i_target = vimap_[v_target];
          NT Lij = weight_calculator_(hi);
          tripletList.push_back(Triplet(i_source, i_target, Lij));
          tripletList.push_back(Triplet(i_target, i_source, Lij));
          tripletList.push_back(Triplet(i_source, i_source, -Lij));
          tripletList.push_back(Triplet(i_target, i_target, -Lij));
        }
      }
    }
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
  }

  void update_mesh(Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz)
  {
    update_map(Xx, Xy, Xz);

    // normalize area
    //NT surface_area = area(faces(mesh_), mesh_);
    //std::cout << "area= " << surface_area << nl;
    //std::cout << "area sqrt= " << CGAL::sqrt(surface_area) << nl;
    //normalize_area(Xx, Xy, Xz);
    //update_map(Xx, Xy, Xz);
    //surface_area = area(faces(mesh_), mesh_);
    //std::cout << "surface_area normalized= " << surface_area << nl;
  }

private:
  void calc_mass_matrix(Eigen_matrix& D)
  {
    for(face_descriptor f : frange_)
    {
      double area = face_area(f, mesh_);
      for(vertex_descriptor v : vertices_around_face(halfedge(f, mesh_), mesh_))
      {
        auto idx = vimap_[v];
        if(!is_constrained(v))
          D.coeffRef(idx, idx) += 2.0 * area;
        else
          D.coeffRef(idx, idx) = 1.0;
      }
    }
    D /= 12.0;
  }

  void compute_coeff_matrix(Eigen_matrix& A, const Eigen_matrix& L, const Eigen_matrix& D, const double& time)
  {
    assert(A.rows() != 0);
    assert(A.cols() != 0);
    assert(A.rows() == L.rows());
    assert(A.cols() == L.cols());
    assert(A.rows() == D.rows());
    assert(A.cols() == D.cols());
    A = D - time * L;
  }

  void compute_rhs(Eigen_vector& bx, Eigen_vector& by, Eigen_vector& bz,
                   Eigen_matrix& D)
  {
    for(vertex_descriptor vi : vrange_)
    {
      int index = vimap_[vi];
      Point p = get(vpmap_, vi);
      bx.coeffRef(index) = p.x();
      by.coeffRef(index) = p.y();
      bz.coeffRef(index) = p.z();
    }
    bx = D * bx;
    by = D * by;
    bz = D * bz;
  }

  // update mesh
  // -----------------------------------
  Triangle triangle(face_descriptor f) const
  {
    halfedge_descriptor h = halfedge(f, mesh_);
    vertex_descriptor v1  = target(h, mesh_);
    vertex_descriptor v2  = target(next(h, mesh_), mesh_);
    vertex_descriptor v3  = target(next(next(h, mesh_), mesh_), mesh_);
    return Triangle(get(vpmap_, v1), get(vpmap_, v2), get(vpmap_, v3));
  }

  /*
  void normalize_area(Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz)
  {
    NT surface_area = area(faces(mesh_), mesh_);
    Xx /= CGAL::sqrt(surface_area);
    Xy /= CGAL::sqrt(surface_area);
    Xz /= CGAL::sqrt(surface_area);
  }
  */

  void update_map(Eigen_vector& Xx, Eigen_vector& Xy, Eigen_vector& Xz)
  {
    for (vertex_descriptor v : vertices(mesh_))
    {
      int index = get(vimap_, v);
      NT x_new = Xx.coeffRef(index);
      NT y_new = Xy.coeffRef(index);
      NT z_new = Xz.coeffRef(index);
      put(vpmap_, v, Point(x_new, y_new, z_new));
    }
  }

  // handling constrains
  // -----------------------------------
  bool is_constrained(const vertex_descriptor& v)
  {
    return get(vcmap_, v);
  }

  template<typename FaceRange>
  void check_face_range(const FaceRange& face_range)
  {
    BOOST_FOREACH(face_descriptor f, face_range)
    {
      frange_.insert(f);
      BOOST_FOREACH(vertex_descriptor v, vertices_around_face(halfedge(f, mesh_), mesh_))
        vrange_.insert(v);
    }
  }

  private:
  // data members
  // ------------
  PolygonMesh& mesh_;
  VertexPointMap& vpmap_;
  std::size_t nb_vert_;
  std::set<face_descriptor> frange_;
  std::set<vertex_descriptor> vrange_;
  IndexMap vimap_ = get(boost::vertex_index, mesh_);
  VertexConstraintMap vcmap_;
  Edge_cotangent_weight<PolygonMesh, VertexPointMap> weight_calculator_;
  Incident_area<PolygonMesh> inc_areas_calculator_;

};


} // internal
} // PMP
} // CGAL


#endif // CGAL_POLYGON_MESH_PROCESSING_CURVATURE_FLOW_IMPL_H
