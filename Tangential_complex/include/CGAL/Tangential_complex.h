// Copyright (c) 2014  INRIA Sophia-Antipolis (France)
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Clement Jamin


#ifndef TANGENTIAL_COMPLEX_H
#define TANGENTIAL_COMPLEX_H

#include <CGAL/Tangential_complex/config.h>

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Dimension.h>

#include <CGAL/Epick_d.h>
#include <CGAL/Regular_triangulation_euclidean_traits.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/Tangential_complex/utilities.h>
#include <CGAL/Tangential_complex/Point_cloud.h>

#ifdef CGAL_TC_PROFILING
# include <CGAL/Mesh_3/Profiling_tools.h>
#endif

#include <CGAL/IO/Triangulation_off_ostream.h> // CJTODO TEMP

#include <Eigen/Core>
#include <Eigen/Eigen>

#include <boost/iterator/transform_iterator.hpp>

#include <vector>
#include <utility>
#include <sstream>
#include <iostream>
#include <limits>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/parallel_for.h>
#endif

namespace CGAL {

using namespace Tangential_complex_;
  
/// The class Tangential_complex represents a tangential complex
template <
  typename Kernel, 
  int Intrinsic_dimension,
  typename Concurrency_tag = CGAL::Parallel_tag,
  typename Tr = Regular_triangulation
  <
    Regular_triangulation_euclidean_traits<
      Epick_d<Dimension_tag<Intrinsic_dimension> > >,

    Triangulation_data_structure
    <
      typename Regular_triangulation_euclidean_traits<
        Epick_d<Dimension_tag<Intrinsic_dimension> > >::Dimension,
      Triangulation_vertex<Regular_triangulation_euclidean_traits<
        Epick_d<Dimension_tag<Intrinsic_dimension> > >, std::size_t >,
      Triangulation_full_cell<Regular_triangulation_euclidean_traits<
        Epick_d<Dimension_tag<Intrinsic_dimension> > > >
    >
  >
>
class Tangential_complex
{
  typedef typename Kernel::FT                         FT;
  typedef typename Kernel::Point_d                    Point;
  typedef typename Kernel::Vector_d                   Vector;

  typedef Tr                                          Triangulation;
  typedef typename Triangulation::Geom_traits         Tr_traits;
  typedef typename Triangulation::Point               Tr_point;
  typedef typename Triangulation::Bare_point          Tr_bare_point;
  typedef typename Triangulation::Vertex_handle       Tr_vertex_handle;
  typedef typename Triangulation::Full_cell_handle    Tr_full_cell_handle;
  
  typedef std::vector<Vector>                         Tangent_space_basis;

  typedef std::vector<Point>                          Points;
  typedef Point_cloud_data_structure<Points>          Points_ds;
  typedef typename Points_ds::KNS_range               KNS_range;
  typedef typename Points_ds::KNS_iterator            KNS_iterator;
  typedef typename Points_ds::INS_range               INS_range;
  typedef typename Points_ds::INS_iterator            INS_iterator;

  // Store a local triangulation and a handle to its center vertex
  struct Tr_and_VH
  {
  public:
    Tr_and_VH()
      : m_tr(NULL) {}
    Tr_and_VH(int dim)
      : m_tr(new Triangulation(dim)) {}
    
    ~Tr_and_VH() { delete m_tr; }

    Triangulation & construct_triangulation(int dim)
    { 
      m_tr = new Triangulation(dim);
      return tr();
    }

    Triangulation &      tr()       { return *m_tr; }
    Triangulation const& tr() const { return *m_tr; }
    

    Tr_vertex_handle const& center_vertex() const { return m_center_vertex; }
    Tr_vertex_handle & center_vertex() { return m_center_vertex; }

  private:
    Triangulation* m_tr;
    Tr_vertex_handle m_center_vertex;
  };

  typedef typename std::vector<Tr_and_VH>             Tr_container;
  typedef typename std::vector<Tangent_space_basis>   TS_container;

public:
  /// Constructor
  Tangential_complex(const Kernel &k = Kernel())
  : m_k(k){}
  
  /// Constructor for a range of points
  template <typename InputIterator>
  Tangential_complex(InputIterator first, InputIterator last, 
                     const Kernel &k = Kernel())
  : m_k(k), m_points(first, last), m_points_ds(m_points, k) {}

  /// Destructor
  ~Tangential_complex() {}

  void compute_tangential_complex()
  {
#ifdef CGAL_TC_PROFILING
    Wall_clock_timer t;
#endif

    // We need to do that because we don't want the container to copy the
    // already-computed triangulations (while resizing) since it would
    // invalidate the vertex handles stored beside the triangulations
    m_triangulations.resize(m_points.size());
    m_tangent_spaces.resize(m_points.size());
    
#ifdef CGAL_LINKED_WITH_TBB
    // Parallel
    if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
    {
      // Apply moves in triangulation
      tbb::parallel_for(tbb::blocked_range<size_t>(0, m_points.size()),
        Compute_tangent_triangulation(*this)
      );
    }
    // Sequential
    else
#endif // CGAL_LINKED_WITH_TBB
    {
      for (std::size_t i = 0 ; i < m_points.size() ; ++i)
        compute_tangent_triangulation(i);
    }
    
#ifdef CGAL_TC_PROFILING
    std::cerr << "Tangential complex computed in " << t.elapsed() 
              << " seconds." << std::endl;
#endif
  }

  std::ostream &export_to_off(std::ostream & os, 
                              bool color_inconsistencies = false)
  {
    const int ambient_dim = Ambient_dimension<Point>::value;
    if (ambient_dim < 2)
    {
      std::cerr << "Error: export_to_off => ambient dimension should be >= 2."
                << std::endl;
      os << "Error: export_to_off => ambient dimension should be >= 2." 
         << std::endl;
      return os;
    }
    if (ambient_dim > 3)
    {
      std::cerr << "Warning: export_to_off => ambient dimension should be "
                   "<= 3. Only the first 3 coordinates will be exported." 
                << std::endl;
    }

    int num_coords = min(ambient_dim, 3);

    if (Intrinsic_dimension < 1 || Intrinsic_dimension > 3)
    {
      std::cerr << "Error: export_to_off => intrinsic dimension should be "
                   "between 1 and 3."
                << std::endl;
      os << "Error: export_to_off => intrinsic dimension should be "
            "between 1 and 3."
         << std::endl;
      return os;
    }

    std::stringstream output;

    //******** VERTICES ************

    Points::const_iterator it_p = m_points.begin();
    Points::const_iterator it_p_end = m_points.end();
    // For each point p
    for ( ; it_p != it_p_end ; ++it_p)
    {
      int i = 0;
      for ( ; i < num_coords ; ++i)
        output << (*it_p)[i] << " ";
      if (i == 2)
        output << "0";
      output << std::endl;
    }

    //******** CELLS ************

    std::size_t num_cells = 0;
    Tr_container::const_iterator it_tr = m_triangulations.begin();
    Tr_container::const_iterator it_tr_end = m_triangulations.end();
    // For each triangulation
    for ( ; it_tr != it_tr_end ; ++it_tr)
    {
      Triangulation const& tr    = it_tr->tr();
      Tr_vertex_handle center_vh = it_tr->center_vertex();

      std::vector<Tr_full_cell_handle> incident_cells;
      tr.incident_full_cells(center_vh, std::back_inserter(incident_cells));

      std::vector<Tr_full_cell_handle>::const_iterator it_c = incident_cells.begin();
      std::vector<Tr_full_cell_handle>::const_iterator it_c_end= incident_cells.end();
      // For each cell
      for ( ; it_c != it_c_end ; ++it_c)
      {
        output << Intrinsic_dimension + 1 << " ";
        
        if (color_inconsistencies)
        {
          std::set<std::size_t> c;
          for (int i = 0 ; i < Intrinsic_dimension + 1 ; ++i)
          {
            std::size_t data = (*it_c)->vertex(i)->data();
            output << data << " ";
            c.insert(data);
          }
          if (is_simplex_consistent(c))
            output << "200 200 200";
          else
            output << "255 0 0";
        }
        else
        {
          for (int i = 0 ; i < Intrinsic_dimension + 1 ; ++i)
            output << (*it_c)->vertex(i)->data() << " ";
        }

        output << std::endl;
        ++num_cells;
      }
    }

    os << "OFF \n"
       << m_points.size() << " " 
       << num_cells << " "
       << "0 \n"
       << output.str();

    return os;
  }

private:

  class Compare_distance_to_ref_point
  {
  public:
    Compare_distance_to_ref_point(Point const& ref, Kernel const& k)
      : m_ref(ref), m_k(k) {}

    bool operator()(Point const& p1, Point const& p2)
    {
      Kernel::Squared_distance_d sqdist = m_k.squared_distance_d_object();
      return sqdist(p1, m_ref) < sqdist(p2, m_ref);
    }

  private:
    Point const& m_ref;
    Kernel const& m_k;
  };

  struct Tr_vertex_to_global_point
  {
    typedef typename Tr_vertex_handle	argument_type;
    typedef typename Point	result_type;

    Tr_vertex_to_global_point(Points const& points)
      : m_points(points) {}

    result_type operator()(argument_type const& vh) const
    {
      return m_points[vh->data()];
    }

  private:
    Points const& m_points;
  };

  struct Tr_vertex_to_bare_point
  {
    typedef typename Tr_vertex_handle	argument_type;
    typedef typename Tr_bare_point	result_type;

    Tr_vertex_to_bare_point(Tr_traits const& traits)
      : m_traits(traits) {}

    result_type operator()(argument_type const& vh) const
    {
      typename Tr_traits::Point_drop_weight_d pdw =
        m_traits.point_drop_weight_d_object();
      return pdw(vh->point());
    }

  private:
    Tr_traits const& m_traits;
  };

#ifdef CGAL_LINKED_WITH_TBB
  // Functor for compute_tangential_complex function
  class Compute_tangent_triangulation
  {
    Tangential_complex & m_tc;

  public:
    // Constructor
    Compute_tangent_triangulation(Tangential_complex &tc)
    : m_tc(tc)
    {}

    // Constructor
    Compute_tangent_triangulation(const Compute_tangent_triangulation &ctt)
    : m_tc(ctt.m_tc)
    {}

    // operator()
    void operator()( const tbb::blocked_range<size_t>& r ) const
    {
      for( size_t i = r.begin() ; i != r.end() ; ++i)
        m_tc.compute_tangent_triangulation(i);
    }
  };
#endif // CGAL_LINKED_WITH_TBB

  void compute_tangent_triangulation(std::size_t i)
  {
    //std::cerr << "***********************************************" << std::endl;
    Triangulation &local_tr =
      m_triangulations[i].construct_triangulation(Intrinsic_dimension);
    const Tr_traits &local_tr_traits = local_tr.geom_traits();
    Tr_vertex_handle &center_vertex = m_triangulations[i].center_vertex();

    // Kernel functor & objects
    Kernel::Difference_of_points_d k_diff_pts =
      m_k.difference_of_points_d_object();
    Kernel::Squared_distance_d k_sqdist = 
      m_k.squared_distance_d_object();

    // Triangulation's traits functor & objects
    Tr_traits::Squared_distance_d sqdist = 
      local_tr_traits.squared_distance_d_object();
    Tr_traits::Point_drop_weight_d drop_w = 
      local_tr_traits.point_drop_weight_d_object();
    Tr_traits::Center_of_sphere_d center_of_sphere = 
      local_tr_traits.center_of_sphere_d_object();

    // Estimate the tangent space
    const Point &center_pt = m_points[i];
    m_tangent_spaces[i] = compute_tangent_space(center_pt);
      
    //***************************************************
    // Build a minimal triangulation in the tangent space
    // (we only need the star of p)
    //***************************************************

    // Insert p
    Tr_point wp = local_tr_traits.construct_weighted_point_d_object()(
      local_tr_traits.construct_point_d_object()(0, 0),
      0);
    center_vertex = local_tr.insert(wp);
    center_vertex->data() = i;

    //const int NUM_NEIGHBORS = 150;
    //KNS_range ins_range = m_points_ds.query_ANN(center_pt, NUM_NEIGHBORS);
    INS_range ins_range = m_points_ds.query_incremental_ANN(center_pt);
    
    // While building the local triangulation, we keep the radius
    // of the sphere "star sphere" centered at "center_vertex" 
    // and which contains all the
    // circumspheres of the star of "center_vertex"
    FT star_sphere_squared_radius = std::numeric_limits<FT>::max();

    // Insert points until we find a point which is outside "star shere"
    for (INS_iterator nn_it = ins_range.begin() ; 
         nn_it != ins_range.end() ; 
         ++nn_it)
    {
      std::size_t neighbor_point_idx = nn_it->first;

      // ith point = p, which is already inserted
      if (neighbor_point_idx != i)
      {
        const Point &neighbor_pt = m_points[neighbor_point_idx];

        if (k_sqdist(center_pt, neighbor_pt) > star_sphere_squared_radius)
          break;

        Tr_point proj_pt = project_point_and_compute_weight(
          neighbor_pt, center_pt, m_tangent_spaces[i]);

        FT squared_dist_to_tangent_plane = 
          local_tr_traits.point_weight_d_object()(proj_pt);
        FT w =  -squared_dist_to_tangent_plane;
        Tr_point wp = local_tr_traits.construct_weighted_point_d_object()(
          drop_w(proj_pt),
          w);
          
        Tr_vertex_handle vh = local_tr.insert_if_in_star(wp, center_vertex);
        //Tr_vertex_handle vh = local_tr.insert(wp);
        if (vh != Tr_vertex_handle())
        {
          vh->data() = neighbor_point_idx;

          // Let's recompute star_sphere_squared_radius
          if (local_tr.current_dimension() >= Intrinsic_dimension)
          {
            star_sphere_squared_radius = 0.;
            // Get the incident cells and look for the biggest circumsphere
            std::vector<Tr_full_cell_handle> incident_cells;
            local_tr.incident_full_cells(
              center_vertex, 
              std::back_inserter(incident_cells));
            for (auto cell : incident_cells) // CJTODO C++11
            {
              if (local_tr.is_infinite(cell))
              {
                star_sphere_squared_radius = std::numeric_limits<FT>::max();
                break;
              }
              else
              {
                //*********************************
                // We don't compute the circumsphere of the simplex in the
                // local tangent plane since it would involve to take the 
                // weights of the points into account later
                // (which is a problem since the ANN is performed on the
                // points in the ambient dimension)
                // Instead, we compute the subspace defined by the simplex
                // and we compute the circumsphere in this subspace
                // and we extract the diameter
                Tangent_space_basis tsb;
                tsb.reserve(Intrinsic_dimension);
                Point const& orig = m_points[cell->vertex(0)->data()];
                for (int ii = 1 ; ii <= Intrinsic_dimension ; ++ii)
                { 
                  tsb.push_back(k_diff_pts(
                    m_points[cell->vertex(ii)->data()], orig));
                }
                tsb = compute_gram_schmidt_basis(tsb, m_k);

                std::vector<Tr_bare_point> proj_pts;
                proj_pts.reserve(Intrinsic_dimension + 1);
                // For each point p
                for (int ii = 0 ; ii <= Intrinsic_dimension ; ++ii)
                {
                  proj_pts.push_back(project_point(
                    m_points[cell->vertex(ii)->data()], orig, tsb));
                }
                
                Tr_bare_point c = center_of_sphere(
                  proj_pts.begin(), proj_pts.end());

                FT sq_circumdiam = 4.*sqdist(c, proj_pts[0]);
                if (sq_circumdiam > star_sphere_squared_radius)
                  star_sphere_squared_radius = sq_circumdiam;
              }
            }
          }
        }
        //std::cerr << star_sphere_squared_radius << std::endl;
      }
    }


    // CJTODO DEBUG
    //std::cerr << "\nChecking topology and geometry..."
    //          << (local_tr.is_valid(true) ? "OK.\n" : "Error.\n");
    // DEBUG: output the local mesh into an OFF file
    //std::stringstream sstr;
    //sstr << "data/local_tri_" << i << ".off";
    //std::ofstream off_stream_tr(sstr.str());
    //CGAL::export_triangulation_to_off(off_stream_tr, local_tr);
  }

  Tangent_space_basis compute_tangent_space(const Point &p) const
  {
    // Kernel functors
    Kernel::Construct_vector_d      constr_vec = m_k.construct_vector_d_object();
    Kernel::Squared_length_d        sqlen      = m_k.squared_length_d_object();
    Kernel::Scaled_vector_d         scaled_vec = m_k.scaled_vector_d_object();
    Kernel::Scalar_product_d        inner_pdct = m_k.scalar_product_d_object();
    Kernel::Difference_of_vectors_d diff_vec   = m_k.difference_of_vectors_d_object();

    KNS_range kns_range = m_points_ds.query_ANN(
      p, NUM_POINTS_FOR_PCA, false);

    //******************************* PCA *************************************

    const int amb_dim = Ambient_dimension<Point>::value;
    // One row = one point
    Eigen::MatrixXd mat_points(NUM_POINTS_FOR_PCA, amb_dim);
    KNS_iterator nn_it = kns_range.begin();
    for (int j = 0 ; 
         j < NUM_POINTS_FOR_PCA && nn_it != kns_range.end() ; 
         ++j, ++nn_it)
    {
      for (int i = 0 ; i < amb_dim ; ++i)
        mat_points(j, i) = m_points[nn_it->first][i]; // CJTODO: Use kernel functor
    }
    Eigen::MatrixXd centered = mat_points.rowwise() - mat_points.colwise().mean();
    Eigen::MatrixXd cov = centered.adjoint() * centered;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);

    // The eigenvectors are sorted in increasing order of their corresponding
    // eigenvalues
    Tangent_space_basis ts;
    for (int i = amb_dim - 1 ; i >= amb_dim - Intrinsic_dimension ; --i)
    {
      ts.push_back(constr_vec(
        amb_dim, 
        eig.eigenvectors().col(i).data(), 
        eig.eigenvectors().col(i).data() + amb_dim));
    }

    //*************************************************************************

    //Vector n = m_k.point_to_vector_d_object()(p);
    //n = scaled_vec(n, 1./sqrt(sqlen(n)));
    //std::cerr << "IP = " << inner_pdct(n, ts[0]) << " & " << inner_pdct(n, ts[1]) << std::endl;

    return compute_gram_schmidt_basis(ts, m_k);

    /*
    // CJTODO: this is only for a sphere in R^3
    Vector t1(-p[1] - p[2], p[0], p[0]);
    Vector t2(p[1] * t1[2] - p[2] * t1[1],
              p[2] * t1[0] - p[0] * t1[2],
              p[0] * t1[1] - p[1] * t1[0]);
    
    // Normalize t1 and t2
    Kernel::Scaled_vector_d scale = m_k.scaled_vector_d_object();

    Tangent_space_basis ts;
    ts.reserve(Intrinsic_dimension);
    ts.push_back(scale(t1, 1./CGAL::sqrt(sqlen(t1))));
    ts.push_back(scale(t2, 1./CGAL::sqrt(sqlen(t2))));

    return ts;

    // Alternative code (to be used later)
    //Vector n = m_k.point_to_vector_d_object()(p);
    //n = scaled_vec(n, 1./sqrt(sqlen(n)));
    //Vector t1(12., 15., 65.);
    //Vector t2(32., 5., 85.);
    //Tangent_space_basis ts;
    //ts.reserve(Intrinsic_dimension);
    //ts.push_back(diff_vec(t1, scaled_vec(n, inner_pdct(t1, n))));
    //ts.push_back(diff_vec(t2, scaled_vec(n, inner_pdct(t2, n))));
    //return compute_gram_schmidt_basis(ts, m_k);
    */
  }
  
  // Project the point in the tangent space
  // The weight will be the squared distance between p and the projection of p
  Tr_bare_point project_point(const Point &p, const Point &origin, 
                         const Tangent_space_basis &ts) const
  {
    Kernel::Scalar_product_d inner_pdct = m_k.scalar_product_d_object();
    Kernel::Difference_of_points_d diff_points =
      m_k.difference_of_points_d_object();

    std::vector<FT> coords;
    // Ambiant-space coords of the projected point
    coords.reserve(Intrinsic_dimension);
    for (std::size_t i = 0 ; i < Intrinsic_dimension ; ++i)
    {
      // Compute the inner product p * ts[i]
      Vector v = diff_points(p, origin);
      FT coord = inner_pdct(v, ts[i]);
      coords.push_back(coord);
    }

    return Tr_bare_point(Intrinsic_dimension, coords.begin(), coords.end());
  }

  // Project the point in the tangent space
  // The weight will be the squared distance between p and the projection of p
  Tr_point project_point_and_compute_weight(
    const Point &p, const Point &origin, const Tangent_space_basis &ts) const
  {
    Kernel::Scalar_product_d inner_pdct = m_k.scalar_product_d_object();
    Kernel::Difference_of_points_d diff_points =
      m_k.difference_of_points_d_object();

    std::vector<FT> coords;
    // Ambiant-space coords of the projected point
    std::vector<FT> p_proj(origin.cartesian_begin(), origin.cartesian_end()); // CJTODO: use kernel functors?
    coords.reserve(Intrinsic_dimension);
    for (std::size_t i = 0 ; i < Intrinsic_dimension ; ++i)
    {
      // Compute the inner product p * ts[i]
      Vector v = diff_points(p, origin);
      FT coord = inner_pdct(v, ts[i]);
      coords.push_back(coord);

      // p_proj += coord * v;
      for (int j = 0 ; j < Ambient_dimension<Point>::value ; ++j)
        p_proj[i] += coord * ts[i][j];
    }

    Point projected_pt(Ambient_dimension<Point>::value, 
                       p_proj.begin(), p_proj.end());
    return Tr_point(
      Tr_bare_point(Intrinsic_dimension, coords.begin(), coords.end()), 
      m_k.squared_distance_d_object()(p, projected_pt));
  }

  // A simplex here is a list of point indices
  bool is_simplex_consistent(std::set<std::size_t> const& simplex)
  {
    // Check if the simplex is in the stars of all its vertices
    std::set<std::size_t>::const_iterator it_point_idx = simplex.begin();
    // For each point
    for ( ; it_point_idx != simplex.end() ; ++it_point_idx)
    {
      std::size_t point_idx = *it_point_idx;
      Triangulation const& tr = m_triangulations[point_idx].tr();
      Tr_vertex_handle center_vh = m_triangulations[point_idx].center_vertex();

      std::vector<Tr_full_cell_handle> incident_cells;
      tr.incident_full_cells(center_vh, std::back_inserter(incident_cells));

      std::vector<Tr_full_cell_handle>::const_iterator it_c = incident_cells.begin();
      std::vector<Tr_full_cell_handle>::const_iterator it_c_end= incident_cells.end();
      // For each cell
      bool found = false;
      for ( ; !found && it_c != it_c_end ; ++it_c)
      {
        std::set<std::size_t> cell;
        for (int i = 0 ; i < Intrinsic_dimension + 1 ; ++i)
          cell.insert((*it_c)->vertex(i)->data());
        if (cell == simplex)
          found = true;
      }

      if (!found)
        return false;
    }
    
    return true;
  }

private:
  const Kernel        m_k;
  Points              m_points;
  Points_ds           m_points_ds;
  TS_container        m_tangent_spaces;
  Tr_container        m_triangulations; // Contains the triangulations 
                                        // and their center vertex

}; // /class Tangential_complex

}  // end namespace CGAL

#endif // TANGENTIAL_COMPLEX_H
