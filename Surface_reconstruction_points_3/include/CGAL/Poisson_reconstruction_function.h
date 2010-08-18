// Copyright (c) 2007-09  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Laurent Saboret, Pierre Alliez

#ifndef CGAL_POISSON_RECONSTRUCTION_FUNCTION_H
#define CGAL_POISSON_RECONSTRUCTION_FUNCTION_H

#include <vector>
#include <deque>
#include <algorithm>
#include <cmath>

#include <CGAL/trace.h>
#include <CGAL/Reconstruction_triangulation_3.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Taucs_solver_traits.h>
#include <CGAL/centroid.h>
#include <CGAL/property_map.h>
#include <CGAL/surface_reconstruction_points_assertions.h>
#include <CGAL/Memory_sizer.h>
#include <CGAL/Peak_memory_sizer.h>
#include <CGAL/poisson_refine_triangulation.h>
#include <CGAL/Robust_circumcenter_filtered_traits_3.h>

#include <boost/shared_ptr.hpp>

namespace CGAL {


/// Given a set of 3D points with oriented normals sampled on the boundary of a 3D solid,
/// the Poisson Surface Reconstruction method [Kazhdan06] solves for an approximate indicator function
/// of the inferred solid, whose gradient best matches the input normals.
/// The output scalar function, represented in an adaptive octree, is then iso-contoured
/// using an adaptive marching cubes.
///
/// Poisson_reconstruction_function implements a variant of this algorithm which solves
/// for a piecewise linear function on a 3D Delaunay triangulation instead of an adaptive octree.
///
/// @heading Is Model for the Concepts:
/// Model of the 'ImplicitFunction' concept.
///
/// @heading Parameters:
/// @param Gt Geometric traits class.

template <class Gt>
class Poisson_reconstruction_function
{
// Public types
public:

  typedef Gt Geom_traits; ///< Geometric traits class

  // Geometric types
  typedef typename Geom_traits::FT FT; ///< typedef to Geom_traits::FT
  typedef typename Geom_traits::Point_3 Point; ///< typedef to Geom_traits::Point_3
  typedef typename Geom_traits::Vector_3 Vector; ///< typedef to Geom_traits::Vector_3
  typedef typename Geom_traits::Sphere_3 Sphere; ///< typedef to Geom_traits::Sphere_3

// Private types
private:

  /// Internal 3D triangulation, of type Reconstruction_triangulation_3.
  // Note: poisson_refine_triangulation() requires a robust circumcenter computation.
  typedef Reconstruction_triangulation_3<Robust_circumcenter_filtered_traits_3<Gt> >
                                                   Triangulation;

  // Repeat Triangulation types
  typedef typename Triangulation::Triangulation_data_structure Triangulation_data_structure;
  typedef typename Geom_traits::Ray_3 Ray;
  typedef typename Geom_traits::Plane_3 Plane;
  typedef typename Geom_traits::Segment_3 Segment;
  typedef typename Geom_traits::Triangle_3 Triangle;
  typedef typename Geom_traits::Tetrahedron_3 Tetrahedron;
  typedef typename Triangulation::Cell_handle   Cell_handle;
  typedef typename Triangulation::Vertex_handle Vertex_handle;
  typedef typename Triangulation::Cell   Cell;
  typedef typename Triangulation::Vertex Vertex;
  typedef typename Triangulation::Facet  Facet;
  typedef typename Triangulation::Edge   Edge;
  typedef typename Triangulation::Cell_circulator  Cell_circulator;
  typedef typename Triangulation::Facet_circulator Facet_circulator;
  typedef typename Triangulation::Cell_iterator    Cell_iterator;
  typedef typename Triangulation::Facet_iterator   Facet_iterator;
  typedef typename Triangulation::Edge_iterator    Edge_iterator;
  typedef typename Triangulation::Vertex_iterator  Vertex_iterator;
  typedef typename Triangulation::Point_iterator   Point_iterator;
  typedef typename Triangulation::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Triangulation::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename Triangulation::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename Triangulation::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename Triangulation::All_cells_iterator       All_cells_iterator;
  typedef typename Triangulation::Locate_type Locate_type;

// Data members.
// Warning: the Surface Mesh Generation package makes copies of implicit functions,
// thus this class must be lightweight and stateless.
private:

  // operator() is pre-computed on vertices of *m_tr by solving
  // the Poisson equation Laplacian(f) = divergent(normals field).
  boost::shared_ptr<Triangulation> m_tr;

  // contouring and meshing
  Point m_sink; // Point with the minimum value of operator()
  mutable Cell_handle m_hint; // last cell found = hint for next search

// Public methods
public:

  /// Creates a Poisson implicit function from the [first, beyond) range of points.
  ///
  /// @commentheading Template Parameters:
  /// @param InputIterator iterator over input points.
  /// @param PointPMap is a model of boost::ReadablePropertyMap with a value_type = Point_3.
  ///        It can be omitted if InputIterator value_type is convertible to Point_3.
  /// @param NormalPMap is a model of boost::ReadablePropertyMap with a value_type = Vector_3.

  // This variant requires all parameters.
  template <typename InputIterator,
            typename PointPMap,
            typename NormalPMap
  >
  Poisson_reconstruction_function(
    InputIterator first,  ///< iterator over the first input point.
    InputIterator beyond, ///< past-the-end iterator over the input points.
    PointPMap point_pmap, ///< property map to access the position of an input point.
    NormalPMap normal_pmap) ///< property map to access the *oriented* normal of an input point.
  : m_tr(new Triangulation)
  {
    CGAL::Timer task_timer; task_timer.start();
    CGAL_TRACE_STREAM << "Creates Poisson triangulation...\n";

    // Inserts points in triangulation
    m_tr->insert(
      first,beyond,
      point_pmap,
      normal_pmap);

    // Prints status
    CGAL_TRACE_STREAM << "Creates Poisson triangulation: " << task_timer.time() << " seconds, "
                                                           << (CGAL::Memory_sizer().virtual_size()>>20) << " Mb allocated"
                                                           << std::endl;
  }

  /// @cond SKIP_IN_MANUAL
  // This variant creates a default point property map = Dereference_property_map.
  template <typename InputIterator,
            typename NormalPMap
  >
  Poisson_reconstruction_function(
    InputIterator first,  ///< iterator over the first input point.
    InputIterator beyond, ///< past-the-end iterator over the input points.
    NormalPMap normal_pmap) ///< property map to access the *oriented* normal of an input point.
  : m_tr(new Triangulation)
  {
    CGAL::Timer task_timer; task_timer.start();
    CGAL_TRACE_STREAM << "Creates Poisson triangulation...\n";

    // Inserts points in triangulation
    m_tr->insert(
      first,beyond,
      normal_pmap);

    // Prints status
    CGAL_TRACE_STREAM << "Creates Poisson triangulation: " << task_timer.time() << " seconds, "
                                                           << (CGAL::Memory_sizer().virtual_size()>>20) << " Mb allocated"
                                                           << std::endl;
  }
  /// @endcond

  /// Returns a sphere bounding the inferred surface.
  Sphere bounding_sphere() const
  {
    return m_tr->input_points_bounding_sphere();
  }

  /// The function compute_implicit_function() must be called
  /// after the insertion of oriented points.
  /// It computes the piecewise linear scalar function operator() by:
  /// - applying Delaunay refinement,
  /// - solving for operator() at each vertex of the triangulation with a sparse linear solver,
  /// - and shifting and orienting operator() such that it is 0 at all input points and negative inside the inferred surface.
  ///
  /// @commentheading Template parameters:
  /// @param SparseLinearAlgebraTraits_d Symmetric definite positive sparse linear solver.
  /// The default solver is TAUCS Multifrontal Supernodal Cholesky Factorization.
  ///
  /// @return false if the linear solver fails.

  // This variant requires all parameters.
  template <class SparseLinearAlgebraTraits_d>
  bool compute_implicit_function(
    SparseLinearAlgebraTraits_d solver = SparseLinearAlgebraTraits_d()) ///< sparse linear solver
  {
    CGAL::Timer task_timer; task_timer.start();
    CGAL_TRACE_STREAM << "Delaunay refinement...\n";

    // Delaunay refinement
    const FT radius_edge_ratio_bound = 2.5;
    const unsigned int max_vertices = (unsigned int)1e7; // max 10M vertices
    const FT enlarge_ratio = 1.5;
    const FT radius = sqrt(bounding_sphere().squared_radius()); // get triangulation's radius
    const FT cell_radius_bound = radius/5.; // large
    unsigned int nb_vertices_added = delaunay_refinement(radius_edge_ratio_bound,cell_radius_bound,max_vertices,enlarge_ratio);

    // Prints status
    CGAL_TRACE_STREAM << "Delaunay refinement: " << "added " << nb_vertices_added << " Steiner points, "
                                                 << task_timer.time() << " seconds, "
                                                 << (CGAL::Memory_sizer().virtual_size()>>20) << " Mb allocated"
                                                 << std::endl;
    task_timer.reset();

    CGAL_TRACE_STREAM << "Solve Poisson equation with normalized divergence...\n";

    // Computes the Poisson indicator function operator()
    // at each vertex of the triangulation.
    double lambda = 0.1;
    if ( ! solve_poisson(solver, lambda) )
    {
      std::cerr << "Error: cannot solve Poisson equation" << std::endl;
      return false;
    }

    // Shift and orient operator() such that:
    // - operator() = 0 on the input points,
    // - operator() < 0 inside the surface.
    set_contouring_value(median_value_at_input_vertices());

    // Prints status
    CGAL_TRACE_STREAM << "Solve Poisson equation: " << task_timer.time() << " seconds, "
                                                    << (CGAL::Memory_sizer().virtual_size()>>20) << " Mb allocated"
                                                    << std::endl;
    task_timer.reset();

    return true;
  }

  /// @cond SKIP_IN_MANUAL
  // This variant provides the default sparse linear traits class = Taucs_symmetric_solver_traits.
  bool compute_implicit_function()
  {
    return compute_implicit_function< Taucs_symmetric_solver_traits<double> >();
  }
  /// @endcond

  /// 'ImplicitFunction' interface: evaluates the implicit function at a given 3D query point.
  FT operator()(const Point& p) const
  {
    m_hint = m_tr->locate(p,m_hint);

    if(m_hint == NULL)
      return 1e38;

    if(m_tr->is_infinite(m_hint))
      return 1e38;

    FT a,b,c,d;
    barycentric_coordinates(p,m_hint,a,b,c,d);
    return a * m_hint->vertex(0)->f() +
           b * m_hint->vertex(1)->f() +
           c * m_hint->vertex(2)->f() +
           d * m_hint->vertex(3)->f();
  }

  /// Returns a point located inside the inferred surface.
  Point get_inner_point() const
  {
    // Gets point / the implicit function is minimum
    return m_sink;
  }

// Private methods:
private:

  /// Delaunay refinement (break bad tetrahedra, where
  /// bad means badly shaped or too big). The normal of
  /// Steiner points is set to zero.
  /// Returns the number of vertices inserted.
  unsigned int delaunay_refinement(FT radius_edge_ratio_bound, ///< radius edge ratio bound (ignored if zero)
                                   FT cell_radius_bound, ///< cell radius bound (ignored if zero)
                                   unsigned int max_vertices, ///< number of vertices bound
                                   FT enlarge_ratio) ///< bounding box enlarge ratio
  {
    CGAL_TRACE("Calls delaunay_refinement(radius_edge_ratio_bound=%lf, cell_radius_bound=%lf, max_vertices=%u, enlarge_ratio=%lf)\n",
               radius_edge_ratio_bound, cell_radius_bound, max_vertices, enlarge_ratio);

    Sphere elarged_bsphere = enlarged_bounding_sphere(enlarge_ratio);
    unsigned int nb_vertices_added = poisson_refine_triangulation(*m_tr,radius_edge_ratio_bound,cell_radius_bound,max_vertices,elarged_bsphere);

    CGAL_TRACE("End of delaunay_refinement()\n");

    return nb_vertices_added;
  }

  /// Poisson reconstruction.
  /// Returns false on error.
  ///
  /// @commentheading Template parameters:
  /// @param SparseLinearAlgebraTraits_d Symmetric definite positive sparse linear solver.
  template <class SparseLinearAlgebraTraits_d>
  bool solve_poisson(
    SparseLinearAlgebraTraits_d solver, ///< sparse linear solver
    double lambda)
  {
    CGAL_TRACE("Calls solve_poisson()\n");

    double time_init = clock();

    double duration_assembly = 0.0;
    double duration_solve = 0.0;

    long old_max_memory = CGAL::Peak_memory_sizer().peak_virtual_size();

    CGAL_TRACE("  %ld Mb allocated\n", long(CGAL::Memory_sizer().virtual_size()>>20));
    CGAL_TRACE("  Creates matrix...\n");

    // get #variables
    unsigned int nb_variables = m_tr->index_unconstrained_vertices();

    // at least one vertex must be constrained
    if(nb_variables == m_tr->number_of_vertices())
    {
      constrain_one_vertex_on_convex_hull();
      nb_variables = m_tr->index_unconstrained_vertices();
    }

    // Assemble linear system A*X=B
    typename SparseLinearAlgebraTraits_d::Matrix A(nb_variables); // matrix is symmetric definite positive
    typename SparseLinearAlgebraTraits_d::Vector X(nb_variables), B(nb_variables);

    Finite_vertices_iterator v, e;
    for(v = m_tr->finite_vertices_begin(),
        e = m_tr->finite_vertices_end();
        v != e;
        ++v)
    {
      if(!v->constrained())
      {
        B[v->index()] = div_normalized(v); // rhs -> divergent
        assemble_poisson_row<SparseLinearAlgebraTraits_d>(A,v,B,lambda);
      }
    }
    
    duration_assembly = (clock() - time_init)/CLOCKS_PER_SEC;
    CGAL_TRACE("  Creates matrix: done (%.2lf s)\n", duration_assembly);

    CGAL_TRACE("  %ld Mb allocated\n", long(CGAL::Memory_sizer().virtual_size()>>20));
    CGAL_TRACE("  Solve sparse linear system...\n");

    // Solve "A*X = B". On success, solution is (1/D) * X.
    time_init = clock();
    double D;
    if(!solver.linear_solver(A, B, X, D))
      return false;
    CGAL_surface_reconstruction_points_assertion(D == 1.0);
    duration_solve = (clock() - time_init)/CLOCKS_PER_SEC;

    // Prints peak memory (Windows only)
    long max_memory = CGAL::Peak_memory_sizer().peak_virtual_size();
    if (max_memory <= 0) { // if peak_virtual_size() not implemented
        CGAL_TRACE("  Sorry. Cannot get solver max allocation on this system.\n");
    } else {
      if (max_memory > old_max_memory) {
        CGAL_TRACE("  Max allocation in solver = %ld Mb\n", max_memory>>20);
      } else {
        CGAL_TRACE("  Sorry. Failed to get solver max allocation.\n");
        CGAL_TRACE("  Max allocation since application start = %ld Mb\n", max_memory>>20);
      }
    }

    CGAL_TRACE("  Solve sparse linear system: done (%.2lf s)\n", duration_solve);

    // copy function's values to vertices
    unsigned int index = 0;
    for (v = m_tr->finite_vertices_begin(), e = m_tr->finite_vertices_end(); v!= e; ++v)
      if(!v->constrained())
        v->f() = X[index++];

    CGAL_TRACE("  %ld Mb allocated\n", long(CGAL::Memory_sizer().virtual_size()>>20));
    CGAL_TRACE("End of solve_poisson()\n");

    return true;
  }

  /// Shift and orient the implicit function such that:
  /// - the implicit function = 0 for points / f() = contouring_value,
  /// - the implicit function < 0 inside the surface.
  ///
  /// Returns the minimum value of the implicit function.
  FT set_contouring_value(FT contouring_value)
  {
    // median value set to 0.0
    shift_f(-contouring_value);

    // check value on convex hull (should be positive)
    Vertex_handle v = any_vertex_on_convex_hull();
    if(v->f() < 0.0)
      flip_f();

    // Update m_sink
    FT sink_value = find_sink();
    return sink_value;
  }


/// Gets median value of the implicit function over input vertices.
  FT median_value_at_input_vertices() const
  {
    std::deque<FT> values;
    Finite_vertices_iterator v, e;
    for(v = m_tr->finite_vertices_begin(),
        e= m_tr->finite_vertices_end();
        v != e; 
        v++)
      if(v->type() == Triangulation::INPUT)
        values.push_back(v->f());

    int size = values.size();
    if(size == 0)
    {
      std::cerr << "Contouring: no input points\n";
      return 0.0;
    }

    std::sort(values.begin(),values.end());
    int index = size/2;
    // return values[size/2];
    return 0.5 * (values[index] + values[index+1]); // avoids singular cases
  }

  void barycentric_coordinates(const Point& p,
                               Cell_handle cell,
                               FT& a,
                               FT& b,
                               FT& c,
                               FT& d) const
  {
    const Point& pa = cell->vertex(0)->point();
    const Point& pb = cell->vertex(1)->point();
    const Point& pc = cell->vertex(2)->point();
    const Point& pd = cell->vertex(3)->point();

    FT v = volume(pa,pb,pc,pd);
    a = std::fabs(volume(pb,pc,pd,p) / v);
    b = std::fabs(volume(pa,pc,pd,p) / v);
    c = std::fabs(volume(pb,pa,pd,p) / v);
    d = std::fabs(volume(pb,pc,pa,p) / v);
  }

  FT find_sink()
  {
    m_sink = CGAL::ORIGIN;
    FT min_f = 1e38;
    Finite_vertices_iterator v, e;
    for(v = m_tr->finite_vertices_begin(),
        e= m_tr->finite_vertices_end();
        v != e;
        v++)
    {
      if(v->f() < min_f)
      {
        m_sink = v->point();
        min_f = v->f();
      }
    }
    return min_f;
  }

  void shift_f(const FT shift)
  {
    Finite_vertices_iterator v, e;
    for(v = m_tr->finite_vertices_begin(),
        e = m_tr->finite_vertices_end();
        v!= e;
        v++)
      v->f() += shift;
  }

  void flip_f()
  {
    Finite_vertices_iterator v, e;
    for(v = m_tr->finite_vertices_begin(),
          e = m_tr->finite_vertices_end();
        v != e;
        v++)
      v->f() = -v->f();
  }

  Vertex_handle any_vertex_on_convex_hull()
  {
    // TODO: return NULL if none and assert
    std::vector<Vertex_handle> vertices;
    vertices.reserve(32);
    m_tr->incident_vertices(m_tr->infinite_vertex(),std::back_inserter(vertices));
    typename std::vector<Vertex_handle>::iterator it = vertices.begin();
    return *it;
  }

  void constrain_one_vertex_on_convex_hull(const FT value = 0.0)
  {
    Vertex_handle v = any_vertex_on_convex_hull();
    v->constrained() = true;
    v->f() = value;
  }

  // divergent
  FT div_normalized(Vertex_handle v)
  {
    std::vector<Cell_handle> cells;
    cells.reserve(32);
    m_tr->incident_cells(v,std::back_inserter(cells));
    if(cells.size() == 0)
      return 0.0;

    FT length = 100000;
    int counter = 0;
    FT div = 0.0;
    typename std::vector<Cell_handle>::iterator it;
    for(it = cells.begin(); it != cells.end(); it++)
    {
      Cell_handle cell = *it;
      if(m_tr->is_infinite(cell))
        continue;

      // compute average normal per cell
      Vector n = cell_normal(cell);

      // zero normal - no need to compute anything else
      if(n == CGAL::NULL_VECTOR)
        continue;

      // compute n'
      int index = cell->index(v);
      const Point& x = cell->vertex(index)->point();
      const Point& a = cell->vertex((index+1)%4)->point();
      const Point& b = cell->vertex((index+2)%4)->point();
      const Point& c = cell->vertex((index+3)%4)->point();
      Vector nn = (index%2==0) ? CGAL::cross_product(b-a,c-a) : CGAL::cross_product(c-a,b-a);
      nn = nn / std::sqrt(nn*nn); // normalize
      Vector p = a - x;
      Vector q = b - x;
      Vector r = c - x;
      FT p_n = std::sqrt(p*p);
      FT q_n = std::sqrt(q*q);
      FT r_n = std::sqrt(r*r);
      FT solid_angle = p*(CGAL::cross_product(q,r));
      solid_angle = std::abs(solid_angle * 1.0 / (p_n*q_n*r_n + (p*q)*r_n + (q*r)*p_n + (r*p)*q_n));
      Triangle face(a,b,c);
      FT area = std::sqrt(face.squared_area());
      length = std::sqrt((x-a)*(x-a)) + std::sqrt((x-b)*(x-b)) + std::sqrt((x-c)*(x-c));
      counter++;
      div += n * nn * area * 3 / length ;
    }
    return div;
  }

  Vector cell_normal(Cell_handle cell)
  {
    const Vector& n0 = cell->vertex(0)->normal();
    const Vector& n1 = cell->vertex(1)->normal();
    const Vector& n2 = cell->vertex(2)->normal();
    const Vector& n3 = cell->vertex(3)->normal();
    Vector n = n0 + n1 + n2 + n3;
    FT sq_norm = n*n;
    if(sq_norm != 0.0)
      return n / std::sqrt(sq_norm); // normalize
    else
      return CGAL::NULL_VECTOR;
  }

  // cotan formula as area(voronoi face) / len(primal edge)
  FT cotan_geometric(Edge& edge)
  {
    Cell_handle cell = edge.first;
    Vertex_handle vi = cell->vertex(edge.second);
    Vertex_handle vj = cell->vertex(edge.third);

    // primal edge
    const Point& pi = vi->point();
    const Point& pj = vj->point();
    Vector primal = pj - pi;
    FT len_primal = std::sqrt(primal * primal);
    return area_voronoi_face(edge) / len_primal;
  }

  // spin around edge
  // return area(voronoi face)
  FT area_voronoi_face(Edge& edge)
  {
    // circulate around edge
    Cell_circulator circ = m_tr->incident_cells(edge);
    Cell_circulator done = circ;
    std::vector<Point> voronoi_points;
    do
    {
      Cell_handle cell = circ;
      if(!m_tr->is_infinite(cell))
        voronoi_points.push_back(m_tr->dual(cell));
      else // one infinite tet, switch to another calculation
        return area_voronoi_face_boundary(edge);
      circ++;
    }
    while(circ != done);

    if(voronoi_points.size() < 3)
    {
      CGAL_surface_reconstruction_points_assertion(false);
      return 0.0;
    }

    // sum up areas
    FT area = 0.0;
    const Point& a = voronoi_points[0];
    unsigned int nb_triangles = voronoi_points.size() - 2;
    for(unsigned int i=1;i<nb_triangles;i++)
    {
      const Point& b = voronoi_points[i];
      const Point& c = voronoi_points[i+1];
      Triangle triangle(a,b,c);
      area += std::sqrt(triangle.squared_area());
    }
    return area;
  }

  // approximate area when a cell is infinite
  FT area_voronoi_face_boundary(Edge& edge)
  {
    FT area = 0.0;
    Vertex_handle vi = edge.first->vertex(edge.second);
    Vertex_handle vj = edge.first->vertex(edge.third);

    const Point& pi = vi->point();
    const Point& pj = vj->point();
    Point m = CGAL::midpoint(pi,pj);

    // circulate around each incident cell
    Cell_circulator circ = m_tr->incident_cells(edge);
    Cell_circulator done = circ;
    do
    {
      Cell_handle cell = circ;
      if(!m_tr->is_infinite(cell))
      {
        // circumcenter of cell
        Point c = m_tr->dual(cell);
        Tetrahedron tet = m_tr->tetrahedron(cell);

        int i = cell->index(vi);
        int j = cell->index(vj);
        int k = -1, l = -1;
        other_two_indices(i,j, &k,&l);
        Vertex_handle vk = cell->vertex(k);
        Vertex_handle vl = cell->vertex(l);

        const Point& pk = vk->point();
        const Point& pl = vl->point();

        // if circumcenter is outside tet
        // pick barycenter instead
        if(tet.has_on_unbounded_side(c))
        {
          Point cell_points[4] = {pi,pj,pk,pl};
          c = CGAL::centroid(cell_points, cell_points+4);
        }

        Point ck = CGAL::circumcenter(pi,pj,pk);
        Point cl = CGAL::circumcenter(pi,pj,pl);

        Triangle mcck(m,m,ck);
        Triangle mccl(m,m,cl);

        area += std::sqrt(mcck.squared_area());
        area += std::sqrt(mccl.squared_area());
      }
      circ++;
    }
    while(circ != done);
    return area;
  }

  // Gets indices different from i and j
  void other_two_indices(int i, int j, int* k, int* l)
  {
    CGAL_surface_reconstruction_points_assertion(i != j);
    bool k_done = false;
    bool l_done = false;
    for(int index=0;index<4;index++)
    {
      if(index != i && index != j)
      {
        if(!k_done)
        {
          *k = index;
          k_done = true;
        }
        else
        {
          *l = index;
          l_done = true;
        }
      }
    }
    CGAL_surface_reconstruction_points_assertion(k_done);
    CGAL_surface_reconstruction_points_assertion(l_done);
  }

  /// Assemble vi's row of the linear system A*X=B
  ///
  /// @commentheading Template parameters:
  /// @param SparseLinearAlgebraTraits_d Symmetric definite positive sparse linear solver.
  template <class SparseLinearAlgebraTraits_d>
  void assemble_poisson_row(typename SparseLinearAlgebraTraits_d::Matrix& A,
                            Vertex_handle vi,
                            typename SparseLinearAlgebraTraits_d::Vector& B,
                            double lambda)
  {
    // for each vertex vj neighbor of vi
    std::vector<Edge> edges;
    m_tr->incident_edges(vi,std::back_inserter(edges));

    double diagonal = 0.0;

  for(typename std::vector<Edge>::iterator it = edges.begin();
        it != edges.end();
        it++)
    {
      Vertex_handle vj = it->first->vertex(it->third);
      if(vj == vi){
        vj = it->first->vertex(it->second);
      }
      if(m_tr->is_infinite(vj))
        continue;

      // get corresponding edge
      Edge edge( it->first, it->first->index(vi), it->first->index(vj));
      if(vi->index() < vj->index()){
        std::swap(edge.second,  edge.third);
      }

      double cij = cotan_geometric(edge);
      if(vj->constrained())
        B[vi->index()] -= cij * vj->f(); // change rhs
      else
        A.set_coef(vi->index(),vj->index(), -cij, true /*new*/); // off-diagonal coefficient

      diagonal += cij;
    }
    // diagonal coefficient
    if (vi->type() == Triangulation::INPUT)
      A.set_coef(vi->index(),vi->index(), diagonal + lambda, true /*new*/) ;
    else
      A.set_coef(vi->index(),vi->index(), diagonal, true /*new*/);

  }


  /// Computes enlarged geometric bounding sphere of the embedded triangulation.
  Sphere enlarged_bounding_sphere(FT ratio) const
  {
    Sphere bsphere = bounding_sphere(); // triangulation's bounding sphere
    return Sphere(bsphere.center(), bsphere.squared_radius() * ratio*ratio);
  }

}; // end of Poisson_reconstruction_function


} //namespace CGAL

#endif // CGAL_POISSON_RECONSTRUCTION_FUNCTION_H
