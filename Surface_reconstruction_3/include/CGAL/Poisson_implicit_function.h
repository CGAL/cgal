// Copyright (c) 2007  INRIA (France).
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
//
// Author(s)     : Laurent Saboret, Pierre Alliez


#ifndef CGAL_POISSON_IMPLICIT_FUNCTION_H
#define CGAL_POISSON_IMPLICIT_FUNCTION_H

#include <queue>
#include <list>
#include <algorithm>

#include <CGAL/Implicit_fct_delaunay_triangulation_3.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/taucs_solver.h>
#include <CGAL/k_nearest_neighbor.h>
#include <CGAL/surface_reconstruction_assertions.h>
#include <CGAL/estimate_normals_pca_3.h>
#include <CGAL/estimate_normals_jet_fitting_3.h>

CGAL_BEGIN_NAMESPACE


// functor for priority queue
template<class Candidate>
struct less // read more priority
{
  bool operator()(const Candidate& c1,
                  const Candidate& c2) const
  {
    return (c1.score() < c2.score());
  }
};

// functor for priority queue
template<class Candidate>
struct more // read more priority
{
  bool operator()(const Candidate& c1,
                  const Candidate& c2) const
  {
    return (c1.score() > c2.score());
  }
};

template <class Handle, class Point>
class Candidate
{
private:
  Handle m_v0;
  Handle m_v1;
  Handle m_v2;
  Handle m_v3;
  float m_score;

public:

	Candidate(Handle v0,
		        Handle v1,
					  Handle v2,
					  Handle v3,
            const float score)
  {
		m_v0 = v0;
		m_v1 = v1;
		m_v2 = v2;
		m_v3 = v3;
		m_score = score;
  }
  ~Candidate() {}

public:

	const float& score() const { return m_score; }
  float& score() { return m_score; }

	Handle v0() { return m_v0; }
	Handle v1() { return m_v1; }
	Handle v2() { return m_v2; }
	Handle v3() { return m_v3; }

};

/// Poisson_implicit_function computes an indicator function f() piecewise-linear
/// over the tetrahedra. We solve the Poisson equation
/// Laplacian(f) = divergent(normals field) at each vertex
/// of the triangulation via the TAUCS sparse linear
/// solver. One vertex must be constrained.
///
/// @heading Is Model for the Concepts: Model of the Reconstruction_implicit_function concept.
///
/// @heading Design Pattern:
/// Poisson_implicit_function is a
/// Strategy [GHJV95]: it implements a strategy of surface mesh reconstruction.
///
/// @heading Parameters:
/// @param ImplicitFctDelaunayTriangulation_3 3D Delaunay triangulation, 
///        model of ImplicitFctDelaunayTriangulation_3 concept.

template <class Gt, class ImplicitFctDelaunayTriangulation_3>
class Poisson_implicit_function
{
// Public types
public:

  typedef ImplicitFctDelaunayTriangulation_3 Triangulation;

  typedef Gt Geom_traits; ///< Kernel's geometric traits
  typedef typename Geom_traits::FT FT;
	typedef typename Geom_traits::Point_3 Point;
	typedef typename Geom_traits::Iso_cuboid_3 Iso_cuboid_3;
	typedef typename Geom_traits::Sphere_3 Sphere;

  typedef typename Triangulation::Point_with_normal Point_with_normal; ///< Model of PointWithNormal_3
  typedef typename Triangulation::Normal Normal; ///< Model of OrientedNormal_3 concept.

// Private types
private:

  // Repeat ImplicitFctDelaunayTriangulation_3 types
  typedef typename Triangulation::Triangulation_data_structure Triangulation_data_structure;
	typedef typename Geom_traits::Ray_3 Ray;
	typedef typename Geom_traits::Plane_3 Plane;
	typedef typename Geom_traits::Vector_3 Vector;
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

// Data members
private:

	Triangulation& m_dt; // f() is pre-computed on vertices of m_dt by solving
                       // the Poisson equation Laplacian(f) = divergent(normals field)

	// neighbor search
	typedef typename CGAL::K_nearest_neighbor<Geom_traits,Vertex_handle> K_nearest_neighbor;
	typedef typename CGAL::Point_vertex_handle_3<Vertex_handle> Point_vertex_handle_3;
	K_nearest_neighbor m_nn_search;

	// delaunay refinement
	typedef typename CGAL::Candidate<Vertex_handle,Point> Candidate;
	typedef typename std::priority_queue< Candidate,
                                        std::vector<Candidate>,
                                        less<Candidate> > Refinement_pqueue;

	// contouring and meshing
	Point m_sink; // Point with the minimum value of f()
	Cell_handle m_hint; // last cell found = hint for next search

// Public methods
public:

  /// Create a Poisson indicator function f() piecewise-linear
  /// over the tetrahedra of pdt.
  ///
  /// @param pdt ImplicitFctDelaunayTriangulation_3 base of the Poisson indicator function.
	Poisson_implicit_function(Triangulation& pdt)
	: m_dt(pdt)
	{
	}

 // /// Insert the first...beyond point set into pdt and
 // /// create a Poisson indicator function f() piecewise-linear
 // /// over the tetrahedra of pdt.
 // ///
 // /// Precondition: the value type of InputIterator must be 'Point' or Point_with_normal.
 // ///
 // /// @param pdt ImplicitFctDelaunayTriangulation_3 base of the Poisson indicator function.
 // /// @param first First point to add.
 // /// @param beyond Past-the-end point to add.
 // template < class InputIterator >
	//Poisson_implicit_function(Triangulation& pdt,
 //                           InputIterator first, InputIterator beyond)
	//: m_dt(pdt)
	//{
 //   insert(first, beyond);
	//}

 // /// Insert points.
 // ///
 // /// Precondition: the value type of InputIterator must be 'Point' or Point_with_normal.
 // ///
 // /// @param first First point to add to pdt.
 // /// @param beyond Past-the-end point to add to pdt.
 // /// @return the number of inserted points.
 // template < class InputIterator >
 // int insert(InputIterator first, InputIterator beyond)
	//{
 //   return m_dt.insert(first, beyond);
	//}

  /// Get embedded triangulation.
  Triangulation& triangulation()
  {
    return m_dt;
  }
  const Triangulation& triangulation() const
  {
    return m_dt;
  }

  /// Get the bounding box.
	Iso_cuboid_3 bounding_box() const
	{
		return m_dt.bounding_box();
	}

  /// Get bounding sphere.
	Sphere bounding_sphere() const
	{
		return m_dt.bounding_sphere();
	}

  /// Get the region of interest, ignoring the outliers.
  /// This method is used to define the OpenGL arcball sphere.
	Sphere region_of_interest() const
	{
    // A good candidate is a sphere containing the dense region of the point cloud:
    // - center point is barycenter
    // - Radius is 2 * standard deviation
    Point barycenter = m_dt.barycenter();
		float radius = 2.f * (float)m_dt.diameter_standard_deviation();

    return Sphere(barycenter, radius*radius);
	}

  /// You should call compute_implicit_function() once when points insertion is over.
  /// It computes the Poisson indicator function f()
  /// at each vertex of the triangulation by:
  /// - applying a Delaunay refinement to define the function
  ///   inside and outside the surface.
  /// - solving the Poisson equation
  ///   Laplacian(f) = divergent(normals field) at each vertex
  ///   of the triangulation via the TAUCS sparse linear
  ///   solver. One vertex must be constrained.
  /// - shifting and orienting f() such that f() = 0 on the input points,
  ///   and f() < 0 inside the surface.
  ///
  /// Return false on error.
  /// TODO: add parameters to compute_implicit_function()?
	bool compute_implicit_function()
	{
		// Delaunay refinement
	  const FT quality = 2.5;
	  const unsigned int max_vertices = (unsigned int)1e7; // max 10M vertices
	  const FT enlarge_ratio = 1.5;
	  delaunay_refinement(quality,max_vertices,enlarge_ratio,50000);

	  // Smooth normals field.
	  // Commented out as it shrinks the reconstructed model.
	  //extrapolate_normals();

	  // Solve Poisson equation
	  double duration_assembly, duration_factorization, duration_solve;
	  bool success = solve_poisson(&duration_assembly, &duration_factorization, &duration_solve);

	  // Shift and orient f() such that:
    // - f() = 0 on the input points,
    // - f() < 0 inside the surface.
    set_contouring_value(median_value_at_input_vertices());

    return success;
	}

  /// Estimate normal directions using linear least
  /// squares fitting of a plane on the k nearest neighbors.
	void estimate_normals_pca(unsigned int k)
	{
	  CGAL::estimate_normals_pca_3(m_dt.points_begin(), m_dt.points_end(), m_dt.normals_begin(), k);
	}

  /// Estimate normal directions using jet fitting on the k nearest
  /// neighbors.
	void estimate_normals_jet_fitting(unsigned int k)
	{
	  CGAL::estimate_normals_jet_fitting_3(m_dt.points_begin(), m_dt.points_end(), m_dt.normals_begin(), k);
	}

  /// Delaunay refinement (break bad tetrahedra, where
  /// bad means badly shaped or too big). The normal of
  /// Steiner points is set to zero.
  /// Return the number of vertices inserted.
	unsigned int delaunay_refinement(const FT threshold,
		                               const unsigned int maximum,
																	 const FT enlarge_ratio,
																	 const unsigned int restart_each)
	{

		// create enlarged bounding box
		Iso_cuboid_3 enlarged_bbox = enlarged_bounding_box(enlarge_ratio);

		// push all cells to the queue
		Refinement_pqueue queue;

		// init queue
		init_queue(queue,threshold,enlarged_bbox);

		unsigned int nb = 0;
		while(!queue.empty())
		{
			if(nb%restart_each == 0)
			{
				reset_queue(queue,threshold,enlarged_bbox);
			  if(queue.empty())
					break;
			}

			Candidate candidate = queue.top();
			queue.pop();
			Vertex_handle v0 = candidate.v0();
			Vertex_handle v1 = candidate.v1();
			Vertex_handle v2 = candidate.v2();
			Vertex_handle v3 = candidate.v3();
			Cell_handle cell = NULL;
			if(m_dt.is_cell(v0,v1,v2,v3,cell))
			{
				Point point = m_dt.dual(cell);
				Vertex_handle v = m_dt.insert(point,Triangulation::STEINER,cell);

				if(nb++ > maximum)
					break; // premature ending

				// iterate over incident cells and feed queue
				std::list<Cell_handle> cells;
				m_dt.incident_cells(v,std::back_inserter(cells));
				typename std::list<Cell_handle>::iterator it;
				for(it = cells.begin();
						it != cells.end();
						it++)
				{
					Cell_handle c = *it;
					if(m_dt.is_infinite(c))
						continue;

					FT rer = radius_edge_ratio(c);
					Point point = m_dt.dual(c);
					bool inside = enlarged_bbox.has_on_bounded_side(point);
					if(inside && rer > threshold)
					{
						Vertex_handle v0 = c->vertex(0);
						Vertex_handle v1 = c->vertex(1);
						Vertex_handle v2 = c->vertex(2);
						Vertex_handle v3 = c->vertex(3);
						float score = (float)max_edge_len(cell);
						queue.push(Candidate(v0,v1,v2,v3,score));
					}
				}
			}
		}
	  m_dt.invalidate_bounding_box();
		return nb;
	}

	unsigned int delaunay_refinement_shell(FT size_shell,
		                                     FT sizing,
		                                     const unsigned int maximum)
	{
		// make parameters relative to size
    Sphere bounding_sphere = m_dt.bounding_sphere();
		FT size = sqrt(bounding_sphere.squared_radius());
		size_shell *= size;
		sizing *= size;

		init_nn_search_shell();

		typedef typename CGAL::Candidate<Vertex_handle,Point> Candidate;
		typedef typename std::priority_queue<Candidate,
		                                     std::vector<Candidate>,
		                                     more<Candidate> > PQueue;

		// push all cells to the queue
		PQueue queue;
		Finite_cells_iterator c;
		for(c = m_dt.finite_cells_begin();
			  c != m_dt.finite_cells_end();
				c++)
		{
			Point p;
			FT size = 0.0;
			if(is_refinable(c,size_shell,sizing,size,p))
			{
				Vertex_handle v0 = c->vertex(0);
				Vertex_handle v1 = c->vertex(1);
				Vertex_handle v2 = c->vertex(2);
				Vertex_handle v3 = c->vertex(3);
				queue.push(Candidate(v0,v1,v2,v3,(float)size));
			}
		}

		unsigned int nb = 0;
		while(!queue.empty())
		{
			Candidate candidate = queue.top();
			queue.pop();
			Vertex_handle v0 = candidate.v0();
			Vertex_handle v1 = candidate.v1();
			Vertex_handle v2 = candidate.v2();
			Vertex_handle v3 = candidate.v3();

			Cell_handle cell = NULL;
			if(m_dt.is_cell(v0,v1,v2,v3,cell))
			{
				Point point = m_dt.dual(cell);
				Vertex_handle v = m_dt.insert(point, Triangulation::STEINER);

				if(nb++ > maximum)
					return nb; // premature ending

				// iterate over incident cells and feed queue
				std::list<Cell_handle> cells;
				m_dt.incident_cells(v,std::back_inserter(cells));
				typename std::list<Cell_handle>::iterator it;
				for(it = cells.begin();
						it != cells.end();
						it++)
				{
					Cell_handle c = *it;
					if(m_dt.is_infinite(c))
						continue;

					Point p;
					FT size = 0.0;
					if(is_refinable(c,size_shell,sizing,size,p))
					{
						Vertex_handle v0 = c->vertex(0);
						Vertex_handle v1 = c->vertex(1);
						Vertex_handle v2 = c->vertex(2);
						Vertex_handle v3 = c->vertex(3);
						queue.push(Candidate(v0,v1,v2,v3,(float)size));
					}
				}
			}
		}
		m_nn_search = K_nearest_neighbor();
		return nb;
	}

	/// Extrapolate the normals field:
	/// compute null normals by averaging neighbour normals.
	void extrapolate_normals()
	{
	  // Compute extrapolated normals and store them in extrapolated_normals[]
		std::map<Vertex_handle,Normal> extrapolated_normals; // vector + orientation
		Finite_vertices_iterator v;
		for(v = m_dt.finite_vertices_begin();
				v != m_dt.finite_vertices_end();
				v++)
		{
			if(v->normal().get_vector() != CGAL::NULL_VECTOR)
				continue;

			Vector normal = CGAL::NULL_VECTOR;  // normal vector to compute
      bool oriented_normal = true;        // normal orientation to compute

			std::list<Vertex_handle> vertices;
			m_dt.incident_vertices(v,std::back_inserter(vertices));
			for(typename std::list<Vertex_handle>::iterator it = vertices.begin();
					it != vertices.end();
					it++)
			{
				Vertex_handle nv = *it;
				normal = normal + nv->normal().get_vector();
				oriented_normal &=  nv->normal().is_normal_oriented();
			}

			FT sq_norm = normal * normal;
			if(sq_norm > 0.0)
        normal = normal / std::sqrt(sq_norm);

      extrapolated_normals[v] = Normal(normal, oriented_normal);
		}

		// set normals
		for(v = m_dt.finite_vertices_begin();
				v != m_dt.finite_vertices_end();
				v++)
		{
			if(v->normal().get_vector() != CGAL::NULL_VECTOR)
				continue;

			typename std::map<Vertex_handle,Normal>::iterator it = extrapolated_normals.find(v);
			if(it != extrapolated_normals.end())
				v->normal() = extrapolated_normals[v];
		}
	}

	/// Poisson reconstruction.
	/// Return false on error.
	bool solve_poisson(double* duration_assembly,
		                 double* duration_factorization,
			  						 double* duration_solve)
	{
		double time_init = clock();

	  *duration_assembly = 0.0;
	  *duration_factorization = 0.0;
	  *duration_solve = 0.0;

		// get #variables
		unsigned int nb_variables = m_dt.index_unconstrained_vertices();

		// at least one vertex must be constrained
		if(nb_variables == m_dt.number_of_vertices())
		{
			constrain_one_vertex_on_convex_hull();
			nb_variables = m_dt.index_unconstrained_vertices();
		}

		// Assemble linear system
		Taucs_solver solver;
		std::vector<double> X(nb_variables);
		std::vector<double> B(nb_variables);

		Finite_vertices_iterator v;
		for(v = m_dt.finite_vertices_begin();
			  v != m_dt.finite_vertices_end();
			  v++)
		{
			if(!v->constrained())
			{
				B[v->index()] = div(v); // rhs -> divergent
				assemble_poisson_row(solver,v,B);
			}
		}

		*duration_assembly = (clock() - time_init)/CLOCKS_PER_SEC;

		/*
		time_init = clock();
		if(!solver.solve_conjugate_gradient(B,X,10000,1e-15))
			return false;
		*duration_solve = (clock() - time_init)/CLOCKS_PER_SEC;
		*/

		// Choleschy factorization M = L L^T
		time_init = clock();
		if(!solver.factorize_ooc())
			return false;
		*duration_factorization = (clock() - time_init)/CLOCKS_PER_SEC;

		// direct solve by forward and backward substitution
		time_init = clock();
		if(!solver.solve_ooc(B,X))
			return false;
		*duration_solve = (clock() - time_init)/CLOCKS_PER_SEC;

		/*
		// Choleschy factorization M = L L^T
		time_init = clock();
		if(!solver.factorize(true))
			return false;
		*duration_factorization = (clock() - time_init)/CLOCKS_PER_SEC;

		// direct solve by forward and backward substitution
		time_init = clock();
		if(!solver.solve(B,X,1))
			return false;
		*duration_solve = (clock() - time_init)/CLOCKS_PER_SEC;
		*/

		// set values to vertices
		unsigned int index = 0;
		for(v = m_dt.finite_vertices_begin();
			  v != m_dt.finite_vertices_end();
			  v++)
			if(!v->constrained())
				v->f() = X[index++];

		return true;
	}

  /// Shift and orient the implicit function such that:
  /// - the implicit function = 0 for points / f() = contouring_value,
  /// - the implicit function < 0 inside the surface.
  ///
  /// Return the minimum value of the implicit function.
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

  /// Evaluate implicit function for any 3D point.
	FT f(const Point& p)
	{
		m_hint = m_dt.locate(p,m_hint);

		if(m_hint == NULL)
			return 1e38;

		if(m_dt.is_infinite(m_hint))
			return 1e38;

		FT a,b,c,d;
		barycentric_coordinates(p,m_hint,a,b,c,d);
		return a * m_hint->vertex(0)->f() +
			     b * m_hint->vertex(1)->f() +
				   c * m_hint->vertex(2)->f() +
					 d * m_hint->vertex(3)->f();
	}

  /// [ImplicitFunction interface]
  ///
  /// Evaluate implicit function for any 3D point.
	FT operator() (Point p)
	{
		return f(p);
	}

  /// Get point / the implicit function is minimum.
	const Point& sink() const { return m_sink; }

  /// Get average value of the implicit function over input vertices.
	FT average_value_at_input_vertices() const
	{
		FT sum = 0.0;
		unsigned int nb = 0;
		Finite_vertices_iterator v;
		for(v = m_dt.finite_vertices_begin();
				v != m_dt.finite_vertices_end();
				v++)
		{
			if(v->type() == Triangulation::INPUT)
			{
				sum += v->f();
				nb++;
			}
		}
		if(nb > 0)
			return sum / (FT)nb;
		else
		{
			std::cerr << "Contouring: no input points\n";
			return (FT)0.0;
		}
	}

  /// Get median value of the implicit function over input vertices.
	FT median_value_at_input_vertices() const
	{
		std::vector<FT> values;
		Finite_vertices_iterator v;
		for(v = m_dt.finite_vertices_begin();
				v != m_dt.finite_vertices_end();
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

  /// Get min value of the implicit function over input vertices.
	FT min_value_at_input_vertices() const
	{
		FT min_value = 1e38;
		Finite_vertices_iterator v;
		for(v = m_dt.finite_vertices_begin();
				v != m_dt.finite_vertices_end();
				v++)
		{
			if(v->type() == Triangulation::INPUT)
				min_value = (std::min)(min_value, v->f());
		}
		if (m_dt.number_of_vertices() > 0)
		{
			return min_value;
		}
		else
		{
			std::cerr << "Contouring: no input points\n";
			return (FT)0.0;
		}
	}

  /// Get max value of the implicit function over input vertices.
	FT max_value_at_input_vertices() const
	{
		FT max_value = -1e38;
		Finite_vertices_iterator v;
		for(v = m_dt.finite_vertices_begin();
				v != m_dt.finite_vertices_end();
				v++)
		{
			if(v->type() == Triangulation::INPUT)
				max_value = (std::max)(max_value, v->f());
		}
		if (m_dt.number_of_vertices() > 0)
		{
			return max_value;
		}
		else
		{
			std::cerr << "Contouring: no input points\n";
			return (FT)0.0;
		}
	}

  /// Get median value of the implicit function over convex hull vertices.
	FT median_value_at_convex_hull() const
	{
	  // Get convex hull vertices
		std::list<Vertex_handle> convex_hull_vertices;
		m_dt.incident_vertices(m_dt.infinite_vertex(),std::back_inserter(convex_hull_vertices));

		// Get values of the implicit function over convex hull vertices
		std::vector<FT> values;
		typename std::list<Vertex_handle>::iterator it;
		for(it = convex_hull_vertices.begin();
			  it != convex_hull_vertices.end();
				it++)
		{
			Vertex_handle v = *it;
  		values.push_back(v->f());
  	}

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

  /// Get average value of the implicit function over convex hull vertices.
	FT average_value_at_convex_hull() const
	{
		std::list<Vertex_handle> convex_hull_vertices;
		m_dt.incident_vertices(m_dt.infinite_vertex(),std::back_inserter(convex_hull_vertices));

		FT sum = 0.0;
		unsigned int nb = 0;
		typename std::list<Vertex_handle>::iterator it;
		for(it = convex_hull_vertices.begin();
			  it != convex_hull_vertices.end();
				it++,nb++)
		{
			Vertex_handle v = *it;
			sum += v->f();
		}
		if(nb != 0)
			return sum / (FT)nb;
		else
			return 0.0;
	}

// Private methods:
private:

	// PA: todo change type (FT)
	// check if this is in CGAL already
	void barycentric_coordinates(const Point& p,
		                           Cell_handle cell,
															 double& a,
															 double& b,
															 double& c,
															 double& d)
	{
		const Point& pa = cell->vertex(0)->point();
		const Point& pb = cell->vertex(1)->point();
		const Point& pc = cell->vertex(2)->point();
		const Point& pd = cell->vertex(3)->point();
		Tetrahedron ta(pb,pc,pd,p);
		Tetrahedron tb(pa,pc,pd,p);
		Tetrahedron tc(pb,pa,pd,p);
		Tetrahedron td(pb,pc,pa,p);
		Tetrahedron tet(pa,pb,pc,pd);
		double v = tet.volume();
		a = std::fabs(ta.volume() / v);
		b = std::fabs(tb.volume() / v);
		c = std::fabs(tc.volume() / v);
		d = std::fabs(td.volume() / v);
	}

	// radius-edge ratio (the ratio of the circumradius to
	// the shortest edge length of tetrahedron)
	// check template type
	double radius_edge_ratio(Cell_handle c)
	{
		double r = circumradius(c);
		double l = len_shortest_edge(c);
		if(l != 0.0)
			return r/l;
		else
			return 1e38;
	}

	FT len_shortest_edge(Cell_handle c)
	{
		FT d1 = distance(c->vertex(0),c->vertex(1));
		FT d2 = distance(c->vertex(0),c->vertex(2));
		FT d3 = distance(c->vertex(0),c->vertex(3));
		FT d4 = distance(c->vertex(1),c->vertex(2));
		FT d5 = distance(c->vertex(1),c->vertex(3));
		FT d6 = distance(c->vertex(2),c->vertex(3));
		return (std::min)((std::min)((std::min)(d1,d2),d3),
			                (std::min)((std::min)(d4,d5),d6));
	}

	FT distance(Vertex_handle v1,
		          Vertex_handle v2)
	{
		const Point& a = v1->point();
		const Point& b = v2->point();
		return std::sqrt(CGAL::squared_distance(a,b));
	}

	FT find_sink()
	{
		m_sink = CGAL::ORIGIN;
		FT min_f = 1e38;
		Finite_vertices_iterator v;
		for(v = m_dt.finite_vertices_begin();
			  v != m_dt.finite_vertices_end();
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
		Finite_vertices_iterator v;
		for(v = m_dt.finite_vertices_begin();
			  v != m_dt.finite_vertices_end();
			  v++)
		  v->f() += shift;
	}

	void flip_f()
	{
		Finite_vertices_iterator v;
		for(v = m_dt.finite_vertices_begin();
			  v != m_dt.finite_vertices_end();
			  v++)
		  v->f() = -v->f();
	}

	Vertex_handle any_vertex_on_convex_hull()
	{
		// TODO: return NULL if none and assert
		std::list<Vertex_handle> vertices;
		m_dt.incident_vertices(m_dt.infinite_vertex(),std::back_inserter(vertices));
		typename std::list<Vertex_handle>::iterator it = vertices.begin();
		return *it;
	}

	void constrain_one_vertex_on_convex_hull(const FT value = 0.0)
	{
		Vertex_handle v = any_vertex_on_convex_hull();
		v->constrained() = true;
		v->f() = value;
	}

	// divergent
	FT div(Vertex_handle v)
	{
	  std::list<Cell_handle> cells;
    m_dt.incident_cells(v,std::back_inserter(cells));
		if(cells.size() == 0)
		  return 0.0;

		FT div = 0.0;
		typename std::list<Cell_handle>::iterator it;
		for(it = cells.begin();
				it != cells.end();
				it++)
		{
			Cell_handle cell = *it;
			if(m_dt.is_infinite(cell))
				continue;

			// compute average normal per cell
			Vector n = cell_normal(cell);

			// zero normal - no need to compute anything else
			if(n == CGAL::NULL_VECTOR)
				continue;

			// compute n'
			int index = cell->index(v);
			const Point& a = cell->vertex((index+1)%4)->point();
			const Point& b = cell->vertex((index+2)%4)->point();
			const Point& c = cell->vertex((index+3)%4)->point();
			Vector nn = (index%2==0) ? CGAL::cross_product(b-a,c-a) : CGAL::cross_product(c-a,b-a);
			nn = nn / std::sqrt(nn*nn); // normalize

			Triangle face(a,b,c);
			FT area = std::sqrt(face.squared_area());

			div += n * nn * area;
		}
		return div;
	}

	Vector cell_normal(Cell_handle cell)
	{
		const Vector& n0 = cell->vertex(0)->normal().get_vector();
		const Vector& n1 = cell->vertex(1)->normal().get_vector();
		const Vector& n2 = cell->vertex(2)->normal().get_vector();
		const Vector& n3 = cell->vertex(3)->normal().get_vector();
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
		Cell_circulator circ = m_dt.incident_cells(edge);
		Cell_circulator done = circ;
		std::vector<Point> voronoi_points;
		do
		{
			Cell_handle cell = circ;
			if(!m_dt.is_infinite(cell))
				voronoi_points.push_back(m_dt.dual(cell));
			else // one infinite tet, switch to another calculation
				return area_voronoi_face_boundary(edge);
			circ++;
		}
		while(circ != done);

		if(voronoi_points.size() < 3)
		{
			CGAL_surface_reconstruction_assertion(false);
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
		Cell_circulator circ = m_dt.incident_cells(edge);
		Cell_circulator done = circ;
		do
		{
			Cell_handle cell = circ;
			if(!m_dt.is_infinite(cell))
			{
				// circumcenter of cell
				Point c = m_dt.dual(cell);
				Tetrahedron tet = m_dt.tetrahedron(cell);

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

  // Get indices different from i and j
	void other_two_indices(int i, int j, int* k, int* l)
	{
		CGAL_surface_reconstruction_assertion(i != j);
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
		CGAL_surface_reconstruction_assertion(k_done);
		CGAL_surface_reconstruction_assertion(l_done);
	}

	void assemble_poisson_row(Taucs_solver& solver,
		                        Vertex_handle vi,
														std::vector<double>& B)
	{
		// assemble new row
		solver.begin_row();

		std::list<Vertex_handle> vertices;
		m_dt.incident_vertices(vi,std::back_inserter(vertices));
		double diagonal = 0.0;
		for(typename std::list<Vertex_handle>::iterator it = vertices.begin();
				it != vertices.end();
				it++)
		{
			Vertex_handle vj = *it;
			if(m_dt.is_infinite(vj))
				continue;

			// get corresponding edge
			Edge edge = sorted_edge(vi,vj);

			double cij = cotan_geometric(edge);
			if(vj->constrained())
				B[vi->index()] -= cij * vj->f(); // change rhs
			else
				// off-diagonal coefficient
				solver.add_value(vj->index(),-cij);

			diagonal += cij;
		}

		// diagonal coefficient
		solver.add_value(vi->index(),diagonal);

		// end matrix row
		solver.end_row();
	}

	Edge sorted_edge(Vertex_handle vi,
		               Vertex_handle vj)
	{
		int i1 = 0;
		int i2 = 0;
		Cell_handle cell = NULL;
		if(vi->index() > vj->index())
		{
			bool success = m_dt.is_edge(vi,vj,cell,i1,i2);
			CGAL_surface_reconstruction_assertion(success);
			CGAL_surface_reconstruction_assertion(cell->vertex(i1) == vi);
			CGAL_surface_reconstruction_assertion(cell->vertex(i2) == vj);
		}
		else
		{
			bool success = m_dt.is_edge(vj,vi,cell,i1,i2);
			CGAL_surface_reconstruction_assertion(success);
			CGAL_surface_reconstruction_assertion(cell->vertex(i1) == vj);
			CGAL_surface_reconstruction_assertion(cell->vertex(i2) == vi);
		}
		return Edge(cell,i1,i2);
	}

  /// Compute enlarged geometric bounding box of the embedded triangulation
	Iso_cuboid_3 enlarged_bounding_box(FT ratio) const
	{
    // Get triangulation's bounding box
    Iso_cuboid_3 bbox = m_dt.bounding_box();

    // Its center point is:
		FT mx = 0.5 * (bbox.xmax() + bbox.xmin());
		FT my = 0.5 * (bbox.ymax() + bbox.ymin());
		FT mz = 0.5 * (bbox.zmax() + bbox.zmin());
		Point c(mx,my,mz);

    // Compute enlarged bounding box
		FT sx = 0.5 * ratio * (bbox.xmax() - bbox.xmin());
		FT sy = 0.5 * ratio * (bbox.ymax() - bbox.ymin());
		FT sz = 0.5 * ratio * (bbox.zmax() - bbox.zmin());
		Point p(c.x() - sx, c.y() - sy, c.z() - sz);
		Point q(c.x() + sx, c.y() + sy, c.z() + sz);
		return Iso_cuboid_3(p,q);
	}

	void reset_queue(Refinement_pqueue& queue,
		               const FT threshold,
									 Iso_cuboid_3& enlarged_bbox)
	{
		// clear queue
		while(!queue.empty())
			queue.pop();

		// fill up
		init_queue(queue,threshold,enlarged_bbox);
	}

	void init_queue(Refinement_pqueue& queue,
		              const FT threshold,
									Iso_cuboid_3& enlarged_bbox)
	{
		Finite_cells_iterator c;
		for(c = m_dt.finite_cells_begin();
			  c != m_dt.finite_cells_end();
				c++)
		{
			FT rer = radius_edge_ratio(c);
			Point point = m_dt.dual(c);
			bool inside = enlarged_bbox.has_on_bounded_side(point);
			if(inside && rer > threshold)
			{
				Vertex_handle v0 = c->vertex(0);
				Vertex_handle v1 = c->vertex(1);
				Vertex_handle v2 = c->vertex(2);
				Vertex_handle v3 = c->vertex(3);
				float score = (float)max_edge_len(c);
				queue.push(Candidate(v0,v1,v2,v3,score));
			}
		}
	}

	void init_nn_search_shell()
	{
		std::list<Vertex_handle> kvertices;
		for(Finite_vertices_iterator v = m_dt.finite_vertices_begin();
        v != m_dt.finite_vertices_end();
        v++)
		{
			if(v->type() != Triangulation::INPUT)
				continue;
			kvertices.push_back(v);
		}
		m_nn_search = K_nearest_neighbor(kvertices.begin(), kvertices.end());
	}

	bool is_refinable(Cell_handle cell,
		                const FT size_shell,
										const FT sizing,
										FT& size,
										Point& p)
	{
		size = circumradius(cell);
		if(size <= sizing)
			return false;

		// try circumcenter
		p = m_dt.dual(cell);
		if(distance_to_input_points(p) < size_shell)
			return true;

		return false;
	}

	FT distance_to_input_points(const Point& p)
	{
		std::list<Vertex_handle> nearest_kvertices;
		Point_vertex_handle_3 query(p.x(),p.y(),p.z());
		m_nn_search.get_k_nearest_neighbors(query,1,nearest_kvertices);

		typename std::list<Vertex_handle>::iterator it = nearest_kvertices.begin();
		Vertex_handle nv = *it;
		if(nv != NULL)
			return distance(nv->point(),p);
		return 0.0; // default
	}

	FT distance(const Point& a, const Point& b)
	{
		return std::sqrt(CGAL::squared_distance(a,b));
	}

	FT max_edge_len(Cell_handle cell)
	{
		const Point& a = cell->vertex(0)->point();
		const Point& b = cell->vertex(1)->point();
		const Point& c = cell->vertex(2)->point();
		const Point& d = cell->vertex(3)->point();
		FT ab = (a-b)*(a-b);
		FT ac = (a-c)*(a-c);
		FT bc = (c-b)*(c-b);
		FT ad = (a-d)*(a-d);
		FT bd = (d-b)*(d-b);
		FT cd = (c-d)*(c-d);
		FT sq_max = (std::max)((std::max)((std::max)(ab,ac),
			          (std::max)(bc,ad)), (std::max)(bd,cd));
		return std::sqrt(sq_max);
	}

	FT circumradius(Cell_handle c)
	{
		Point center = m_dt.dual(c);
		const Point& p = c->vertex(0)->point();
		return std::sqrt((p-center)*((p-center)));
	}

}; // end of Poisson_implicit_function


CGAL_END_NAMESPACE

#endif // CGAL_POISSON_IMPLICIT_FUNCTION_H
