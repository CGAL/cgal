// Copyright (c) 2013-2014  The University of Western Sydney, Australia.
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
// Authors: Weisheng Si, Quincy Tse

/** @file Cone_spanners_2.h
 *
 * This header implements the abstract base class Cone_spanners_2,
 * from which different kinds of cone-based spanner graphs such as
 * Yao graphs and Theta graphs can derive.
 */

#ifndef CGAL_CONE_SPANNERS_2_H
#define CGAL_CONE_SPANNERS_2_H

// if leda::real is used, pls modify the following definition
#define CGAL_USE_CORE 1

#include <CGAL/Cone_spanners_2/_cxx0x_hack.h>

#include <iostream>
#include <cstdlib>
#include <utility>
#include <CGAL/Polynomial.h>
#include <CGAL/number_utils.h>
#include <CGAL/enum.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/Aff_transformation_2.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace CGAL {

/** \ingroup PkgConeBasedSpanners
 * @brief An abstract base class for different cone-based spanner graphs with a given set of
 *  2D points.
 *
 *  Directed,undirected and bidirectional graphs are supported. For differences among these
 *  three types of graphs, please see the documentations of BGL.
 *
 *  The constructor of this class will compute the rays that divide the cones. This computation can
 *  be either inexact by simply dividing Pi by the number of cones (which is quick), or exact by using roots of 
 *  polynomials (entailing number types such as `CORE::Expr` or `LEDA::Real`, which are slow). 
 *  The inexact computation is done by the general class definition of `Cone_spanners_2` below. 
 *  For exact computation, a partial specialization of this class is defined.
 *  If the template parameter `Kernel` is `Exact_predicates_exact_constructions_kernel_with_sqrt`,
 *  this partial specialization class will be invoked. 
 *
 *  In the construction of Yao graph and Theta graph implemented by this package, 
 *  all predicates and construction functions are from \cgal. 
 *  Based on the previous paragraph, if the kernel `Exact_predicates_exact_constructions_kernel_with_sqrt` is used, 
 *  the Yao or Theta graph will be constructed exactly, otherwise inexactly.
 *
 *  Also, in the computation of rays, the direction of the first ray can be specified by passing a parameter
 *  to the constructor, which allows the first ray to start at any direction.
 *
 */
template <typename Kernel, typename Directedness, typename EdgeProperty>
class Cone_spanners_2 {
  public:
    typedef Kernel                          kernel_type;
    typedef Directedness                    directed_selector;

    typedef typename Kernel::Direction_2             Direction_2;
    typedef typename Kernel::Point_2                 Point_2;
    typedef typename Kernel::Line_2                 Line_2;
	typedef typename Kernel::Aff_transformation_2    Transformation;
    typedef boost::adjacency_list<boost::setS,
                                  boost::vecS,
                                  Directedness,
                                  Point_2,
                                  EdgeProperty
                                 > Graph;
	
    /** @brief Constructor.
	 *
     *  Constructs a `Cone_spanners_2` graph object.
     *
     * @param k     Number of cones to divide space into
     * @param start An iterator pointing to the first point (vertex) in the graph.
     *              (default: nullptr)
     * @param end   An iterator pointing to the place that passes the last point. 
	                (default: nullptr)
     * @param ray0  A direction denoting one of the rays deviding the
     *              cones. This allows arbitary rotations of the rays that divide
     *              the plane.
	 *              (default: positive x-axis) 
     */
#ifdef GXX11
    template <typename PointInputIterator=Point_2*>
#else
    template <typename PointInputIterator>
#endif
    Cone_spanners_2 (const unsigned int k,
                    const PointInputIterator& start=nullptr, 
					const PointInputIterator& end=nullptr,
                    const Direction_2& ray0 = Direction_2(1,0)
				   ) 
		           : num_cones(k), g() 
	{
		if (num_cones<2) {
			std::cout << "The number of cones should be larger than 1!" << std::endl;
			std::exit(1);
		}

		populate_vertices(start, end);

		rays.push_back(ray0);
		const double cone_angle = 2*CGAL_PI/num_cones;
		double sin_value, cos_value;
		for (unsigned int i = 1; i < num_cones; i++) {
			sin_value = std::sin(i*cone_angle);
			cos_value = std::cos(i*cone_angle);
			Direction_2 ray_i = Transformation(cos_value, -sin_value, sin_value, cos_value)(ray0);
		    rays.push_back(ray_i);
		}
    }

    /** @brief copy constructor
     *  @param x  another Cone spanner object to copy from.
	 */
    Cone_spanners_2 (const Cone_spanners_2& x) : rays(x.rays), num_cones(x.num_cones), g(x.g) {}

    /** @brief removes all vertices and edges from the graph.
     */
    void clear () {
      g.clear();
    }

    /** @brief inserts a new point into the graph.
     *  
     *  @param p  The point to add.
     */
    void insert (const Point_2& p) 
	{
      g[boost::add_vertex(g)] = p;
    }

    /** @brief inserts the new points contained in the iterator range `[start, end)`
     *  into the graph. 
     *
     *  @param start  The iterator pointing to the first point to be added.
     *  @param end    The iterator pointing to the place just passing the end of the point list.
     *                
     */
#ifdef GXX11
    template <typename PointInputIterator=Point_2*>
#else
    template <typename PointInputIterator>
#endif
    void insert (PointInputIterator start=nullptr,
                 const PointInputIterator& end=nullptr)
    {
      populate_vertices(start, end);
    }

    /** @brief inserts the points in the iterator range `[start, end)`into the graph as vertices.
     *
     * @param start The start iterator.
     * @param end   The end iterator.
     *
     * @return The updated graph.
     */
    template <typename PointInputIterator>
    Graph& populate_vertices(const PointInputIterator& start, const PointInputIterator& end) 
	{
      for (PointInputIterator curr = start; curr != end; ++curr) 
	  {
        g[boost::add_vertex(g)] = *curr;
      }
      return this->g;
    }

    /** @brief returns the cone spanner graph as a `boost::adjacency_list`.
     */
    Graph get_graph() {
      return this->g;
    }

    /** @brief returns the number of cones configured. 
     */
    const unsigned int& get_num_cones() const {
      return num_cones;
    }

    /** @brief returns the vector of directions. 
     */
    const std::vector<Direction_2>& get_directions() const {
      return rays;
    }

    /** casts the cone spanner graph into a `boost::adjacency_list`. 
     */
    operator Graph() const {
      return get_graph();
    }

    /** Function object that orders 2D graph vertex_descriptors based on the "order
     *  induced by the direction D".
     *
     *  This function object is based on the function `CGAL::compare_signed_distance_to_line_2()`, 
	 *  which orders two points according to their signed distance to a line.
     *
     *  @see  `CGAL::compare_signed_distance_to_line_2()`
     */
    struct  vertex_smaller_2
#ifndef GXX11
        : public std::binary_function <typename Graph::vertex_descriptor,
                                       typename Graph::vertex_descriptor, bool>
#endif
    {
      // typedef for C++11 - doesn't hurt to also have for C++98
      typedef typename Graph::vertex_descriptor first_argument_type;
      typedef typename Graph::vertex_descriptor second_argument_type;
      typedef bool     result_type;

      // constructor
	  vertex_smaller_2(const Graph& g, const Direction_2& d)
	  : graph(g), base_line(Point_2(0,0), d) {}

	  // destructor
      ~vertex_smaller_2(){}

      bool operator() (const typename Graph::vertex_descriptor& p,
                       const typename Graph::vertex_descriptor& q) const 
	  {
		  Comparison_result outcome;
          outcome = compare_signed_distance_to_line(base_line, graph[p], graph[q]);
		  if (outcome == SMALLER)
			  return true;
		  else {
			  if (outcome == LARGER)
				  return false;
		  }

          /* otherwise, outcome == CGAL::EQUAL, 
		   *    tie will be broken by a second order according to the ccw90(base_line) direction. */
          // define a rotation of counter clockwise 90
	      Transformation ccw90(0, -1, 1,  0);
	      // rotate 
	      Line_2 ccw90_line = ccw90(base_line);
          outcome = compare_signed_distance_to_line(ccw90_line, graph[p], graph[q]);
		  if (outcome == SMALLER)
			  return true;
		  else 
			  return false;
      }

      private:
        const Graph& graph;
		const Line_2 base_line;
    };

    /** Pure virtual function to be implemented.
	 * Different cone-based spanners will have different implementations for this function.
	 */
	virtual Graph& build_edges() = 0;

#ifndef DOXYGEN_RUNNING
protected:

	/** Store the rays to divide the plane  */
    std::vector<Direction_2>    rays;

	/** Store the number of cones.  */
    const unsigned int  num_cones;

	/** The boost::adjacency_list data structure to store the graph.  */
    Graph g;
#endif

};      // class Cone_spanners_2


/** @brief A partial specialization of `Cone_spanners_2` for exact computation of cones with 
 * `CORE::Expr` (or `leda::real`).
 */
template <typename Directedness, typename EdgeProperty>
class Cone_spanners_2 <Exact_predicates_exact_constructions_kernel_with_sqrt, 
	                   Directedness, 
					   EdgeProperty>
{
  public:
    typedef Exact_predicates_exact_constructions_kernel_with_sqrt   Kernel;
    typedef Kernel                          kernel_type;
    typedef Directedness                    directed_selector;

    typedef typename Kernel::Direction_2             Direction_2;
    typedef typename Kernel::Point_2                 Point_2;
    typedef typename Kernel::Line_2                  Line_2;
    typedef typename Kernel::FT                      FT;
	typedef typename Kernel::Aff_transformation_2    Transformation;

    typedef boost::adjacency_list<boost::setS,
                                  boost::vecS,
                                  Directedness,
                                  Point_2,
                                  EdgeProperty
                                 > Graph;
	
    /** @brief Constructor.
	 *
     *  Constructs a Cone_spanners_2 Graph object.
     *
     * @param k     Number of cones to divide space into
     * @param start An iterator pointing to the first point (vertex) in the graph.
     *              (default: nullptr)
     * @param end   An iterator pointing to the place that passes the last point. 
	                (default: nullptr)
     * @param ray0  A direction denoting one of the rays dividing the
     *              cones. This allows arbitary rotations of the rays 
     *              that divide the plane.
	 *              (default: positive x-axis) 
     */
#ifdef GXX11
    template <typename PointInputIterator=Point_2*>
#else
    template <typename PointInputIterator>
#endif
    Cone_spanners_2 (const unsigned int k,
                    const PointInputIterator& start=nullptr, 
					const PointInputIterator& end=nullptr,
                    const Direction_2& ray0 = Direction_2(1,0)
				   ) 
		           : num_cones(k), g() 
	{
		if (num_cones<2) {
			std::cout << "The number of cones should be larger than 1!" << std::endl;
			std::exit(1);
		}

		//std::cout << "Specialization is called!" << std::endl;

		populate_vertices(start, end);

        // We actually use -x instead of x since CGAL::root_of() will give the k-th
	    //     smallest root but we want the second largest one without counting.
        Polynomial<FT> x(CGAL::shift(Polynomial<FT>(-1), 1));
	    Polynomial<FT> twox(2*x);
	    Polynomial<FT> a(1), b(x);
        for (unsigned int i = 2; i <= num_cones; ++i) {
		     Polynomial<FT> c = twox*b - a;
			 a = b;
			 b = c;
		}
	    a = b - 1;

		unsigned int m, i;
		if (num_cones % 2 == 0)
			m = num_cones/2;
		else 
			m= num_cones/2 + 1;

		FT cos_value, sin_value; 
		Direction_2 ray_i;
		for (i = 1; i <= m; i++) {
            cos_value = - root_of(i, a.begin(), a.end());
	        sin_value = sqrt(FT(1) - cos_value*cos_value);
			ray_i = Transformation(cos_value, -sin_value, sin_value, cos_value)(ray0);
		    rays.push_back(ray_i);
		}

		// add the remaining half number of rays
		if (num_cones % 2 == 0) {
		    for (i = 0; i < m; i++) {
		        rays.push_back(-rays[i]);
		    }
		}
		else {
		    for (i = 0; i < m-1; i++) {
                cos_value = - root_of(m-i, a.begin(), a.end());
	            sin_value = - sqrt(FT(1) - cos_value*cos_value);
			    ray_i = Transformation(cos_value, -sin_value, sin_value, cos_value)(ray0);
		        rays.push_back(ray_i);
		    }
		}
    }

    /** @brief copy constructor.
     *  @param x  another Cone_spanners_2 object to copy from.
	 */
    Cone_spanners_2 (const Cone_spanners_2& x) : rays(x.rays), num_cones(x.num_cones), g(x.g) {}

    /** @brief removes all vertices and edges from the graph.
     */
    void clear () {
      g.clear();
    }

    /** @brief inserts a new point into the graph.

     *  @param p  The point to add.
     */
    void insert (const Point_2& p) 
	{
      g[boost::add_vertex(g)] = p;
    }

    /** @brief inserts the points in the iterator range  `[start, end)`
     *  into the graph. 
     *
     *  @param start  The iterator pointing to the first point to be added.
     *  @param end    The iterator pointing to the place just passing the end of the point list.
     *                
     */
#ifdef GXX11
    template <typename PointInputIterator=Point_2*>
#else
    template <typename PointInputIterator>
#endif
    void insert (PointInputIterator start=nullptr,
                 const PointInputIterator& end=nullptr)
    {
      populate_vertices(start, end);
    }

    /** @brief inserts the points in the iterator range `[start, end)` into the graph as vertices.
     *
     * @param start The start iterator.
     * @param end   The end iterator.
     *
     * @return The updated graph.
     */
    template <typename PointInputIterator>
    Graph& populate_vertices(const PointInputIterator& start, const PointInputIterator& end) 
	{
      for (PointInputIterator curr = start; curr != end; ++curr) 
	  {
        g[boost::add_vertex(g)] = *curr;
      }
      return this->g;
    }

    /** @brief returns the cone spanner graph as a `boost::adjacency_list`.
     */
    Graph get_graph() {
      return this->g;
    }

    /** @brief returns the number of cones configured. 
     */
    const unsigned int& get_num_cones() const {
      return num_cones;
    }

    /** @brief returns the vector of directions.
     */
    const std::vector<Direction_2>& get_directions() const {
      return rays;
    }

    /** Casts the cone_spanner graph into a `boost::adjacency_list`. 
     */
    operator Graph() const {
      return get_graph();
    }

    /** Function object that orders 2D graph vertex_descriptors based on the "order
     *  induced by the direction D".
     *
     *  This function object is based on the function object of directionally_smaller_2 
	 *  which orders two points according to the "order
     *  induced by the direction D".
     *
     *  @see  directionally_smaller_2
     */
    struct  vertex_smaller_2
#ifndef GXX11
        : public std::binary_function <typename Graph::vertex_descriptor,
                                       typename Graph::vertex_descriptor, bool>
#endif
    {
      // typedef for C++11 - doesn't hurt to also have for C++98
      typedef typename Graph::vertex_descriptor first_argument_type;
      typedef typename Graph::vertex_descriptor second_argument_type;
      typedef bool     result_type;

      // constructor
	  vertex_smaller_2(const Graph& g, const Direction_2& d)
	  : graph(g), base_line(Point_2(0,0), d) {}

	  // destructor
      ~vertex_smaller_2(){}

      bool operator() (const typename Graph::vertex_descriptor& p,
                       const typename Graph::vertex_descriptor& q) const 
	  {
		  Comparison_result outcome;
          outcome = compare_signed_distance_to_line(base_line, graph[p], graph[q]);
		  if (outcome == SMALLER)
			  return true;
		  else {
			  if (outcome == LARGER)
				  return false;
		  }

          /* otherwise, outcome == CGAL::EQUAL, 
		   *    tie will be broken by a second order according to the ccw90(base_line) direction. */
          // define a rotation of counter clockwise 90
	      Transformation ccw90(0, -1, 1,  0);
	      // rotate 
	      Line_2 ccw90_line = ccw90(base_line);
          outcome = compare_signed_distance_to_line(ccw90_line, graph[p], graph[q]);
		  if (outcome == SMALLER)
			  return true;
		  else 
			  return false;
      }

      private:
        const Graph& graph;
		const Line_2 base_line;
    };

    /** Pure virtual function to be implemented.
	 * Different cone-based spanners will have different implementations for this function.
	 */
	virtual Graph& build_edges() = 0;

protected:

	/** Store rays to divide the plane  */
    std::vector<Direction_2>    rays;

	/** Indicate the number of cones.  */
    const unsigned int  num_cones;

	/** The boost::adjacency_list data structure to store the graph.  */
    Graph g;

};      // end of specialization: Cone_spanners_2

/* serialize, to be implemented in future
template < typename Kernel, typename Directedness, typename EdgeProperty >
std::istream& operator>> (std::istream& is, Cone_spanners_2<Kernel, Directedness, EdgeProperty>& cone_spanner);

template < typename Kernel, typename Directedness, typename EdgeProperty > 
std::ostream& operator<< (std::ostream& os, const Cone_spanners_2<Kernel, Directedness, EdgeProperty>& cone_spanner);
*/

}  // namespace CGAL

#ifdef GXX11
#undef GXX11
#endif

#endif
