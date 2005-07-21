// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Efi Fogel          <efif@post.tau.ac.il>

/*!
 * Cubical_gaussian_map is a data dtructure that represents a Gaussinal map
 * on a cube. The representation consists of 6 arrangements that correspond to
 * the 6 faces of the unit cube (an axes parallel cube whose edges are of length
 * 2 centered at the origin).
 *
 * The arrangements are indexed 0 to 5, where 0, 1, and 2 refer to the
 * arrangements of the faces that reside in the negative X, Y, and Z halfspaces
 * respectively. In a similar way 3, 4, and 5 refer to the arrangements of the
 * faces that reside in the positive X, Y, and Z halfspaces respectively.

 * When reducing the 3D coordinate system onto a 2D coordinate system attached
 * to a paerticular cubic face, we project the 3D coordinate system as follows:
 *
 * Cubic Face Projection
 * ---------- ----------
 * 0          X,Y,Z => Y,Z
 * 1          X,Y,Z => Z,X
 * 2          X,Y,Z => X,Y
 * 3          X,Y,Z => Z,Y
 * 4          X,Y,Z => X,Z
 * 5          X,Y,Z => Y,X
 *
 * As mentioned above each arrangement is bounded by a 2-unit square centered at
 * the origin of a 2D Cartesian coordinate system. Each 2D vertex and 2D edge of
 * the boundary is indexed as depicted below:
 *
 * 1:(-1,1)           3:(1,1)
 *     ----------------|
 *     |       1       |
 *     |               |
 *     |               |
 *     |0             3|
 *     |               |
 *     |               |
 *     |       2       |
 *     |----------------
 * 0:(-1,-1)         2:(1,-1)
 *
 *
 *                           x
 *                           ^
 *                           |
 *                           |
 *           ----------------|-->y
 *      y   /0             1/|
 *      ^  /        4      /2|
 *      | /               /  |
 *      |/2             3/   |
 *      |----------------3   |
 *     /|1             3|    |
 *   |/_|               |  3 |
 *  x   |               |    |
 *      |        5      |   0/
 *      |               |   /
 *      |               |  /
 *      |               |1/
 *      |0             2|/
 *      ----------------|-->x
 *                     /
 *                   |/_
 *                  y
 *    
 *          y,x
 *           ^
 *           |
 *           |
 *           |----------------
 *          /|2             3|
 *         /1|               |
 *        /  |               |
 *       /   |       2       |
 *      /3   |               |
 *      |    |               |
 *      |  0 |               |
 *      |    |0             1|
 *      |   0/------------------>y,x
 *      |   /0             2/
 *      |  /        1      /
 *      |2/               /
 *      |/1             3/
 *      |---------------/
 *     /
 *   |/_
 *  x,y
 *
 *
 * The table below lists the adjacent halfedges on adjacent cube faces
 *
 * Face | Vertex | ad. face | ad. vertex
 * -------------------------------------
 *   0  |    0   |   1      |    1
 *   0  |    2   |   5      |    1
 *   1  |    0   |   2      |    1
 *   1  |    2   |   3      |    1
 *   2  |    0   |   0      |    1
 *   2  |    2   |   4      |    1
 *   3  |    2   |   4      |    3
 *   3  |    0   |   2      |    3
 *   4  |    2   |   5      |    3
 *   4  |    0   |   0      |    3
 *   5  |    2   |   3      |    3
 *   5  |    0   |   1      |    3
 *
 *
 * The table below lists the incident edges of each face as indexed in the
 * 2 incident faces. For example, the common edge of faces 0 and 1 is index
 * 2 in face 0 and 0 in face 1
 *
 * Face |  0  |  1  |  2  |  3  |  4  |  5
 * ---- ------------------------------------
 *  0   | -   | 2,0 | 0,3 | -   | 1,2 | -
 *  1   | 0,2 | -   | 2,0 | 3,0 | -   | 1,2
 *  2   | 3,0 | 0,2 | -   | 1,2 | 2,0 | -
 *  3   | -   | 0,3 | 1,2 | -   | 3,1 | 1,3
 *  4   | 2,1 | -   | 0,2 | 1,3 | -   | 3,1
 *  5   | 0,3 | 2,1 | -   | 3,1 | 1,3 | -
 *
 * The table below lists the circular transition of corner vertices. Each pair
 * of indices consists of a arrangement id followed by a vertex id.
 * 
 * (0,0) => (1,0) | (1,0) => (2,0) | (2,0) => (0,0)
 * (0,1) => (2,2) | (1,1) => (0,2) | (2,1) => (1,2)
 * (0,2) => (5,0) | (1,2) => (3,0) | (2,2) => (4,0)
 * (0,3) => (4,2) | (1,3) => (5,2) | (2,3) => (3,2)
 *
 * (3,0) => (2,1) | (4,0) => (0,1) | (5,0) => (1,1)
 * (3,1) => (1,3) | (4,1) => (2,3) | (5,1) => (0,3)
 * (3,2) => (4,1) | (4,2) => (5,1) | (5,2) => (3,1)
 * (3,3) => (5,3) | (4,3) => (3,3) | (5,3) => (4,3)
 *
 * Let cur_arr_id, cur_ver_id demote the ids of a arrangement and corner vertex.
 * Then, the following procedure compute the transion:
 */

#ifndef CGAL_CUBICAL_GAUSSIAN_MAP_3_H
#define CGAL_CUBICAL_GAUSSIAN_MAP_3_H

#define CGAL_ARR_SEGMENT_TRAITS                 0
#define CGAL_ARR_NON_CACHING_SEGMENT_TRAITS     1

#include <CGAL/Cartesian.h>
#include <CGAL/Arrangement_2.h>
// #include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/HalfedgeDS_face_base.h>
#include <CGAL/Point_3.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include "CGAL/Cgm_arr_dcel.h"
#if CGAL_CGM_TRAITS == CGAL_ARR_SEGMENT_TRAITS
#include <CGAL/Arr_segment_traits_2.h>
#elif CGAL_CGM_TRAITS == CGAL_ARR_NON_CACHING_SEGMENT_TRAITS
#include <CGAL/Arr_non_caching_segment_traits_2.h>
#endif

#include "CGAL/Cgm_plane_3.h"

CGAL_BEGIN_NAMESPACE

template <class FT, class Point_2>
class Cgm_normalize {
public:
  void operator()(Point_2 & p) {}
};

#if defined(CGAL_USE_LEDA)
template <class Point_2>
class Cgm_normalize<leda::rational, Point_2> {
public:
  void operator()(Point_2 & p)
  {
    leda::rational x = p.x();
    x.normalize();
    leda::rational y = p.y();
    y.normalize();
    p = Point_2(x, y);
  }
};
#endif

template <class T_Kernel,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T>
#endif
          class T_Dcel = Cgm_arr_dcel>
class Cubical_gaussian_map_3 {
private:
  typedef Cubical_gaussian_map_3<T_Kernel, T_Dcel>        Self;

public:
  // Arrangement traits and types:
#if CGAL_CGM_TRAITS == CGAL_ARR_SEGMENT_TRAITS
  typedef CGAL::Arr_segment_traits_2<T_Kernel>          Arr_traits;
#else
  typedef CGAL::Arr_non_caching_segment_traits_2<T_Kernel>   Arr_traits;
#endif
  
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
  typedef T_Dcel<Arr_traits>                            Arr_dcel;
#else
  typedef typename T_Dcel::template Dcel<Arr_traits>    Arr_dcel;
#endif

  typedef T_Kernel                                      Kernel;

private:
  typedef CGAL::Arrangement_2<Arr_traits,Arr_dcel>      Arr;

  // Point_location;
  // typedef typename CGAL::Arr_walk_along_line_point_location<Arr>
  typedef typename CGAL::Arr_naive_point_location<Arr>  Point_location;

  typedef CGAL::Cgm_plane_3<Kernel>                     Cgm_plane_3;
public:
  
  typedef typename Kernel::FT                           FT;
  typedef typename Kernel::Point_3                      Point_3;
  typedef typename Kernel::Vector_3                     Vector_3;
  typedef typename Kernel::Plane_3                      Plane_3;

  typedef typename Arr_traits::Point_2                  Point_2;
  typedef typename Arr_traits::Curve_2                  Curve_2;
  typedef typename Arr_traits::X_monotone_curve_2       X_monotone_curve_2;

  enum {NUM_CORNERS = 4};
  enum {NUM_EDGES = 12};
  
  /*! The indices of the 6 faces */
  enum {MINX = 0, MINY, MINZ, MAXX, MAXY, MAXZ, NUM_FACES};

  typedef Arr                                           Arrangement;
  typedef typename Arr::Vertex_iterator                 Arr_vertex_iterator;
  typedef typename Arr::Vertex_const_iterator
    Arr_vertex_const_iterator;
  typedef typename Arr::Halfedge_iterator               Arr_halfedge_iterator;
  typedef typename Arr::Halfedge_const_iterator
    Arr_halfedge_const_iterator;
  typedef typename Arr::Face_iterator                   Arr_face_iterator;
  typedef typename Arr::Face_const_iterator             Arr_face_const_iterator;
  typedef typename Arr::Edge_iterator                   Arr_edge_iterator;
  typedef typename Arr::Edge_const_iterator             Arr_edge_const_iterator;

  typedef typename Arr::Halfedge_around_vertex_circulator
    Arr_halfedge_around_vertex_circulator;
  typedef typename Arr::Halfedge_around_vertex_const_circulator
    Arr_halfedge_around_vertex_const_circulator;
  typedef typename Arr::Ccb_halfedge_circulator
    Arr_ccb_halfedge_circulator;
  typedef typename Arr::Ccb_halfedge_const_circulator
    Arr_ccb_halfedge_const_circulator;
  typedef typename Arr::Holes_iterator                  Arr_holes_iterator;
  
  typedef typename Arr::Vertex_handle                   Arr_vertex_handle;
  typedef typename Arr::Vertex_const_handle             Arr_vertex_const_handle;
  typedef typename Arr::Halfedge_handle                 Arr_halfedge_handle;
  typedef typename Arr::Halfedge_const_handle
    Arr_halfedge_const_handle;
  typedef typename Arr::Face_handle                     Arr_face_handle;
  typedef typename Arr::Face_const_handle               Arr_face_const_handle;

  typedef typename Arr::Vertex                          Arr_vertex;
  typedef typename Arr_vertex::Vertex_location          Vertex_location;
  
public:
  /*! This class represents a normal projected onto the unit cube */
  class Projected_normal {
  private:
    /*! The central projection of the normal onto the unit cube */
    Point_3 m_central_projection;

    /*! The arrangement vertex handle of the projected noraml. When a normal is
     * projected onto the a of a unit-cube face, it is projected onto more than
     * one face. A normal can be projected onto at most three faces of the unit
     * cube (when it is projected onto a corner). This is not space-efficient,
     * because most normals project onto the interior of a single face
     */
    Arr_vertex_handle m_vertices[3];
    
    /*! Indicates what faces of the box the normal project onto */ 
    unsigned int m_faces_mask;

    /*! Indicates what arrangement vertices have been set already */ 
    unsigned int m_vertices_mask;

    /*! The number of faces of the box the normal projects onto (1,2, or 3) */
    unsigned int m_num_faces;

    /*! Indicates whether the vertex is real */
    bool m_is_real;

  public:
    /*! Parameter-less constructor */
    Projected_normal() :
      m_faces_mask(0x0), m_vertices_mask(0x0), m_num_faces(0),
      m_is_real(false) {}

    /*! Constructor */
    Projected_normal(const Point_3 & central_projection,
                     unsigned int faces_mask, unsigned int num_faces) :
      m_central_projection(central_projection),
      m_faces_mask(faces_mask), m_vertices_mask(0x0), m_num_faces(num_faces),
      m_is_real(false)
    {}

    /*! Constructor */
    Projected_normal(const Vector_3 & normal) :
      m_vertices_mask(0x0), m_is_real(false)
    {
      compute_projection(normal);
    }
    
    /*! Obtain the projected normal onto the box faces */
    const Point_3 & get_projected_normal() const
    { return m_central_projection; }

    /*! Obtain the mask that indicates what faces of the box the nromal
     * projects onto
     */
    unsigned int get_faces_mask() const { return m_faces_mask; }

    /*! Return true if a arrangement vertex has been set already */
    bool is_vertex_set(unsigned int i) const
    {
      return m_vertices_mask & get_mask(i);
    }
    
    /*! sets a arrangement vertex */
    void set_vertex(Arr_vertex_handle vertex, unsigned int i)
    {
      m_vertices[i] = vertex;
      m_vertices_mask |= get_mask(i);
    }

    /*! obtains the vertices mask  */
    unsigned int get_vertices_mask() const { return m_vertices_mask; }

    /*! obtains a arrangement vertex */
    Arr_vertex_handle & get_vertex(unsigned int i) { return m_vertices[i]; }

    /*! \brief obtains the number of box faces the normal projects onto */
    unsigned int get_num_faces() const { return m_num_faces; }

    /*! Obtain the flag that indicates whether the projected normal is real */
    bool get_is_real() const { return m_is_real; }
    
    /*! Set the flag that indicates whether the projected normal is real */
    void set_is_real(bool value) { m_is_real = value; }

    /*! Computes the projection of a normal vector onto the unit cube..
     * \param normal the vector to project 
     *
     *  This function computes the central projection of a normal onto the unit
     * cube. In particular, it computes the 3D projected normal onto the unit
     * cube. If the projected normal coincides with a corner vertex of the cube,
     * it is included in 3 cube faces. Otherwise, if it coincides with an edge
     * of the cube, it is included in 2 cube faces. Otherwise, it is included in
     * 1 cube face. We record the number of unit-cube faces the projected
     * normal is included in, and a mask that represents the including unit-cube
     * faces. 
     */
    void compute_projection(const Vector_3 & normal)
    {
      unsigned int i;
      for (i = 0; i < 3; ++i) {
        unsigned int m;
    
        if (normal[i] < 0) m = 0;
        else if (normal[i] > 0) m = 1;
        else continue;
    
        unsigned int j = (i + 1) % 3;
        unsigned int k = (i + 2) % 3;
        unsigned int m_i = 3 * m + i;

        // Store the central projection of the normal:
        FT vec[3];
        vec[i] = get_extreme_coordinate(m);
        FT scale = vec[i] / normal[i];
        vec[j] = normal[j] * scale;
        vec[k] = normal[k] * scale;
        if ((get_extreme_coordinate(0) > vec[j]) ||
            (vec[j] > get_extreme_coordinate(1)) ||
            (get_extreme_coordinate(0) > vec[k]) ||
            (vec[k] > get_extreme_coordinate(1)))
          continue;

        m_central_projection = Point_3(vec[0], vec[1], vec[2]);
        m_faces_mask = get_mask(m_i);
        m_num_faces = 1;

        // Handle degenerate cases:
        if (get_extreme_coordinate(0) == vec[j]) {
          m_faces_mask |= get_mask(j);
          m_num_faces++;
        } else if (get_extreme_coordinate(1) == vec[j]) {
          m_faces_mask |= get_mask(3+j);
          m_num_faces++;
        }
        if (get_extreme_coordinate(0) == vec[k]) {
          m_faces_mask |= get_mask(k);
          m_num_faces++;
        } else if (get_extreme_coordinate(1) == vec[k]) {
          m_faces_mask |= get_mask(3+k);
          m_num_faces++;
        }
        break;
      }
    }
  };

public:
  // Arrangement public methods:
  
  /*! Obtain the (const) arrangement associated with a unit-cube face
   * \param id the id of the unit-cube face
   */
  const Arrangement & get_arrangement(unsigned int id) const
  { return m_arrangements[id]; }

  /*! Obtain the (mutable) arrangement associated with a unit-cube face
   * \param id the id of the unit-cube face
   */
  Arrangement & get_arrangement(unsigned int id) { return m_arrangements[id]; }
  
  /*! Calculate the normal vector that projects onto the given point of the
   * given unit-cube face
   * \param id the index of the unit cube face
   * \param p the 2D point that the normal projects onto
   * \return the normal vector
   */
  Vector_3 calculate_normal(unsigned int id, Point_2 & p) const
  {
    switch (id) {
     case 0: return Vector_3(get_extreme_coordinate(0), p.y(), p.x()); break;
     case 1: return Vector_3(p.x(), get_extreme_coordinate(0), p.y()); break;
     case 2: return Vector_3(p.y(), p.x(), get_extreme_coordinate(0)); break;
     case 3: return Vector_3(get_extreme_coordinate(1), p.x(), p.y()); break;
     case 4: return Vector_3(p.y(), get_extreme_coordinate(1), p.x()); break;
     case 5: return Vector_3(p.x(), p.y(), get_extreme_coordinate(1)); break;
    }
    return Vector_3();
  }
  
  /*! Calculate the 2D point, which is the projection of a normal, and is
   * contained in the arrangement associated with a given unit-cube face
   * \param projected_normat the projected normal
   * \param id the id of the unit-cube face
   */
  Point_2 construct_point(const Point_3 & p3, unsigned int id) const
  {
    unsigned int j = (id + 1) % 3;
    unsigned int k = (id + 2) % 3;    

    Point_2 p2 = (id < 3) ? Point_2(p3[k], p3[j]) : Point_2(p3[j], p3[k]);
    Cgm_normalize<FT, Point_2> normalize;
    normalize(p2);
    return p2;
  }

  /*! A circulator over the real halfedges incident to a real vertex.
   * This is an adapter that traverses the incident arrangement halfedges
   * incident to a given vertex, skipping the non-real halfedges, and jumping
   * from one unit-cube face to the adjacent face, when the given vertex is
   * a corner vertex or lies on an edge.
   */
  class Halfedge_around_vertex_const_circulator :
    public Arr_halfedge_around_vertex_const_circulator
  {
  private:
    typedef Halfedge_around_vertex_const_circulator             Self;
    typedef Arr_halfedge_around_vertex_const_circulator         Base;

    /*! \brief advances the base pointer until it points to a real halfedge */
    void advance()
    {
      Base * base = this;
      while (!(*base)->get_is_real() || (*base)->face()->is_unbounded()) {
        // Move to the halfedge of the neighboring arrangement, :
        if ((*base)->face()->is_unbounded()) {
          void * tmp = (*base)->target()->get_adjacent_vertex();
          Arr_vertex_const_handle vh = *((Arr_vertex_const_handle *) (&tmp));
          
          // Find the first halfedge after the one incident to the unbounded
          // face:
          *base = vh->incident_halfedges();
          Arr_halfedge_around_vertex_const_circulator begin_hec = *base;
          do {
            if ((*base)->face()->is_unbounded()) break;
            ++*base;
          } while ((*base) != begin_hec);
          ++*base;
        }

        // Increment the halfedge if it is artificial
        if (!(*base)->get_is_real()) {
          ++*base;
        }
      }
    }
  
  public:
    /*! Default constructor */
    Halfedge_around_vertex_const_circulator() : Base() {}

    /*! Constructor */
    Halfedge_around_vertex_const_circulator(Base base) : Base(base)
    {
      advance();
    }

    /*! pre-increment operator */
    Self & operator++()
    {
      CGAL_assertion(!(*this)->face()->is_unbounded());
      CGAL_assertion((*this)->get_is_real());
    
      // ++Base::(*this);
      Base * base = this;
      ++*base;
      advance();
      return *this;
    }

    /*! Post increment */
    Self operator++(int)
    {
      Self tmp = *this;
      ++*this;
      return tmp;
    }

    /*! Equality */
    bool operator==(const Self & circulator)
    {
      Base * base_p = this;
      return (*base_p == circulator);;
    }

    /*! Inequality */
    bool operator!=(const Self & circulator)
    {
      Base * base_p = this;
      return (*base_p != circulator);;
    }

    /*! Assignment operator from Base */
    Self & operator=(Base & other_base)
    {
      Base * base = this;
      *base = other_base;
      advance();
    }
  };
  
  /*! Return the mask of a face given by the face id */
  static unsigned int get_mask(unsigned int id)
  {
    /*! A mask for each face */
    static unsigned int s_mask[NUM_FACES] = {
      0x1 << MINX, 0x1 << MINY, 0x1 << MINZ,
      0x1 << MAXX, 0x1 << MAXY, 0x1 << MAXZ
    };
    return s_mask[id];
  }

protected:
  typedef unsigned int * Coord_index_iter;

  /*! Return an extreme (1 or -1) coordinate of the unit cube in the user
   * defined number-type.
   * \param i an index.
   * \return -1 if id is 0, and 1 if id is 1 depending on whether i is equal to
   * 0 or 1 respectively.
   */
  static const FT & get_extreme_coordinate(unsigned int i)
  {
    /*! -1,+1 - the unit-cube extreme values */
    static FT s_extreme[2] = {-1, 1};
    return s_extreme[i];
  }

  /*! Generate and return the corner point given by its index according to the
   * chart below.
   * The index is specified on the left, and the coordinates on the right
   * 1:(-1,1)           3:(1,1)
   *     ----------------|
   *     |       1       |
   *     |               |
   *     |               |
   *     |0             3|
   *     |               |
   *     |               |
   *     |       2       |
   *     |----------------
   * 0:(-1,-1)         2:(1,-1)
   * \param index the (local) index of the point
   * \return the corner point
   */
  static const Point_2 & get_corner_point(unsigned int index)
  {
    /*! The corners */
    static Point_2 s_corners[NUM_CORNERS] = {
      Point_2(get_extreme_coordinate(0), get_extreme_coordinate(0)),
      Point_2(get_extreme_coordinate(0), get_extreme_coordinate(1)),
      Point_2(get_extreme_coordinate(1), get_extreme_coordinate(0)),
      Point_2(get_extreme_coordinate(1), get_extreme_coordinate(1))
    };
    return s_corners[index];
  }

  /*! Obtain the index of a given corner point
   * \param point the corner point
   */
  static unsigned int get_corner_index(Point_2 & point)
  {
    for (unsigned int i = 0; i < NUM_CORNERS; ++i)
      if (get_corner_point(i) == point) return i;
    CGAL_assertion(0);
    return static_cast<unsigned int>(-1);
  }
  
  /*! A specific corner-vertex is identified by the id of the unit-cube face
   * associated with the arrangement containing the vertex, and the id of the
   * corner within that face.
   */
  typedef std::pair<unsigned int, unsigned int> Corner_id;

  /*! Obtain the handle of a corner vertex given by its id
   * \param the id of the corner vertex
   * \return the handle of the corner vertex
   */
  Arr_vertex_handle get_corner_vertex_handle(Corner_id corner_id)
  {
    return m_corner_vertices[corner_id.first][corner_id.second];
  }
  
  /*! Obtain the id of the next corner vertex in the cyclic chain of corner
   * vertices.
   * \param face_id the id of the unit-cube face (projection)
   * \param corner_index the index of the corner within the cube face
   */
  static Corner_id get_next_corner_id(Corner_id & corner_id)
  { return get_next_corner_id(corner_id.first, corner_id.second); }
  
  /*! Obtain the id of the next corner vertex in the cyclic chain of corner
   * vertices.
   * \param face_id the id of the unit-cube face (projection)
   * \param corner_index the index of the corner within the cube face
   */
  static Corner_id get_next_corner_id(unsigned int face_id,
                                      unsigned int corner_index)
  {
    /*! The next corner */
    static Corner_id s_next_corner[NUM_FACES][NUM_CORNERS] = {
      {Corner_id(1,0), Corner_id(2,2), Corner_id(5,0), Corner_id(4,2)},
      {Corner_id(2,0), Corner_id(0,2), Corner_id(3,0), Corner_id(5,2)},
      {Corner_id(0,0), Corner_id(1,2), Corner_id(4,0), Corner_id(3,2)},
      {Corner_id(2,1), Corner_id(1,3), Corner_id(4,1), Corner_id(5,3)},
      {Corner_id(0,1), Corner_id(2,3), Corner_id(5,1), Corner_id(3,3)},
      {Corner_id(1,1), Corner_id(0,3), Corner_id(3,1), Corner_id(4,3)}
    };
    return s_next_corner[face_id][corner_index];
  }

  /*! Obtain the id of one of the two corner vertices that are the target
   * vertices of 2 boundary halfedges that represent the same unit-cube edge.
   * Each edge is associated with two halfedges on the (outer) connected
   * component boundary (CCB) of two adjacent unit-cube faces respectively.
   * The corner vertex whose id is returned is the target vertex of such a
   * halfedge.
   * \param edge_index the index of one of the dozen edges of the unit-cube.
   * \param index the (arbitrary) index that indicates which one of the two
   * incident vertices is seeked (0 and 1 indicate first and second
   * respectively).
   * \return the id of the corner vertex
   */
  static Corner_id get_incident_corner_id(unsigned int edge_index,
                                          unsigned int index)
  {
    /*! The adjacent halfedges */
    static Corner_id s_adjacent_halfedges[NUM_EDGES][2] = {
      {Corner_id(0, 2), Corner_id(1, 0)},
      {Corner_id(0, 3), Corner_id(5, 0)},
      {Corner_id(1, 2), Corner_id(2, 0)},
      {Corner_id(1, 3), Corner_id(3, 0)},
      {Corner_id(2, 2), Corner_id(0, 0)},
      {Corner_id(2, 3), Corner_id(4, 0)},
      {Corner_id(3, 3), Corner_id(4, 1)},
      {Corner_id(3, 2), Corner_id(2, 1)},
      {Corner_id(4, 3), Corner_id(5, 1)},
      {Corner_id(4, 2), Corner_id(0, 1)},
      {Corner_id(5, 3), Corner_id(3, 1)},
      {Corner_id(5, 2), Corner_id(1, 1)}
    };
    return s_adjacent_halfedges[edge_index][index];
  }

  Arr m_arrangements[NUM_FACES];

  /*! The corner vertices */
  Arr_vertex_handle m_corner_vertices[NUM_FACES][NUM_CORNERS];

  /*! Convert number of incident faces to vertex location
   * \param num_faces the number of faces
   */
  Vertex_location num_faces_2_vertex_location(unsigned int num_faces)
  {
    CGAL_assertion(num_faces);
    return (num_faces == 1) ? Arr_vertex::Interior :
      ((num_faces == 2) ? Arr_vertex::Edge : Arr_vertex::Corner);
  }

  /*! Convert a  face mask to the index of the face resresented by the mask
   * \param mask the face mask
   * \return the face index
   * \pre the mask represents a single face. That is, only one bit among the 6
   * least significant bits is turned on
   */
  unsigned int face_mask_2_face_index(unsigned int mask)
  {
    unsigned int id;
    for (id = 0; id < NUM_FACES; id++)
    if (mask & get_mask(id)) return id;
    return (unsigned int) -1;
  }
  
  /*! Splits an edge of a arrangement. The embedding of the edge is a curve
   * containing a given point. The curve is split into two sub-curves, and
   * re-inserted into the arrangement.
   * \param arr_id the arrangement identifier
   * \param p the split point
   * \return handle to a halfedge e. The embedding of e is one of the
   * sub-curves, and the embedding of the target vertex of e is the given
   * point p.
   * \pre p is on the boundary of the box
   */
  Arr_halfedge_handle split(unsigned int id, const Point_2 & p,
                            Projected_normal & dual)
  {
    Arr & arr = m_arrangements[id];
    Point_location pl(arr);
    CGAL::Object obj = pl.locate(p);
    if (const Arr_vertex_const_handle * const_vertex_ptr =
        CGAL::object_cast<Arr_vertex_const_handle>(&obj))
    {
      Arr_vertex_handle vertex = arr.non_const_handle(*const_vertex_ptr);
      return vertex->incident_halfedges();
    }
    else if (const Arr_halfedge_const_handle * const_edge_ptr =
             CGAL::object_cast<Arr_halfedge_const_handle>(&obj))
    {
      Arr_halfedge_handle edge = arr.non_const_handle(*const_edge_ptr);
      const X_monotone_curve_2 & curve = edge->curve();
      X_monotone_curve_2 cv1, cv2;
      Arr_traits traits;
      traits.split_2_object()(curve, p, cv1, cv2);
      /*! \todo
       * edge->set_is_real(false);
       * edge>twin()->set_is_real(false);
       */
      Arr_halfedge_handle he = arr.split_edge(edge, cv1, cv2);
      return he;
    }
    CGAL_assertion_msg(0, "Failed to locate an edge!");
    return arr.non_const_handle(Arr_halfedge_const_handle());
  }

  /*! Inserts a curve (segment) into a arrangement between 2 projected normals.
   * The endpoints of the segment are central projections.
   * \param proj_normal1 represents the source projected normal
   * \param proj_normal2 represents the target projected normal
   * \param id the index of the unit-cube face, which is also the index of the
   * arrangement to insert the segment into
   * \param boundary_flag indicates whether insert() should be called
   * recursively in case the segment to be added lies on the boundary. A
   * recursive call is invoced with the flag turned on. This indicates that the
   * recursion should stop.
   */
  void insert(Projected_normal & proj_normal1, Projected_normal & proj_normal2,
              unsigned int id, bool unique = false)
  {
    unsigned int i = (id + 0) % 3;
    
    unsigned int faces_mask1 = proj_normal1.get_faces_mask();
    unsigned int faces_mask2 = proj_normal2.get_faces_mask();

#if 0
    std::cout << "inserting:" << std::endl
              << "projected normal 1: "
              << proj_normal1.get_projected_normal() << ", 0x"
              << std::hex << proj_normal1.get_faces_mask() << ", "
              << proj_normal1.get_num_faces() << std::endl
              << "projected normal 2: "
              << proj_normal2.get_projected_normal() << ", 0x"
              << std::hex << proj_normal2.get_faces_mask() << ", "
              << proj_normal2.get_num_faces() << std::endl
              << "id: " << id << std::endl;
#endif
    
    // construct the 2d point that is the projection of normal onto the face:
    const Point_3 & point1 = proj_normal1.get_projected_normal();
    Point_2 p1 = construct_point(point1, id);

    // Construct the 2d point that is the projection of normal onto the face:
    const Point_3 & point2 = proj_normal2.get_projected_normal();
    Point_2 p2 = construct_point(point2, id);

    /* If the projected normals are not uniqe, we need to locate the vertices
     * in case they already exists.
     */
    if (!unique) {
      Arr & arr = m_arrangements[id];
      Point_location pl(arr);
      if (!proj_normal1.is_vertex_set(i)) {
        if (proj_normal1.get_num_faces() == 3) {
          // Set 1st:
          Corner_id corner_id(id, get_corner_index(p1));
          proj_normal1.set_vertex(get_corner_vertex_handle(corner_id),
                                  corner_id.first % 3);
          // Set 2nd:
          corner_id = get_next_corner_id(corner_id);
          proj_normal1.set_vertex(get_corner_vertex_handle(corner_id),
                                  corner_id.first % 3);
          // Set 3rd:
          corner_id = get_next_corner_id(corner_id);
          proj_normal1.set_vertex(get_corner_vertex_handle(corner_id),
                                  corner_id.first % 3);
        } else if (proj_normal1.get_num_faces() == 2) {
          CGAL::Object obj = pl.locate(p1);
          if (const Arr_vertex_const_handle * const_vertex_ptr =
              CGAL::object_cast<Arr_vertex_const_handle>(&obj))
          {
            Arr_vertex_handle vertex = arr.non_const_handle(*const_vertex_ptr);
            proj_normal1.set_vertex(vertex, i);
            // Adjacent face:
            unsigned int adj_id =
              face_mask_2_face_index(faces_mask1 & !get_mask(id));
            void * tmp = vertex->get_adjacent_vertex();
            Arr_vertex_handle adj_vertex = *((Arr_vertex_handle *) (&tmp));
            proj_normal1.set_vertex(adj_vertex, adj_id % 3);
          }
        } else {
          CGAL::Object obj = pl.locate(p1);
          if (const Arr_vertex_const_handle * const_vertex_ptr =
              CGAL::object_cast<Arr_vertex_const_handle>(&obj))
          {
            Arr_vertex_handle vertex = arr.non_const_handle(*const_vertex_ptr);
            proj_normal1.set_vertex(vertex, i);
          }
        }
      }

      if (!proj_normal2.is_vertex_set(i)) {
        if (proj_normal2.get_num_faces() == 3) {
          // Set 1st:
          Corner_id corner_id(id, get_corner_index(p2));
          proj_normal2.set_vertex(get_corner_vertex_handle(corner_id),
                                  corner_id.first % 3);
          // Set 2nd:
          corner_id = get_next_corner_id(corner_id);
          proj_normal2.set_vertex(get_corner_vertex_handle(corner_id),
                                  corner_id.first % 3);
          // Set 3rd:
          corner_id = get_next_corner_id(corner_id);
          proj_normal2.set_vertex(get_corner_vertex_handle(corner_id),
                                  corner_id.first % 3);
        } else if (proj_normal2.get_num_faces() == 2) {
          CGAL::Object obj = pl.locate(p2);
          if (const Arr_vertex_const_handle * const_vertex_ptr =
              CGAL::object_cast<Arr_vertex_const_handle>(&obj))
          {
            Arr_vertex_handle vertex = arr.non_const_handle(*const_vertex_ptr);
            proj_normal2.set_vertex(vertex, i);
            // Adjacent face:
            unsigned int adj_id =
              face_mask_2_face_index(faces_mask2 & !get_mask(id));
            void * tmp = vertex->get_adjacent_vertex();
            Arr_vertex_handle adj_vertex = *((Arr_vertex_handle *) (&tmp));
            proj_normal2.set_vertex(adj_vertex, adj_id % 3);
          }
        } else {
          CGAL::Object obj = pl.locate(p2);
          if (const Arr_vertex_const_handle * const_vertex_ptr =
              CGAL::object_cast<Arr_vertex_const_handle>(&obj))
          {
            Arr_vertex_handle vertex = arr.non_const_handle(*const_vertex_ptr);
            proj_normal2.set_vertex(vertex, i);
          }
        }
      }
    }

    /* If the computed point is on the boundary of the box face, but not a
     * corner point, locate it. If the locate returns an edge, split it at the
     * point:
     */
    if (proj_normal1.get_num_faces() == 2) {
      if (!proj_normal1.is_vertex_set(i)) {
        Arr_halfedge_handle edge = split(id, p1, proj_normal1);
        const Arr_vertex_handle & v = edge->target();
        proj_normal1.set_vertex(v, i);
        v->set_face_id(id);
        v->set_location(Arr_vertex::Edge);
        if (proj_normal1.get_is_real()) v->set_is_real(true);
        
        // Find adjacent face:
        unsigned int adjacent_mask = faces_mask1 & ~get_mask(id);
        unsigned int adjacent_id = face_mask_2_face_index(adjacent_mask);

        // Find adjacent vertex and update incidence relation:
        unsigned int adjacent_i = (adjacent_id + 0) % 3;
        CGAL_assertion(!proj_normal1.is_vertex_set(adjacent_i));
        Point_2 adjacent_p = construct_point(point1, adjacent_id);
        Arr_halfedge_handle adjacent_edge = split(adjacent_id, adjacent_p,
                                                 proj_normal1);
        const Arr_vertex_handle & adjacent_v = adjacent_edge->target();
        adjacent_v->set_location(Arr_vertex::Edge);
        if (proj_normal1.get_is_real()) adjacent_v->set_is_real(true);
        proj_normal1.set_vertex(adjacent_v, adjacent_i);
        v->set_adjacent_vertex(*((void **) (&adjacent_v)));
        adjacent_v->set_adjacent_vertex(*((void **) (&v)));
      }
    } else if (proj_normal1.get_num_faces() == 3) {
      unsigned index = get_corner_index(faces_mask1, id);
      CGAL_assertion(index < NUM_CORNERS);
      Arr_vertex_handle v = m_corner_vertices[id][index];
      /*! \todo Set the entire cyclic chain of 3 vertices, and not just the
       * single current vertex
       */
      if (proj_normal1.get_is_real()) v->set_is_real(true);
      Arr_halfedge_handle edge = v->incident_halfedges();
      if (!proj_normal1.is_vertex_set(i)) proj_normal1.set_vertex(v, i);
    }

    /* If the computed point is on the boundary of the box face, but not a
     * corner point, locate it. If the locate returns an edge, split it at the
     * point:
     */
    if (proj_normal2.get_num_faces() == 2) {
      if (!proj_normal2.is_vertex_set(i)) {
        Arr_halfedge_handle edge = split(id, p2, proj_normal2);
        const Arr_vertex_handle & v = edge->target();
        proj_normal2.set_vertex(v, i);
        v->set_face_id(id);
        v->set_location(Arr_vertex::Edge);
        if (proj_normal2.get_is_real()) v->set_is_real(true);

        // Find adjacent face:
        unsigned int adjacent_mask = faces_mask2 & ~get_mask(id);
        unsigned int adjacent_id = face_mask_2_face_index(adjacent_mask);

        // Find adjacent vertex and update incidence relation:
        unsigned int adjacent_i = (adjacent_id + 0) % 3;
        CGAL_assertion(!proj_normal2.is_vertex_set(adjacent_i));
        Point_2 adjacent_p = construct_point(point2, adjacent_id);
        Arr_halfedge_handle adjacent_edge = split(adjacent_id, adjacent_p,
                                                 proj_normal2);
        const Arr_vertex_handle & adjacent_v = adjacent_edge->target();
        adjacent_v->set_location(Arr_vertex::Edge);
        if (proj_normal2.get_is_real()) adjacent_v->set_is_real(true);
        proj_normal2.set_vertex(adjacent_v, adjacent_i);
        v->set_adjacent_vertex(*((void **) (&adjacent_v)));
        adjacent_v->set_adjacent_vertex(*((void **) (&v)));
      }
    } else if (proj_normal2.get_num_faces() == 3) {
      unsigned index = get_corner_index(faces_mask2, id);
      CGAL_assertion(index < NUM_CORNERS);
      Arr_vertex_handle v = m_corner_vertices[id][index];
      /*! \todo Set the entire cyclic chain of 3 vertices, and not just the
       * single current vertex
       */
      if (proj_normal2.get_is_real()) v->set_is_real(true);
      Arr_halfedge_handle edge = v->incident_halfedges();
      if (!proj_normal2.is_vertex_set(i)) proj_normal2.set_vertex(v, i);
    }

    /* If the points are on the same boundary, no need to insert a curve.
     * Update the is_real_flag and return
     */
    if (is_on_same_boundary(faces_mask1, faces_mask2, id)) {
      const Arr_vertex_handle & v1 = proj_normal1.get_vertex(i);
      const Arr_vertex_handle & v2 = proj_normal2.get_vertex(i);

      // Find the edge:
      Arr_halfedge_around_vertex_circulator edge = v2->incident_halfedges();
      Arr_halfedge_around_vertex_circulator begin_hec = edge;
      do {
        if (edge->source() == v1) break;
        ++edge;
      } while (edge != begin_hec);
      CGAL_assertion((v1 == edge->source()) && (v2 == edge->target()));
      
      edge->set_is_real(true);
      edge->twin()->set_is_real(true);

      handle_new_boundary_edge(edge);

      // Insert an edge into the adjacent face:
      Arr_halfedge_handle adjacent_edge;
      if (edge->face()->is_unbounded()) 
        adjacent_edge = get_adjacent_halfedge_handle(edge->twin());
      else 
        adjacent_edge = get_adjacent_halfedge_handle(edge);
      adjacent_edge->set_is_real(true);
      adjacent_edge->twin()->set_is_real(true);
      
      return;
    }

    /* Insert a curve to the arrangement that corresponds to the unit-cube face
     * with index 'id':
     */  
    X_monotone_curve_2 cv(p1, p2);
    Arr_halfedge_handle edge;
    if (proj_normal1.is_vertex_set(i) && proj_normal2.is_vertex_set(i)) {
      const Arr_vertex_handle & v1 = proj_normal1.get_vertex(i);
      const Arr_vertex_handle & v2 = proj_normal2.get_vertex(i);
      // std::cout << "cv: " << cv
      // << ", v1: " << v1->point() << ", v2: " << v2->point()
      // << std::endl;
      edge = m_arrangements[id].insert_at_vertices(cv, v1, v2);
    
      // Set the fields that indicates whehther the vertex is interior or not
      v1->set_location(num_faces_2_vertex_location(proj_normal1.get_num_faces()));
      v2->set_location(num_faces_2_vertex_location(proj_normal2.get_num_faces()));

      // Set the flags that indicates whehther the vertex is real:
      if (proj_normal1.get_is_real()) v1->set_is_real(true);
      if (proj_normal2.get_is_real()) v2->set_is_real(true);
    } else if (proj_normal1.is_vertex_set(i)) {
      const Arr_vertex_handle & v1 = proj_normal1.get_vertex(i);

      //! \todo this can done before cv is constructed!
      // std::cout << "1: " << p1 << "," << p2 << ", v1: " << v1->point()
      // << std::endl;
      Arr_traits traits;
      Comparison_result res = traits.compare_xy_2_object()(p1, p2);
      CGAL_assertion(res != EQUAL);
      if (res == SMALLER)
        edge = m_arrangements[id].insert_from_left_vertex(cv, v1);
      else
        edge = m_arrangements[id].insert_from_right_vertex(cv, v1);
      const Arr_vertex_handle & v2 = edge->target();
      proj_normal2.set_vertex(v2, i);
      v2->set_face_id(id);
      v1->set_location(num_faces_2_vertex_location(proj_normal1.get_num_faces()));
      
      // Set the flags that indicates whehther the vertex is real:
      if (proj_normal1.get_is_real()) v1->set_is_real(true);
      if (proj_normal2.get_is_real()) v2->set_is_real(true);
    } else if (proj_normal2.is_vertex_set(i)) {
      const Arr_vertex_handle & v2 = proj_normal2.get_vertex(i);

      //! \todo this can done before cv is constructed!
      // std::cout << "2: " << p1 << "," << p2 << std::endl;
      Arr_traits traits;
      Comparison_result res = traits.compare_xy_2_object()(p2, p1);
      CGAL_assertion(res != EQUAL);
      if (res == SMALLER)
        edge = m_arrangements[id].insert_from_left_vertex(cv, v2);
      else
        edge = m_arrangements[id].insert_from_right_vertex(cv, v2);
      const Arr_vertex_handle & v1 = edge->target();
      proj_normal1.set_vertex(v1, i);
      v1->set_face_id(id);      
      v2->set_location(num_faces_2_vertex_location(proj_normal2.get_num_faces()));
      // Set the flags that indicates whehther the vertex is real:
      if (proj_normal1.get_is_real()) v1->set_is_real(true);
      if (proj_normal2.get_is_real()) v2->set_is_real(true);     
    } else {
      //! \todo use a default point location
      edge = insert_non_intersecting(m_arrangements[id], cv);
      const Arr_vertex_handle & v1 = edge->source();
      const Arr_vertex_handle & v2 = edge->target();

      Arr_traits traits;
      Comparison_result res = traits.compare_xy_2_object()(p1, p2);
      CGAL_assertion(res != EQUAL);
      if (res == SMALLER) {
        proj_normal1.set_vertex(v1, i);
        proj_normal2.set_vertex(v2, i);
      } else {
        proj_normal1.set_vertex(v2, i);
        proj_normal2.set_vertex(v1, i);
      }
      v1->set_face_id(id);
      v2->set_face_id(id);

      // Set the flags that indicates whehther the vertex is real:
      if (proj_normal1.get_is_real()) v1->set_is_real(true);
      if (proj_normal2.get_is_real()) v2->set_is_real(true);
    }
    edge->set_is_real(true);
    edge->twin()->set_is_real(true);

    handle_new_edge(edge);
  }

  /*! Connect two projected normals. In particular, insert segments that connect
   * the source projected normal, which projects onto the unit-cube face with
   * a given id to, the target projected normal, which projects onto the
   * unit-cube face with another given id.
   * \param id1 the index of the unit-cube face the source normal project onto
   * \param id2 the index of the unit-cube face the target normal project onto
   * \param proj_normal1 represents the source projected normal
   * \param proj_normal2 represents the target projected normal
   * \param unique indicates whether projected normals are unique. If true,
   * then two different projected normals that are passed as parameters since
   * the last clear, must represents two different normals.
   */
  void connect(Projected_normal & proj_normal1, unsigned int id1,
               Projected_normal & proj_normal2, unsigned int id2,
               const Cgm_plane_3 & plane, bool unique = false)
  {
#if 0
    std::cout << "connect from " << id1 << " to " << id2 << std::endl;
#endif
    
    /* If the normals projects onto the same box face f, generate a curve
     * between the prjection points, and insert it into the
     * arrangement that corresponds to the box f:
     */
    if (id1 == id2) {
      insert(proj_normal1, proj_normal2, id1, unique);
      return;
    }

    unsigned int faces_mask1 = proj_normal1.get_faces_mask();
    unsigned int faces_mask2 = proj_normal2.get_faces_mask();
    unsigned int faces_mask = faces_mask1 & faces_mask2;
    if (faces_mask) {
      unsigned int id;
      for (id = 0; id < NUM_FACES; id++)
        if (faces_mask & get_mask(id)) break;
      insert(proj_normal1, proj_normal2, id, unique);
      return;
    }

    /* If the normals projects onto different box faces, we compute the plane
     * induced by the origin and the 2 points, advance to the
     * next box face, and try to connect again
     */

    // Extract the data of the first normal:
    const Point_3 & point1 = proj_normal1.get_projected_normal();
    
    unsigned int i = (id1 + 0) % 3;
    unsigned int j = (id1 + 1) % 3;
    unsigned int k = (id1 + 2) % 3;
  
    // Extract the data of the second normal:
    const Point_3 & point2 = proj_normal2.get_projected_normal();
  
    // We consider the boundary of box face with index 'id1':
    unsigned int l = (id1 < 3) ? 0 : 1;
    unsigned int m =
      ((faces_mask1 & get_mask(j+3)) && (faces_mask2 & get_mask(j))) ?
      0 : ((faces_mask1 & get_mask(j)) && (faces_mask2 & get_mask(j+3))) ?
      1 : ((point1[j] + point2[j]) < 0) ? 0 : 1;
    unsigned int n =
      ((faces_mask1 & get_mask(k+3)) && (faces_mask2 & get_mask(k))) ?
      0 : ((faces_mask1 & get_mask(k)) && (faces_mask2 & get_mask(k+3))) ?
      1 : ((point1[k] + point2[k]) < 0) ? 0 : 1;

    FT vec[3];
    vec[i] = get_extreme_coordinate(l);
    vec[j] = get_extreme_coordinate(m);
    vec[k] = get_extreme_coordinate(n);

    // Check whether the point does not escape the box. There are 2 options:
    faces_mask = 0x0;
    unsigned int id = 0;
    unsigned int num_faces = 0;
    Point_3 point;
    Point_2 origin(CGAL::ORIGIN);
    /* If the k-th coefficient is equal to 0, then the plane is prallel to 
     * the j-k plane or to the k-i plane. In either case, the projection of the
     * 2d point is at infinity. In this case we move to the 2nd option
     */
    bool flag = false;
    bool is_set = false;
    
    if (plane[k] != 0) {
      id = j + 3 * m;
      faces_mask = get_mask(id1) | get_mask(id);
      num_faces = 2;

      // Verify that the new point is different than point1:
      if ((faces_mask1 & faces_mask) != faces_mask) {
        Point_2 p_2(vec[i], vec[j]);
        Point_2 p1_2(point1[i], point1[j]);
        Point_2 p2_2(point2[i], point2[j]);

        typedef typename T_Kernel::Left_turn_2 Left_turn_2;
        T_Kernel kernel;
        Left_turn_2 leftturn = kernel.left_turn_2_object();
        flag = ((leftturn(p1_2, origin, p2_2) &&
                 leftturn(p_2, origin, p2_2) && leftturn(p1_2, origin, p_2)) ||
                (leftturn(p2_2, origin, p1_2) &&
                 leftturn(p_2, origin, p1_2) && leftturn(p2_2, origin, p_2)));
        point = plane.to_3d(p_2, k);

        // Check whether the new point is not outside the cube:
        if (flag && (((n == 0) && (point[k] >= get_extreme_coordinate(0))) ||
                     ((n == 1) && (point[k] <= get_extreme_coordinate(1)))))
        {
          is_set = true;

          // Check whether the new point is a cube corner-point:
          if (point[k] == get_extreme_coordinate(n)) {
            faces_mask |= get_mask(k + 3 * n);
            // If the new point coincides with the destination cube face:
            if ((k + 3 * n) == id2) id = id2;
            num_faces++;
          }
        }
      }
    }

    // If not set, try the other index:
    if (!is_set) {
      CGAL_assertion(plane[j] != 0);

      id = k + 3 * n;
      faces_mask = get_mask(id1) | get_mask(id);
      num_faces = 2;

      // Verify that the new point is different than point1:
      if ((faces_mask1 & faces_mask) != faces_mask) {
        Point_2 p_2(vec[k], vec[i]);
        point = plane.to_3d(p_2, j);

        // Check whether the new point is a cube corner-point:
        if (point[j] == get_extreme_coordinate(m)) {
          faces_mask |= get_mask(j + 3 * m);
          // If the new point coincides with the destination cube face:
          if ((j + 3 * m) == id2) id = id2;
          num_faces++;
        }

        // Check whether the new point is not outside the cube:
        if (((m == 0) && (point[j] >= get_extreme_coordinate(0))) ||
            ((m == 1) && (point[j] <= get_extreme_coordinate(1))))
        {
          is_set = true;
        }
      }
    }

    CGAL_assertion(is_set);
    Projected_normal proj_normal(point, faces_mask, num_faces);
    insert(proj_normal1, proj_normal, id1, unique);
    connect(proj_normal, id, proj_normal2, id2, plane, unique);
  }
  
  /*! \brief returns true if two vertices given by their masks are on the same
   * boundary of the a unit-cube face.
   * \param mask1 the mask of the first point.
   * \param mask2 the mask of the second point.
   * \param id the index of the unit-cube face.
   * \return true iff the vertices are on the same boundary.
   */
  bool is_on_same_boundary(unsigned int mask1, unsigned int mask2,
                           unsigned int id) const
  { return ((mask1 & ~get_mask(id)) & (mask2 & ~get_mask(id))); }
    
  /*! Obtains the index of a corner vertex
   * \param mask the mask of the corner vertex
   * \param id the index of the unit-cube face for which the index of the corner
   * is seeked.
   * \pre mask must have exactly 3 bits turned on that represent the 3 unit-cube
   * faces that share the corner.
   */
  unsigned int get_corner_index(unsigned int mask, unsigned int id) const
  {
    unsigned int j = (id + 1) % 3;
    unsigned int k = (id + 2) % 3;

    if (id < 3) std::swap(j,k);
    return ((mask & get_mask(j)) ? 0 : (mask & get_mask(j+3)) ? 2 : 4) +
    ((mask & get_mask(k)) ? 0 : (mask & get_mask(k+3)) ? 1 : 4);
  }

  /*! Process the corners
   * A corner point is shared by 3 vertices of 3 distinct arrangements.
   * If the vertex is real, we mark the vertex of the arrangement with the
   * lowest id as real and mark the remaining vertices as synthetic.
   * This way we guarantee that we traverse the projected normal exactly once
   * (when we process the arrangement with the lowest id during traversal).
   */
  void process_corners()
  {
    unsigned int arr_id[3], vertex_id[3];

    for (arr_id[0] = 0; arr_id[0] < NUM_FACES; arr_id[0]++) {
      for (vertex_id[0] = 0; vertex_id[0] < NUM_CORNERS; vertex_id[0]++) {
        unsigned int i;
        bool is_real = false;

        // Arr & arr = m_arrangements[arr_id[i]];
        for (i = 0; i < 3; i++) {
          const Arr_vertex_handle & v =
            m_corner_vertices[arr_id[i]][vertex_id[i]];
          if (v->get_is_real()) is_real = true;

          if (i == 2) break;
          Corner_id corner_id = get_next_corner_id(arr_id[i], vertex_id[i]);
          arr_id[i+1] = corner_id.first;
          vertex_id[i+1] = corner_id.second;
        }

        /* Set the flag that indicates that the normal is real for the vertex
         * of the arrangement with the lowest id, and reset the 2 vertices of
         * the remaining 2 respectively
         */
        unsigned int min_arr_id = arr_id[0];
        if (arr_id[1] < min_arr_id) min_arr_id = arr_id[1];
        if (arr_id[2] < min_arr_id) min_arr_id = arr_id[2];
        for (i = 0; i < 3; i++) {
          const Arr_vertex_handle & v =
            m_corner_vertices[arr_id[i]][vertex_id[i]];
          v->set_is_real(is_real && (arr_id[i] == min_arr_id));
          
          unsigned int next_i = (i + 1) % 3;
          Arr_vertex_handle next_v =
            m_corner_vertices[arr_id[next_i]][vertex_id[next_i]];
          // v->set_adjacent_vertex(*((void **) (&next_v)));
        }
      }
    }
  }

  /*! Process the boundary edges.
   * A non-corner point that lies on a boundary edge is shared by 2 vertices of
   * 2 distinct arrangements. This method updates the incidence relations between
   * the non-corner boundary vertices. Each boundary vertex contains a pointer
   * to a vertex handle on an adjacent unit-cube face. For non-corner boundary
   * vertices this forms a cyclic chain of 2 vertices. This method forms these
   * cyclic chains.
   * Recall that when a CGM is constructed incrementaly, the incident relation
   * of the non-corner boundary vertices are updated dynamically, and their
   * correct setting is maintained. When a cgm is constructed using an overlay
   * operation, the incident relation of the non-corner boundary vertices must
   * be updated by this method.
   */
  void process_edges()
  {
    unsigned int i;
    for (i = 0; i < NUM_EDGES; ++i) {
      // find the halfedge of the boundary of arr1:
      Corner_id corner_id1 = get_incident_corner_id(i, 0);
      Arr_vertex_handle v1 = get_corner_vertex_handle(corner_id1);
      Arr_halfedge_around_vertex_circulator tmp1 = v1->incident_halfedges();
      while (!(tmp1->face()->is_unbounded())) ++tmp1;
      Arr_halfedge_handle hec1 = tmp1->next();
      
      // find the halfedge of the boundary of arr2:
      Corner_id corner_id2 = get_incident_corner_id(i, 1);
      Arr_vertex_handle v2 = get_corner_vertex_handle(corner_id2);
      Arr_halfedge_around_vertex_circulator tmp2 = v2->incident_halfedges();
      while (!(tmp2->face()->is_unbounded())) ++tmp2;
      Arr_halfedge_handle hec2 = tmp2->next();

      // Gather the halfedges of the boundary of arr1 in reverse order:
      std::list<Arr_halfedge_handle> l;
      while (hec1->target()->get_location() != Arr_vertex::Corner) {
        l.push_front(hec1);
        hec1 = hec1->next();
      }

      // Update the incident vertices:
      typename std::list<Arr_halfedge_handle>::iterator hei1;
      for (hei1 = l.begin(); hei1 != l.end(); ++hei1) {
        Arr_vertex_handle v1 = (*hei1)->target();
        Arr_vertex_handle v2 = hec2->target();
        CGAL_assertion(v2->get_location() != Arr_vertex::Corner);
        CGAL_assertion(v1->get_location() != Arr_vertex::Corner);
        v1->set_adjacent_vertex(*((void **) (&v2)));
        v2->set_adjacent_vertex(*((void **) (&v1)));
        hec2 = hec2->next();
      }
    }
  }
  
  /*! Update the point of the adjacent faces of a face */
  void update_adjacent_faces(Arr_face_handle face, const Point_3 & point,
                             unsigned int id)
  {
    Arr_ccb_halfedge_circulator hec = face->outer_ccb();
    Arr_ccb_halfedge_circulator hec_begin = hec;
    do {
      Arr_face_handle adjacent_face = hec->twin()->face();
      if (!adjacent_face->is_unbounded() && !hec->is_org_arr(id))
        update_face(adjacent_face, point, id);
      hec++;
    } while(hec != hec_begin);
  }
  
  /*! \brief update the point of a vertex-less arrangement face */
  void update_face(Arr_face_iterator face, const Point_3 & point,
                   unsigned int id)
  {
    // If the point of the face is set already, return
    if (face->get_aux_is_set(id)) return;
    face->set_aux_point(id, point);

    // Traverse the neighbors
    update_adjacent_faces(face, point, id);

    // Traverse the holes:
    Arr_holes_iterator hi;
    for (hi = face->holes_begin(); hi != face->holes_end(); ++hi) {
      Arr_ccb_halfedge_circulator hec = (*hi);
      update_face(hec->face(), point, id);
    }
  }

  /*! \brief update the point of the vertex-less arrangement faces */
  void update_faces()
  {
    for (unsigned int i = 0; i < NUM_FACES; i++) {
      Arr & arr = m_arrangements[i];
      Arr_face_iterator fi = arr.faces_begin();
      // Skip the unbounded face
      for (++fi; fi != arr.faces_end(); ++fi) {
        if (fi->get_aux_is_set(0)) {
          const Point_3 & point = fi->get_aux_point(0);
          update_adjacent_faces(fi, point, 0);
        }
        if (fi->get_aux_is_set(1)) {
          const Point_3 & point = fi->get_aux_point(1);
          update_adjacent_faces(fi, point, 1);
        }
      }
    }
  }

#if 0
  /*! Return the vertex handle that is the projection of a normal and that the
   * normal projects onto the point nearest to the the given point
   * \param point a point on the unit-cube representing a normal.
   */
  Arr_vertex_handle find_nearest_real_vertex(unsigned int arr_id,
                                             const Point_2 & point)
  {
    Arr_vertex_handle vh;
    Arr & arr = m_arrangements[arr_id];
    Arr_locate_type lt;
    Arr_halfedge_handle edge = arr.locate(point);

    Arr_halfedge_handle edge_start = edge;
    do {
      vh = edge->target();
      if (vh->get_is_real()) return vh;
      edge = edge->next();
    } while (edge != edge_start);

    if (vh->get_is_real()) return vh;
  
    void * tmp = vh->get_adjacent_vertex();
    vh = *((Arr_vertex_handle *) (&tmp));

    edge = vh->incident_halfedges();
    edge_start = edge;
    do {
      vh = edge->target();
      if (vh->get_is_real()) return vh;
      edge = edge->next();
    } while (edge != edge_start);

    CGAL_assertion(0);
    return vh;
  }
#endif
  
  /*! Find the planar map on the adjacent cube face given a halfedge
   * \pre The given halfedge is on the boundary of the arrangement
   */
  Arr_face_handle find_adjacent_face(Arr_halfedge_handle heh) const
  {
    void * tmp = heh->source()->get_adjacent_vertex();
    CGAL_assertion(tmp);
    Arr_vertex_handle vh = *((Arr_vertex_handle *) (&tmp));

    // Find the first halfedge incident to the unbounded face:
    Arr_halfedge_around_vertex_circulator hec = vh->incident_halfedges();
    Arr_halfedge_around_vertex_circulator begin_hec = hec;
    do {
      if (hec->face()->is_unbounded()) break;
      ++hec;
    } while (hec != begin_hec);
    return hec->next()->twin()->face();
  }
  
  /*! \brief process the halfedges incident to a corner */
  template <class Halfedge_around_vertex_processor>
  void process_boundary_halfedges(Arr_vertex_const_handle vit,
                                  Halfedge_around_vertex_processor & processor)
  {
    Arr_vertex_const_handle vit_begin = vit;
    do {
      // Find the first halfedge after the one incident to the unbounded face:
      Arr_halfedge_around_vertex_const_circulator hec =
        vit->incident_halfedges();
      Arr_halfedge_around_vertex_const_circulator begin_hec = hec;
      do {
        if (hec->face()->is_unbounded()) break;
        ++hec;
      } while (hec != begin_hec);

      // Increment:
      begin_hec = hec;
      hec++;

      // Check whether the first halfedge is a real boundary edge:
      if (hec->get_is_real()) {
        processor(hec);
        if (processor.terminate()) break;
      }
      hec++;

      // Traverse the remaining halfedges around the vertex:
      bool done = false;
      while (hec != begin_hec) {
        processor(hec);
        done = processor.terminate();
        if (done) break;
        ++hec;
      }
      if (done) break;

      void * tmp = vit->get_adjacent_vertex();
      vit = *((Arr_vertex_handle *) (&tmp));      
    } while (vit != vit_begin);
  }
  
public:
  /*! Parameter-less Constructor */
  Cubical_gaussian_map_3()
  {
    // The m_corner_vertices are set to NULL by their default constructor
  }

  /*! Copy Constructor */
  Cubical_gaussian_map_3(const Cubical_gaussian_map_3 & gaussian_map)
  {
    // Not implemented yet!
    CGAL_assertion(0);
  }

  /*! This function initializes the arrangements with the boundary curves. There
   * are 6 arrangements, one for each face of the axis-parallel unit cube. The 3D
   * unit cube is defined by the 2 twin corners:
   *      (s_extreme[0], s_extreme[0], s_extreme[0]) and
   *      (s_extreme[1], s_extreme[1], s_extreme[1])
   * Each 2D arrangement is initialized to consist of a single axis-parallel
   * square face (beyond the unbounded face) defined by the 2 twin
   * corners:
   *      (s_extreme[0], s_extreme[0]) and
   *      (s_extreme[1], s_extreme[1])
   */
  void init_arrangements()
  {
    const Point_2 & p0 = get_corner_point(0);
    const Point_2 & p1 = get_corner_point(1);
    const Point_2 & p2 = get_corner_point(2);
    const Point_2 & p3 = get_corner_point(3);

    X_monotone_curve_2 cv0(p0, p1);
    X_monotone_curve_2 cv1(p1, p3);
    X_monotone_curve_2 cv3(p3, p2);
    X_monotone_curve_2 cv2(p2, p0);
    
    for (unsigned int i = 0; i < NUM_FACES; ++i) {
      /*! The boundary edges */
      Arr_halfedge_handle boundary_edges[NUM_CORNERS];

      //! \todo use a default point location
      boundary_edges[0] = insert_non_intersecting(m_arrangements[i], cv0);
      boundary_edges[1] =
        m_arrangements[i].insert_from_left_vertex(cv1,
                                                  boundary_edges[0]->target());
      boundary_edges[3] =
        m_arrangements[i].insert_from_right_vertex(cv3,
                                                  boundary_edges[1]->target());
      boundary_edges[2] =
        m_arrangements[i].insert_at_vertices(cv2,
                                             boundary_edges[3]->target(),
                                             boundary_edges[0]->source());

      unsigned int j;
      for (j = 0; j < NUM_CORNERS; j++) {
        Arr_halfedge_handle he = boundary_edges[j];
        he->set_is_real(false);
        he->twin()->set_is_real(false);
        m_corner_vertices[i][j] = boundary_edges[j]->source();
        m_corner_vertices[i][j]->set_location(Arr_vertex::Corner);
        m_corner_vertices[i][j]->set_face_id(i);
      }
    }

    init_corner_incidences();
  }

  /*! Initialize the corners. Typically the init_arrangements() method
   * initializes each arrangement, introducing 4 corner points and 4 boundary
   * segments in each arrangement, and also initializes the global corner-
   * vertex data structure respectively. However, if the CGM is generated in a
   * different way, as in the case of overlay computation, where
   * init_arrangements() is not utilized. The global corner-vertex structire
   * must be initialized separately.
   */
  void init_corners()
  {
    for (unsigned int i = 0; i < NUM_FACES; ++i) {
      Arr & arr = m_arrangements[i];
      Arr_vertex_iterator vit;
      for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
        if (!vit->get_location() ==  Arr_vertex::Corner) continue;
        for (unsigned int j = 0; j < NUM_CORNERS; ++j) {
          if (get_corner_point(j) == vit->point()) {
            m_corner_vertices[i][j] = vit;
            m_corner_vertices[i][j]->set_location(Arr_vertex::Corner);
            break;
          }
        }
      }
      unsigned int j;
      for (j = 0; j < NUM_CORNERS; ++j) {
        CGAL_assertion(*((void **)(&m_corner_vertices[i][j])) != 0);
      }
    }

    init_corner_incidences();
  }
  
  /*! Initialize the incidence relations between the corners. Each boundary
   * vertex contains a pointer to a vertex handle on an adjacent unit-cube face.
   * For corner vertices this forms a cyclic chain of 3 vertices. For non-corner
   * boundary vertices this forms a cyclic chain of 2 vertices. This method
   * forms the 12 corner cyclic chains. Recall that the corner vertices are
   * created once (in init_arrangements()) or obtained once (in init_corners()),
   * and they are retained through out the lifetime of the CGM (along with
   * their cyclic chains.)
   */
  void init_corner_incidences()
  {
    for (unsigned int arr_id = 0; arr_id < NUM_FACES; arr_id++) {
      for (unsigned int ver_id = 0; ver_id < NUM_CORNERS; ver_id++) {
        Arr_vertex_handle v = m_corner_vertices[arr_id][ver_id];
        Corner_id corner_id = get_next_corner_id(arr_id, ver_id);
        unsigned int next_arr_id = corner_id.first;
        unsigned int next_ver_id = corner_id.second;
        Arr_vertex_handle next_v = m_corner_vertices[next_arr_id][next_ver_id];
        v->set_adjacent_vertex(*((void **) (&next_v)));      
      }
    }
  }

  /*! Insert a great arc whose angle is less than Pi and is represented by two
   * normals into the CGM. Each normal defines an end point of the greate arc.
   * The arc is projected and divided into segments, which are inserted into
   * the corresponding arrangements.
   * \param normal1 represents the source normal
   * \param normal2 represents the target normal
   */
  void insert(const Vector_3 & normal1, const Vector_3 & normal2)
  {
    Projected_normal proj_normal1(normal1);
    proj_normal1.set_is_real(true);

    Projected_normal proj_normal2(normal2);
    proj_normal2.set_is_real(true);

    insert(proj_normal1, proj_normal2);
  }
  
  /*! Insert a projection of a great arc represented by two projected normals
   * into the CGM. Each projected normal is the projection of the normal vector
   * that defines an end point of the greate arc. The projected arc is divided
   * into segments, which are inserted into the corresponding arrangements.
   * First, find the two unit-cube faces that the normals project onto. Then,
   * insert a segment or a sequence of segments that connect the projection of
   * the 2 normals.
   * \param proj_normal1 represents the source projected normal
   * \param proj_normal2 represents the target projected normal
   * \param unique indicates whether projected normals are unique. If true,
   * then two different projected normals that are passed as parameters since
   * the last clear, must represents two different normals.
   */
  void insert(Projected_normal & proj_normal1, Projected_normal & proj_normal2,
              bool unique = false)
  {
    unsigned int faces_mask1 = proj_normal1.get_faces_mask();
    unsigned int faces_mask2 = proj_normal2.get_faces_mask();

    /* If the normal vectors project onto the same unot-cube face f,
     * generate a curve between the prjection points, and insert it into the
     * arrangement that corresponds to the box f
     */
    unsigned int faces_mask = faces_mask1 & faces_mask2;
    if (faces_mask) {
      unsigned int id;
      for (id = 0; id < NUM_FACES; id++)
        if (faces_mask & get_mask(id)) break;
      insert(proj_normal1, proj_normal2, id, unique);
      return;
    }

    // Find the match:
    unsigned int num_faces1 = proj_normal1.get_num_faces();
    unsigned int num_faces2 = proj_normal2.get_num_faces();

    /* This is an optional optimization. If the 2 points are on twin cube
     * faces, say face i and i+3, and one of them is a boundary point, say
     * contained in faces i and j, then for this point we choose face j (we
     * elliminate the twin face i, cause it is further than that the other
     * face, j in this scenario, from i+3, where the other point lies in
     */
    if (num_faces1 > 1) {
      for (unsigned int i1 = 0; i1 < NUM_FACES; ++i1) {
        unsigned int i2 = (i1+3) % NUM_FACES;
        if ((get_mask(i1) & faces_mask1) && (get_mask(i2) & faces_mask2)) {
          faces_mask1 &= ~get_mask(i1);
          num_faces1--;
        }
      }
    }
    if (num_faces2 > 1) {
      for (unsigned int i1 = 0; i1 < NUM_FACES; ++i1) {
        unsigned int i2 = (i1+3) % NUM_FACES;
        if ((get_mask(i1) & faces_mask1) && (get_mask(i2) & faces_mask2)) {
          faces_mask2 &= ~get_mask(i2);
          num_faces2--;
        }
      }
    }

    unsigned int id1, id2;

    // Find the (first) index of the box face the first normal projects onto:
    for (id1 = 0; id1 < NUM_FACES; id1++)
      if (faces_mask1 & get_mask(id1)) break;   

    // Find the (first) index of the box face the second normal projects onto:
    for (id2 = 0; id2 < NUM_FACES; id2++)
      if (faces_mask2 & get_mask(id2)) break;    

    CGAL_assertion(id1 != id2);
    
    /* If the normals projects onto two points that are on different unit-cube
     * faces, we compute the plane induced by the origin and the two points,
     * compute all the intersections of the plane with the unit-cube
     * boundaries, and connect all the intersection properly.
     */
    const Point_3 & point1 = proj_normal1.get_projected_normal();
    const Point_3 & point2 = proj_normal2.get_projected_normal();
    Cgm_plane_3 plane(CGAL::ORIGIN, point1, point2);
    connect(proj_normal1, id1, proj_normal2, id2, plane, unique);
  }

  /* Virtual Functions */

  /*! Destructor */
  virtual ~Cubical_gaussian_map_3() { clear(); }

  /*! Handle the introduction of a new boundary edge */
  virtual void handle_new_boundary_edge(Arr_halfedge_handle edge) {}

  /*! Handle the introduction of a new edge */
  virtual void handle_new_edge(Arr_halfedge_handle edge) {}
  
  /*! Clear the internal representation and auxiliary data structures */
  void clear()
  {
    for (unsigned int j = 0; j < NUM_FACES; ++j) m_arrangements[j].clear();
  }
  
  /*! returns true if the representation is empty */
  bool is_empty() const
  {
    for (unsigned int i = 0; i < NUM_FACES; ++i) {
      CGAL_assertion(m_arrangements[i].is_valid());
      if (!m_arrangements[i].is_empty()) return false;
    }
    return true;
  }

  /*! Return the degree of a real vertex
   * \param vh the vertex handle
   * \return the degree
   * \pre vh is a handle of a real vertex
   */
  unsigned int degree(Arr_vertex_const_handle vh) const
  {
    if (vh->get_location() == Arr_vertex::Interior) return vh->degree();

    unsigned int counter = 0;
    Halfedge_around_vertex_const_circulator hec(vh->incident_halfedges());
    Halfedge_around_vertex_const_circulator begin_hec = hec;
    do {
      counter++;
      ++hec;
    } while (hec != begin_hec);
    return counter;
  }

  /*! Return the handle of a halfedge that represents the same segment of a
   * unit-cube edge on the adjacent unit-cube face
   * \param he the halfedge handle that represents a segment of an edge
   */
  Arr_halfedge_handle get_adjacent_halfedge_handle(Arr_halfedge_handle he)
  {
    void * tmp = he->source()->get_adjacent_vertex();
    Arr_vertex_handle v = *((Arr_vertex_handle *) (&tmp));

    Arr_halfedge_around_vertex_circulator hec = v->incident_halfedges();
    Arr_halfedge_around_vertex_circulator begin_hec = hec;
    do {
      if (hec->twin()->face()->is_unbounded()) break;
      ++hec;
    } while (hec != begin_hec);

    CGAL_assertion(hec->twin()->face()->is_unbounded());
    return hec;
  }

  /*! Return the handle of a halfedge that represents the same segment of a
   * unit-cube edge on the adjacent unit-cube face
   * \param he the halfedge handle that represents a segment of an edge
   */
  Arr_halfedge_const_handle
  get_adjacent_halfedge_handle(Arr_halfedge_const_handle he) const
  {
    void * tmp = he->source()->get_adjacent_vertex();
    Arr_vertex_const_handle v = *((Arr_vertex_const_handle *) (&tmp));

    Arr_halfedge_around_vertex_const_circulator hec = v->incident_halfedges();
    Arr_halfedge_around_vertex_const_circulator begin_hec = hec;
    do {
      if (hec->twin()->face()->is_unbounded()) break;
      ++hec;
    } while (hec != begin_hec);

    CGAL_assertion(hec->twin()->face()->is_unbounded());
    return hec;
  }
  
  /*! Print statistics */
  void print_stat()
  {
#if 1
    unsigned int v_num = 0, h_num = 0, f_num = 0;
    unsigned int i;
    for (i = 0; i < NUM_FACES; ++i) {
      Arr & arr = m_arrangements[i];
      std::cout << "Arrangement " << i
                << ", no. vertices: " << arr.number_of_vertices()
                << ",  no. halfedges: " << arr.number_of_halfedges()
                << ",  no. faces: " << arr.number_of_faces()
                << std::endl;
      v_num += arr.number_of_vertices();
      h_num += arr.number_of_halfedges();
      f_num += arr.number_of_faces();
    }
    std::cout << "Cubical Gaussian Map"
              << ", no. vertices: " << v_num
              << ", no. halfedges: " << h_num
              << ", no. faces: " << f_num
              << std::endl;
#endif
  }
};

#if 0
/*! output operator */
template <class T_Kernel,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T>
#endif
          class T_Dcel>
inline std::ostream &
operator<<(std::ostream & os,
           const typename Cubical_gaussian_map_3<T_Kernel,T_Dcel>::Projected_normal &
           proj_normal)
{
  return os << proj_normal.get_projected_normal() << ", "
            << std::hex << proj_normal.get_faces_mask() << ", "
            << proj_normal.get_num_faces();
}
#endif

#if 0
//! \todo This doesn't work for some reasons???

/*! Exporter
 * Recall that default template arguments may not be used in function templates
 */
template <class T_Traits,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T>
#endif
          class T_Dcel>
inline
std::ostream & operator << (std::ostream & out,
                            const Cubical_gaussian_map_3<T_Traits,T_Dcel> & cgm) 
{
  typedef Cubical_gaussian_map_3<T_Traits,T_Dcel> Cgm;
  for (unsigned int i = 0; i < Cubical_gaussian_map_3::NUM_FACES; ++i) {
    const typename Cgm::Arrangement & arr = cgm.get_arrangement(i);
    out << arr;
  }
  return out;
}

/*! Importer
 * Recall that default template arguments may not be used in function templates
 */
template <class T_Traits,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T>
#endif
          class T_Dcel>
inline
std::istream & operator >> (std::istream & in,
                            Cubical_gaussian_map_3<T_Traits,T_Dcel> & cgm)
{
  typedef Cubical_gaussian_map_3<T_Traits,T_Dcel> Cgm;
  for (unsigned int i = 0; i < Cubical_gaussian_map_3::NUM_FACES; ++i) {
    const typename Cgm::Arrangement & arr = cgm.get_arrangement(i);
    arr >> in;
  }
  return in;
}
#endif

CGAL_END_NAMESPACE

#endif
