// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s): Efi Fogel         <efif@post.tau.ac.il>
//            Naama mayer       <naamamay@post.tau.ac.il>

#ifndef CGAL_ARR_POLYHEDRAL_SGM_H
#define CGAL_ARR_POLYHEDRAL_SGM_H

/*! \file
 * Polyhedral _sgm is a data dtructure that represents a 3D convex polyhedron.
 * This representation represents the 2D surface boundary of the shape.
 */

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/type_traits.hpp>

#include <CGAL/basic.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/HalfedgeDS_vector.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_2_algorithms.h>

#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_spherical_gaussian_map_3.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_traits.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_arr_dcel.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_overlay.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_initializer_visitor.h>

namespace CGAL {

/*!
 */
template <typename PolyhedralSgm,
          typename Polyhedron,
          typename Visitor =
            Arr_polyhedral_sgm_initializer_visitor<PolyhedralSgm, Polyhedron> >
class Arr_polyhedral_sgm_initializer :
  public Arr_sgm_initializer<typename PolyhedralSgm::Base>
{
private:
  // Base type:
  typedef Arr_sgm_initializer<typename PolyhedralSgm::Base> Base;
  typedef typename PolyhedralSgm::Geometry_traits_2         Geometry_traits_2;
  typedef typename Geometry_traits_2::Point_2               Point_2;
  typedef typename Geometry_traits_2::X_monotone_curve_2    X_monotone_curve_2;
  typedef typename Geometry_traits_2::Curve_2               Curve_2;  

  typedef typename Geometry_traits_2::Point_3               Point_3;
  typedef typename Geometry_traits_2::Vector_3              Vector_3;
  
  /*! */
  typedef unsigned int *                                    Coord_index_iter;
  
  // Polyhedron types:
  typedef typename Polyhedron::Vertex_const_handle
    Polyhedron_vertex_const_handle;
  typedef typename Polyhedron::Halfedge_const_handle
    Polyhedron_halfedge_const_handle;

  typedef typename Polyhedron::Vertex_iterator
    Polyhedron_vertex_iterator;
  typedef typename Polyhedron::Halfedge_iterator
    Polyhedron_halfedge_iterator;
  typedef typename Polyhedron::Facet_iterator
    Polyhedron_facet_iterator;

  typedef typename Polyhedron::Halfedge_around_vertex_circulator
    Polyhedron_halfedge_around_vertex_circulator;
  typedef boost::is_same<typename Polyhedron::Plane_3, Vector_3>
    Polyhedron_has_normal;

  /*! Transforms a (planar) facet into a normal */
  struct Normal_equation {
    template <class Facet>
    typename Facet::Plane_3 operator()(Facet & f) {
      typename Facet::Halfedge_handle h = f.halfedge();
      return CGAL::cross_product(h->next()->vertex()->point() -
                                 h->vertex()->point(),
                                 h->next()->next()->vertex()->point() -
                                 h->next()->vertex()->point());
    }
  };

  void compute_planes(Polyhedron & polyhedron, boost::true_type)
  {
    std::transform(polyhedron.facets_begin(), polyhedron.facets_end(),
                   polyhedron.planes_begin(), Normal_equation());
  }
  
  /*! Compute the equation of the undelying plane of a facet */
  struct Plane_equation {
    template <typename Facet>
    typename Facet::Plane_3 operator()(Facet & f) {
      typename Facet::Halfedge_handle h = f.halfedge();
      return typename Facet::Plane_3(h->vertex()->point(),
                                     h->next()->vertex()->point(),
                                     h->next()->next()->vertex()->point());
    }
  };

  void compute_planes(Polyhedron & polyhedron, boost::false_type)
  {
    std::transform(polyhedron.facets_begin(), polyhedron.facets_end(),
                   polyhedron.planes_begin(), Plane_equation());
  }
  
  /*! A point adder */
  template <class HDS, class PointIterator_3>
  class Point_adder {
  private:
    typedef Polyhedron_incremental_builder_3<HDS>               Builder;
    Builder & m_B;

  public:
    typedef typename Polyhedron::Vertex_handle
      Polyhedron_vertex_handle;

    /*! Constructor */
    Point_adder(Builder & B) : m_B(B) {}
      
    Polyhedron_vertex_handle operator()(PointIterator_3 pi)
    {
      typedef typename HDS::Vertex      Vertex;
      typedef typename Vertex::Point    Point;
      return m_B.add_vertex(Point((*pi)[0], (*pi)[1], (*pi)[2]));
    }
  };

  /*! Specialized point adder */
  template <class HDS> class Point_adder<HDS, Point_3 *> {
  private:
    typedef Polyhedron_incremental_builder_3<HDS>               Builder;
    Builder & m_B;

  public:
    typedef typename Polyhedron::Vertex_handle
      Polyhedron_vertex_handle;

    /*! Constructor */
    Point_adder(Builder & B) : m_B(B) {}
      
    Polyhedron_vertex_handle operator()(Point_3 * pi)
    { return m_B.add_vertex(*pi); }
  };
  
  /*! */
  template <class PointIterator_3>
  class Build_surface : public Modifier_base<typename Polyhedron::HalfedgeDS>
  {
  private:
    typedef typename Polyhedron::Vertex_handle
      Polyhedron_vertex_handle;
    typedef typename Polyhedron::Facet_handle
      Polyhedron_facet_handle;
    typedef typename Polyhedron::HalfedgeDS             HDS;
    typedef Polyhedron_incremental_builder_3<HDS>       Builder;
    typedef typename Builder::size_type                 size_type;
    typedef unsigned int *                              Coord_index_iter;

    /*! The begin iterator of the points */
    const PointIterator_3 & m_points_begin;

    /*! The end iterator of the points */
    const PointIterator_3 & m_points_end;
    
    /*! The number of points */
    size_type m_num_points;

    /*! The begin iterator of the indices */
    const Coord_index_iter & m_indices_begin;

    /*! The end iterator of the indices */
    const Coord_index_iter & m_indices_end;
    
    /*! The number of facest */
    size_type m_num_facets;

    /*! The index of the marked vertex */
    unsigned int m_marked_vertex_index;

    /*! The index of the marked edge */
    unsigned int m_marked_edge_index;

    /*! The index of the marked face */
    unsigned int m_marked_facet_index;
    
  public:
    /*! Constructor */
    Build_surface(const PointIterator_3 & points_begin,
                  const PointIterator_3 & points_end,
                  unsigned int num_points,
                  const Coord_index_iter & indices_begin,
                  const Coord_index_iter & indices_end,
                  unsigned int num_facets) :
      m_points_begin(points_begin), m_points_end(points_end),
      m_num_points(num_points),
      m_indices_begin(indices_begin), m_indices_end(indices_end),
      m_num_facets(num_facets),
      m_marked_vertex_index(0),
      m_marked_edge_index(0),
      m_marked_facet_index(0)      
    {}

    /*! Destructor */
    virtual ~Build_surface() {}

    /*! Set the marked-vertex index */
    void set_marked_vertex_index(unsigned int id) {m_marked_vertex_index = id;}

    /*! Set the marked-edge index */
    void set_marked_edge_index(unsigned int id) {m_marked_edge_index = id;}

    /*! Set the marked-face index */
    void set_marked_facet_index(unsigned int id) {m_marked_facet_index = id;}

    /*! builds the polyhedron */
    void operator()(HDS & hds)
    {
      // Postcondition: `hds' is a valid polyhedral surface.
      Builder B(hds, true);
      B.begin_surface(m_num_points, m_num_facets);
      // Add the points:
      unsigned int counter = 0;
      Point_adder<HDS, PointIterator_3> add(B);
      for (PointIterator_3 pi = m_points_begin; pi != m_points_end; ++pi) {
        Polyhedron_vertex_handle vh = add(pi);
        if (counter == m_marked_vertex_index) vh->set_marked(true);
        ++counter;
      }
      
      // Add the facets:
      bool facet_ended = true;
      counter = 0;
      for (Coord_index_iter ii = m_indices_begin; ii != m_indices_end; ++ii) {
        int index = *ii;
        if (facet_ended) {
          Polyhedron_facet_handle fh = B.begin_facet();
          if (counter == m_marked_facet_index) fh->set_marked(true);
          B.add_vertex_to_facet(index);
          facet_ended = false;
          continue;
        }
        if (index != -1) {
          B.add_vertex_to_facet(index);
          continue;
        }
        B.end_facet();
        facet_ended = true;
        ++counter;
      }
      B.end_surface();
    }
  };

  /*! A visitor class */
  Visitor * m_visitor;
  
  /*! The index of the marked vertex */
  unsigned int m_marked_vertex_index;

  /*! The index of the marked edge */
  unsigned int m_marked_edge_index;

  /*! The index of the marked face */
  unsigned int m_marked_facet_index;

  /*! */
  Polyhedron_vertex_const_handle m_src_vertex;

  /*! */
  Polyhedron_vertex_const_handle m_trg_vertex;

  /*! */
  Polyhedron_halfedge_const_handle m_halfedge;

  /*! Handle the introduction of a new edge */
  virtual void handle_new_edge(typename Base::Halfedge_handle edge)
  {
    typedef typename Base::Face_handle          Arr_face_handle;
    typedef typename Base::Vertex_handle        Arr_vertex_handle;
    
    Arr_face_handle src_face = edge->twin()->face();
    Arr_face_handle trg_face = edge->face();
    src_face->set_point(m_src_vertex->point());
    trg_face->set_point(m_trg_vertex->point());

    if (m_visitor) {
      m_visitor->update_dual_vertex(m_src_vertex, src_face);
      m_visitor->update_dual_vertex(m_trg_vertex, trg_face);

      m_visitor->update_dual_halfedge(m_halfedge, edge);
      m_visitor->update_dual_halfedge(m_halfedge, edge->twin());

//       m_visitor->update_dual_face(m_halfedge->opposite()->facet(),
//                                   edge->source());
//       m_visitor->update_dual_face(m_halfedge->facet(), edge->target());
    }
  }
  
  /*! Update the polyhedron */
  template <class PointIterator_3>  
  void update_polyhedron(Polyhedron & polyhedron,
                         const PointIterator_3 & points_begin,
                         const PointIterator_3 & points_end,
                         unsigned int num_points,
                         const Coord_index_iter indices_begin,
                         Coord_index_iter indices_end,
                         unsigned int num_facets)
  {
    /*! The builder */
    Build_surface<PointIterator_3>
      surface(points_begin, points_end, num_points,
              indices_begin, indices_end, num_facets);
    surface.set_marked_vertex_index(m_marked_vertex_index);
    surface.set_marked_edge_index(m_marked_edge_index);
    surface.set_marked_facet_index(m_marked_facet_index);
    polyhedron.delegate(surface);

    // Mark the marked (half) edges:
    unsigned int counter = 0;
    typedef typename Polyhedron::Edge_iterator Polyhedron_edge_iterator;
    Polyhedron_edge_iterator ei;
    for (ei = polyhedron.edges_begin(); ei != polyhedron.edges_end(); ++ei) {
      if (counter == m_marked_edge_index) {
        // Mark both halfedges:
        ei->set_marked(true);
        ei->opposite()->set_marked(true);
      }
      ++counter;
    }
  }

  /*! Obtain the normal of a facet of a polyhedron that supports normals */
  template <typename Facet>
  const Vector_3 & get_normal(const Facet & facet, boost::true_type) const
  { return facet->plane(); }
  
  /*! Obtain the normal of a facet of a polyhedron that supports planes */
  template <typename Facet>
  Vector_3 get_normal(const Facet & facet, boost::false_type) const
  { return facet->plane().orthogonal_vector(); }

  /*! Process a polyhedron vertex recursively constructing the Gaussian map
   * of the polyhedron
   * \param src the polyhedron vertex currently processed
   * \param first_time true if the invocation to this function is the first
   * time, and false otherwise
   */
  void process_vertex(Polyhedron_vertex_iterator src, bool first_time)
  {
    m_src_vertex = src;

    typedef typename Base::Vertex_handle                Vertex_handle;
    typedef typename Base::Halfedge_handle              Halfedge_handle;

    Vertex_handle invalid_vertex;
    
    // For each vertex, traverse incident faces:
    Polyhedron_halfedge_around_vertex_circulator hec = src->vertex_begin();

    // If the vertex is not a real vertex of the polyhedron, advance to the
    // next one:
    if (circulator_size(hec) == 0) {
      process_vertex(++src, first_time);
      return;
    }

    CGAL_assertion(circulator_size(hec) >= 3);
    Polyhedron_halfedge_around_vertex_circulator begin_hec = hec;
    Polyhedron_halfedge_around_vertex_circulator next_hec = hec;
    ++next_hec;

    /* If this is not the first invocation, advance the halfedge iterator
     * until its source vertex is processed. It is guaranteed to reach such
     * a halfedge on consecutive invocations.
     */
    if (!first_time) {
      while (!hec->opposite()->vertex()->processed()) {
        hec = next_hec;
        begin_hec = hec;
        ++next_hec;
      }
    }

    // Traverse the incident halfedges:
    do {
      if (!next_hec->processed()) {

        Vector_3 normal1 =
          get_normal(hec->facet(), Polyhedron_has_normal());
        Vector_3 normal2 =
          get_normal(next_hec->facet(), Polyhedron_has_normal());

        m_trg_vertex = next_hec->opposite()->vertex();

#if 0
        std::cout << "process_vertex trg: "
                  << static_cast<float>(todouble(m_trg_vertex->point().x()))
                  << ","
                  << static_cast<float>(todouble(m_trg_vertex->point().y()))
                  << ","
                  << static_cast<float>(todouble(m_trg_vertex->point().z()))
                  << std::endl;
#endif
        
        m_halfedge = next_hec;
#if 0
        Halfedge_handle he = this->insert(normal1, normal2);
#else
        Vertex_handle v1 = hec->facet()->vertex();
        Vertex_handle v2 = next_hec->facet()->vertex();
        /* The arc might be non-x-monotone. In this case, it is broken into 2
         * x-monotone curves. The halfedges of both are obtained.
         */
        typedef typename std::list<Halfedge_handle>     Halfedge_list;
        typedef typename Halfedge_list::iterator        Halfedge_list_iter;
        Halfedge_list hes;
        if (first_time) {
          this->insert(normal1, normal2, std::back_inserter(hes));
          Halfedge_list_iter first = hes.begin();
          Halfedge_list_iter last = hes.end();
          --last;
          hec->facet()->set_vertex((*first)->source());
          next_hec->facet()->set_vertex((*last)->target());
          first_time = false;
        } else {
          if (v1 != invalid_vertex && v2 != invalid_vertex) {
            this->insert(normal1, v1, normal2, v2, std::back_inserter(hes));
          } else if (v1 != invalid_vertex) {
            this->insert(normal1, v1, normal2, std::back_inserter(hes));
            Halfedge_list_iter last = hes.end();
            --last;
            next_hec->facet()->set_vertex((*last)->target());
          } else if (v2 != invalid_vertex) {
            this->insert(normal1, normal2, v2, std::back_inserter(hes));
            Halfedge_list_iter first = hes.begin();
            hec->facet()->set_vertex((*first)->source());
          } else CGAL_error();
        }
#endif
        next_hec->set_processed(true);
        next_hec->opposite()->set_processed(true);
        Halfedge_list_iter first = hes.begin();
        if (v1 != invalid_vertex && v2 != invalid_vertex)
          handle_new_edge(*first);

        /*! \todo use is_valid!
         * this->m_sgm.is_valid();
         */
        if (m_visitor) {
          m_visitor->update_dual_face(m_halfedge->opposite()->facet(),
                                      (*first)->source());
          // m_visitor->update_dual_face(m_halfedge->facet(), (*first)->target());
        }
      }
      hec = next_hec;
      ++next_hec;
    } while (hec != begin_hec);
    src->set_processed(true);

    // Traverse recursively:
    hec = src->vertex_begin();
    begin_hec = hec;
    do {
      if (!(hec->opposite()->vertex()->processed()))
        process_vertex(hec->opposite()->vertex(), false);
      ++hec;
    } while (hec != begin_hec);
  }
  
  /*! Compute the spherical gaussian map of a convex polyhedron
   * \param polyhedron the input polyhedron
   */
  void compute_sgm(Polyhedron & polyhedron)
  {
    typedef typename Base::Vertex_handle                Vertex_handle;

    // Clear the polyhedron:
    Polyhedron_facet_iterator fi;
    for (fi = polyhedron.facets_begin(); fi != polyhedron.facets_end(); ++fi)
      fi->set_vertex(Vertex_handle());

    Polyhedron_halfedge_iterator hei;
    for (hei = polyhedron.halfedges_begin(); hei != polyhedron.halfedges_end();
         ++hei)
      hei->set_processed(false);

    Polyhedron_vertex_iterator vi;
    for (vi = polyhedron.vertices_begin(); vi != polyhedron.vertices_end();
         ++vi)
      vi->set_processed(false);

    // Traverse all verticess recursively:
    process_vertex(polyhedron.vertices_begin(), true);
  }
  
public:
  /*! Constructor */
  Arr_polyhedral_sgm_initializer(PolyhedralSgm & sgm) :
    Base(sgm),
    m_visitor(NULL),
    m_marked_vertex_index(0),
    m_marked_edge_index(0),
    m_marked_facet_index(0)
  {}
    
  /*! Destructor */
  virtual ~Arr_polyhedral_sgm_initializer() {}

  /*! Initialize the Gaussian map
   * \param polyhedron
   * \param visitor
   * \pre The polyhedron polyhedron does not have coplanar facets.
   */
  void operator()(Polyhedron & polyhedron, Visitor * visitor = NULL)
  {
#if 0
    std::copy(polyhedron.points_begin(), polyhedron.points_end(),
              std::ostream_iterator<Point_3>(std::cout, "\n"));
#endif

    m_visitor = visitor;

#if 0
    if (!polyhedron.normalized_border_is_valid())
      polyhedron.normalize_border();
#else
    polyhedron.normalize_border();
#endif
#if 1
    compute_planes(polyhedron, Polyhedron_has_normal());
#endif

    compute_sgm(polyhedron);
  }
  
  /*! Initialize the Spherical Gaussian map */
  template <class PointIterator_3>
  void operator()(const PointIterator_3 & points_begin,
                  const PointIterator_3 & points_end,
                  unsigned int num_points,
                  const Coord_index_iter indices_begin,
                  Coord_index_iter indices_end,
                  unsigned int num_facets,
                  Visitor * visitor = NULL)
  {
    m_visitor = visitor;
 
    Polyhedron polyhedron;
    update_polyhedron(polyhedron, points_begin, points_end, num_points,
                      indices_begin, indices_end, num_facets);

#if 0
    std::copy(polyhedron.points_begin(), polyhedron.points_end(),
              std::ostream_iterator<Point_3>(std::cout, "\n"));
#endif

#if 0
    if (!polyhedron.normalized_border_is_valid())
      polyhedron.normalize_border();
#else
    polyhedron.normalize_border();
#endif
#if 1
    compute_planes(polyhedron, Polyhedron_has_normal());
#endif

    compute_sgm(polyhedron);
    polyhedron.clear();
  }

  /*! Set the marked-vertex index */
  void set_marked_vertex_index(unsigned int id) {m_marked_vertex_index = id;}

  /*! Set the marked-edge index */
  void set_marked_edge_index(unsigned int id) {m_marked_edge_index = id;}

  /*! Set the marked-face index */
  void set_marked_facet_index(unsigned int id) {m_marked_facet_index = id;}
};

/*!
 */
template <class Geometry_traits_2,
          template <class T>
          class T_Dcel = Arr_polyhedral_sgm_arr_dcel>
class Arr_polyhedral_sgm :
  public Arr_spherical_gaussian_map_3<Geometry_traits_2, T_Dcel>
{
private:
  typedef Arr_polyhedral_sgm<Geometry_traits_2, T_Dcel>     Self;
  
public:
  typedef typename Geometry_traits_2::Point_3               Point_3;
  typedef typename Geometry_traits_2::Vector_3              Vector_3;

  typedef T_Dcel<Geometry_traits_2>                         Dcel;
  
  // For some reason MSVC barfs on the friend statement below. Therefore,
  // we declare the Base to be public to overcome the problem.
  typedef Arr_spherical_gaussian_map_3<Geometry_traits_2, T_Dcel>   Base;

  // WE NEED TO ADD THE CGAL NAMESPACE TO PACIFY THE G++ 4.3.3 COMPILER.
  typedef CGAL::Arr_polyhedral_sgm_overlay<Self>
    Arr_polyhedral_sgm_overlay;

#if 0
  /*! Allow the initializer to update the SGM data members */
  template <class Polyhedron, class Visitor>
  friend class Arr_polyhedral_sgm_initializer<Self, Polyhedron, Visitor>;
#endif
  
private:
  /*! The gravity center */
  Point_3 m_center;
  
  /*! Indicated whether the center has been calculated */
  bool m_dirty_center;
  
  /*! Calculate the center of the polyhedron */
  void calculate_center()
  {
    // Count them:
    typename Base::Face_handle fi;
    for (fi = this->faces_begin(); fi != this->faces_end(); fi++) {
      const Point_3 & p = fi->point();
      Vector_3 v = p - CGAL::ORIGIN;
      m_center = m_center + v;
    }

    typename Base::Size num = this->number_of_faces();
    m_center =
      Point_3(m_center.x() / num, m_center.y() / num, m_center.z() / num);

    m_dirty_center = false;
  }
  
public:
  /*! Parameter-less Constructor */
  Arr_polyhedral_sgm() : m_dirty_center(true) {}
  
  /*! Copy Constructor */
  Arr_polyhedral_sgm(const Self & sgm)
  {
    assign(sgm);
  }

  /*! Assign a spherical Gaussian map to this */
  void assign(const Self & sgm)
  {
    // Call the assign of the base class.
    Base::assign(sgm);

    typename Dcel::Face_iterator fit;
    typename Base::Face_const_iterator fit1 = sgm.faces_begin();
	  
    // Set the points of the faces.
    for (fit = this->_dcel().faces_begin(); fit != this->_dcel().faces_end();
         ++fit)
    {
      fit->set_point (fit1->point());
      ++fit1;
    }
	  	  
    typename Dcel::Edge_iterator eit;
    typename Base::Edge_const_iterator eit1 = sgm.edges_begin();

    // Set the arr_mask of the edges
    for (eit = this->_dcel().edges_begin(); eit != this->_dcel().edges_end();
         ++eit)
    {
      eit->set_arr(eit1->arr_mask());
      eit->opposite()->set_arr(eit1->arr_mask());
      ++eit1;
    }
  }
  
  /*! Destructor */
  virtual ~Arr_polyhedral_sgm() { clear(); }

  /*! Clear the internal representation and auxiliary data structures
   */
  void clear()
  {
    m_dirty_center = true;
    Base::clear();
  }
  
 // /*! Compute the minkowski sum of a range of objects of type
 //  * Arr_polyhedral_sgm
 //  */
 // template <class SgmIterator>  
 // void minkowski_sum(SgmIterator begin, SgmIterator end)
 // {
	//typename SgmIterator::value_type * sgm1 = *begin++;
 //   typename SgmIterator::value_type * sgm2 = *begin;
 //   minkowski_sum(sgm1, sgm2);
 // }

 // /*! Compute the minkowski sum of a range of objects of type
 //  * Arr_polyhedral_sgm
 //  */
 // template <typename SgmIterator, typename OverlayTraits>  
 // void minkowski_sum(SgmIterator begin, SgmIterator end,
 //                    OverlayTraits & overlay_traits)
 // {
 //   typename SgmIterator::value_type * sgm1 = *begin++;
 //   typename SgmIterator::value_type * sgm2 = *begin;
 //   minkowski_sum(sgm1, sgm2, overlay_traits);
 // }
 
  /*! Compute the Minkowski sum of 2 objects of type Arr_polyhedral_sgm
   * \param sgm1 the first Arr_polyhedral_sgm object
   * \param sgm2 the second Arr_polyhedral_sgm object
   */
  template <class Arr_polyhedral_sgm>
  void minkowski_sum(const Arr_polyhedral_sgm & sgm1,
                     const Arr_polyhedral_sgm & sgm2)
  {
    // Compute the overlays:
    Arr_polyhedral_sgm_overlay sgm_overlay;
    CGAL::overlay(sgm1, sgm2, *this, sgm_overlay);
    // print_stat();
  }

  /*! Compute the Minkowski sum of 2 objects of type Arr_polyhedral_sgm
   * \param sgm1 the first Arr_polyhedral_sgm object
   * \param sgm2 the second Arr_polyhedral_sgm object
   */
  template <class Arr_polyhedral_sgm, typename OverlayTraits>
  void minkowski_sum(const Arr_polyhedral_sgm & sgm1,
                     const Arr_polyhedral_sgm & sgm2,
                     OverlayTraits & overlay_traits)
  { CGAL::overlay(sgm1, sgm2, *this, overlay_traits); }
  
  /*! Obtain the number of (primal) vertices */
  unsigned int number_of_vertices() const
  { return (static_cast<const Base*>(this))->number_of_faces(); }
  
  /*! Obtain the number of (primal) edges
   * \return the number of (primal) edges.
   * Edges that connect vertices of degree 2 are not counted, as they have
   * been introduced only to make non-x-monotone curves x-monotone.
   *
   */
  unsigned int number_of_edges() const
  {
    unsigned int size = 0;
    typename Base::Vertex_const_iterator vit;
    for (vit = this->vertices_begin(); vit != this->vertices_end(); ++vit)
      if (vit->degree() == 2) size++;
    return (static_cast<const Base*>(this))->number_of_edges() - size;
  }

  /*! Obtain the number of (primal) facets
   * \return the number of (primal) facets.
   * Vertices of degree 2 are not counted, as they have been introduced only
   * to make non-x-monotone curves x-monotone.
   */
  unsigned int number_of_facets() const
  {
    unsigned int size = 0;
    typename Base::Vertex_const_iterator vit;
    for (vit = this->vertices_begin(); vit != this->vertices_end(); ++vit)
      if (vit->degree() > 2) size++;
    return size;
  }
  
  /* Print the sgm vertices */
  void print_vertices()   
  {
    typename Base::Face_const_iterator vit;
    for (vit = this->faces_begin(); vit != this->faces_end(); ++vit)
      std::cout << "vertex of polyhedron = " << vit->point() << std::endl;
  }

  /*! Print statistics */
  void print_stat()
  {
    Base::print_stat();
    
    std::cout << "Polyhedron"
              << ", no. facets: " << number_of_facets()
              << ", no. edges: " << number_of_edges()
              << ", no. vertices: " << number_of_vertices()
              << std::endl;
  }
};

} //namespace CGAL

#endif
