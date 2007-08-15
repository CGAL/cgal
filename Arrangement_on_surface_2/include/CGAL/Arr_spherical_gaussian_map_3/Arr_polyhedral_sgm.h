// Copyright (c) 2006 Tel-Aviv University (Israel).
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
// Author(s)     : Efi Fogel         <efif@post.tau.ac.il>

#ifndef CGAL_ARR_POLYHEDRAL_SGM_H
#define CGAL_ARR_POLYHEDRAL_SGM_H

/*! \file
 * Polyhedral _sgm is a data dtructure that represents a 3D convex polyhedron.
 * This representation represents the 2D surface boundary of the shape.
 */

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
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_polyhedron_3.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_arr_dcel.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_overlay.h>
#include <CGAL/Arr_spherical_gaussian_map_3/Arr_polyhedral_sgm_initializer_visitor.h>

#include <string>
#include <vector>
#include <list>
#include <iostream>

CGAL_BEGIN_NAMESPACE

/*!
 */
template <class PolyhedralSgm,
          class Polyhedron = Arr_polyhedral_sgm_polyhedron_3<PolyhedralSgm>,
          class Visitor = Arr_polyhedral_sgm_initializer_visitor<PolyhedralSgm> >
class Arr_polyhedral_sgm_initializer :
  public Arr_sgm_initializer<typename PolyhedralSgm::Base>
{
private:
  // Base type:
  typedef Arr_sgm_initializer<typename PolyhedralSgm::Base>
                                                          Base;

  typedef typename PolyhedralSgm::Kernel                  Kernel;
  typedef typename Kernel::FT                             FT;
  typedef typename Kernel::Point_3                        Point_3;
  typedef typename Kernel::Vector_3                       Vector_3;

  typedef typename PolyhedralSgm::Geometry_traits_2       Geometry_traits_2;
  typedef typename Geometry_traits_2::Point_2             Point_2;
  typedef typename Geometry_traits_2::X_monotone_curve_2  X_monotone_curve_2;
  typedef typename Geometry_traits_2::Curve_2             Curve_2;  
  
  /*! */
  typedef unsigned int *                                  Coord_index_iter;
  
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

  /*! Transforms a (planar) facet into a normal */
  struct Normal_vector {
    template <class Facet>
    typename Facet::Plane_3 operator()(Facet & f) {
      typename Facet::Halfedge_handle h = f.halfedge();
      // Facet::Plane_3 is the normal vector type. We assume the
      // CGAL Kernel here and use its global functions.
#if 0
      const Point_3 & x = h->vertex()->point();
      const Point_3 & y = h->next()->vertex()->point();
      const Point_3 & z = h->next()->next()->vertex()->point();
#endif 
      Vector_3 normal =
        CGAL::cross_product(h->next()->vertex()->point() -
                            h->vertex()->point(),
                            h->next()->next()->vertex()->point() -
                            h->next()->vertex()->point());
      FT sqr_length = normal.squared_length();
      double tmp = CGAL::to_double(sqr_length);
      return normal / CGAL::sqrt(tmp);
    }
  };

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

#if 0
  /*! Handle the introduction of a new boundary edge */
  virtual void handle_new_boundary_edge(Arr_halfedge_handle edge)
  {
    if (!edge->face()->is_unbounded()) {
      edge->face()->set_point(m_trg_vertex->point());
      if (m_visitor) {
        m_visitor->update_dual_vertex(m_trg_vertex, edge->face());
        m_visitor->update_dual_halfedge(m_halfedge, edge);
      }
    } else {
      edge->twin()->face()->set_point(m_src_vertex->point());
      if (m_visitor) {
        m_visitor->update_dual_vertex(m_src_vertex, edge->twin()->face());
        m_visitor->update_dual_halfedge(m_halfedge, edge->twin());
      }
    }
  }

  /*! Handle the introduction of a new edge */
  virtual void handle_new_edge(Arr_halfedge_handle edge)
  {
    Arr_face_handle src_face = edge->twin()->face();
    Arr_face_handle trg_face = edge->face();
    src_face->set_point(m_src_vertex->point());
    trg_face->set_point(m_trg_vertex->point());

    if (m_visitor) {
      m_visitor->update_dual_vertex(m_src_vertex, src_face);
      m_visitor->update_dual_vertex(m_trg_vertex, trg_face);

      m_visitor->update_dual_halfedge(m_halfedge, edge);
      m_visitor->update_dual_halfedge(m_halfedge, edge->twin());
    }
  }  
#endif
  
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

  /*! Process a polyhedron vertex
   * \param src
   * \param first_time
   */
  void process_vertex(Polyhedron_vertex_iterator src, bool first_time)
  {
    m_src_vertex = src;

#if 0
    CGAL::To_double<typename Kernel::FT> todouble;
    std::cout << "process_vertex src: "
              << static_cast<float>(todouble(m_src_vertex->point().x()))
              << ","
              << static_cast<float>(todouble(m_src_vertex->point().y()))
              << ","
              << static_cast<float>(todouble(m_src_vertex->point().z()))
              << std::endl;
#endif
    
    typedef typename Base::Vertex_handle                Vertex_handle;
    typedef typename Base::Halfedge_handle              Halfedge_handle;

    Vertex_handle invalid_vertex;
    
    // For each vertex, traverse incident faces:
    Polyhedron_halfedge_around_vertex_circulator hec = src->vertex_begin();
    CGAL_assertion(circulator_size(hec) >= 3);
    Polyhedron_halfedge_around_vertex_circulator begin_hec = hec;
    Polyhedron_halfedge_around_vertex_circulator next_hec = hec;
    ++next_hec;

    // If it's not the first time, advance the halfedge iterator until its
    // source vertex is processed:
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
        const Vector_3 & normal1 = hec->facet()->plane();
        const Vector_3 & normal2 = next_hec->facet()->plane();
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
          } else CGAL_assertion(0);
        }
#endif
        next_hec->set_processed(true);
        next_hec->opposite()->set_processed(true);
        Halfedge_list_iter first = hes.begin();
        if (v1 != invalid_vertex && v2 != invalid_vertex) {
#if 0
          const typename Kernel::Point_3 & pps = (*first)->face()->get_point();
          const typename Kernel::Point_3 & ppt = (*first)->twin()->face()->get_point();
          std::cout << "XX pps: "
                    << static_cast<float>(todouble(pps.x()))
                    << ","
                    << static_cast<float>(todouble(pps.y()))
                    << ","
                    << static_cast<float>(todouble(pps.z()))
                    << ","
                    << ", trg: "
                    << static_cast<float>(todouble(ppt.x()))
                    << ","
                    << static_cast<float>(todouble(ppt.y()))
                    << ","
                    << static_cast<float>(todouble(ppt.z()))
                    << std::endl;
#endif     
          (*first)->face()->set_point(m_trg_vertex->point());
          (*first)->twin()->face()->set_point(m_src_vertex->point());
#if 0
          std::cout << "src: "
                    << static_cast<float>(todouble(m_src_vertex->point().x()))
                    << ","
                    << static_cast<float>(todouble(m_src_vertex->point().y()))
                    << ","
                    << static_cast<float>(todouble(m_src_vertex->point().z()))
                    << ","
                    << ", trg: "
                    << static_cast<float>(todouble(m_trg_vertex->point().x()))
                    << ","
                    << static_cast<float>(todouble(m_trg_vertex->point().y()))
                    << ","
                    << static_cast<float>(todouble(m_trg_vertex->point().z()))
                    << std::endl;
#endif
        }
        /*! \todo use is_valid!
         * this->m_sgm.is_valid();
         */
#if 0
        if (m_visitor) {
          for (unsigned int i = 0; i < 3; ++i) {
            if (dual1.is_vertex_set(i)) {
              Arr_vertex_handle & vh = dual1.get_vertex(i);
              m_visitor->update_dual_face(hec->facet(), vh);
            }
            if (dual2.is_vertex_set(i)) {
              Arr_vertex_handle & vh = dual2.get_vertex(i);
              m_visitor->update_dual_face(next_hec->facet(), vh);
            }
          }
        }
#endif
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
  
  /*! Compute the spherical gaussian map */
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

  /*! Initialize the Gaussian map */
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
    std::transform(polyhedron.facets_begin(), polyhedron.facets_end(),
                   polyhedron.planes_begin(), Normal_vector());
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
    std::transform(polyhedron.facets_begin(), polyhedron.facets_end(),
                   polyhedron.planes_begin(), Normal_vector());
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
template <class T_Kernel,
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
          template <class T>
#endif
          class T_Dcel = Arr_polyhedral_sgm_arr_dcel>
class Arr_polyhedral_sgm :
  public Arr_spherical_gaussian_map_3<T_Kernel,T_Dcel>
{
private:
  typedef Arr_polyhedral_sgm<T_Kernel, T_Dcel>              Self;
  
public:
  typedef T_Kernel                                          Kernel;
  
  // For some reason MSVC barfs on the friend statement below. Therefore,
  // we declare the Base to be public to overcome the problem.
  typedef Arr_spherical_gaussian_map_3<T_Kernel, T_Dcel>    Base;

#if 0
  /*! Allow the initializer to update the SGM data members */
  template <class Polyhedron, class Visitor>
  friend class Arr_polyhedral_sgm_initializer<Self, Polyhedron, Visitor>;
#endif
  
public:
  // Arrangement traits and types:
  typedef typename Base::Geometry_traits_2                  Geometry_traits_2;
  
#ifndef CGAL_CFG_NO_TMPL_IN_TMPL_PARAM
  typedef T_Dcel<Geometry_traits_2>                         Dcel;
#else
  typedef typename T_Dcel::template Dcel<Geometry_traits_2> Dcel;
#endif

  typedef Arr_polyhedral_sgm_overlay<Self>
    Arr_polyhedral_sgm_overlay;

  typedef typename Kernel::Point_3                          Point_3;
  typedef typename Kernel::Vector_3                         Vector_3;

private:
  /*! The gravity center */
  Point_3 m_center;
  
  /*! Indicated whether the center has been calculated */
  bool m_dirty_center;
  
  /*! Calculate the center of the polyhedron */
  void calculate_center()
  {
    // Count them:
    unsigned int vertices_num = 0;
    typename Base::Face_handle fi;
    for (fi = this->faces_begin(); fi != this->faces_end(); fi++) {
      vertices_num++;
      const Point_3 & p = fi->get_point();
      Vector_3 v = p - CGAL::ORIGIN;
      m_center = m_center + v;
    }

    typedef typename Kernel::FT FT;
    FT num((int)vertices_num);
    FT x = m_center.x() / num;
    FT y = m_center.y() / num;
    FT z = m_center.z() / num;
    m_center = Point_3(x, y, z);

    m_dirty_center = false;
  }
  
public:
  /*! Parameter-less Constructor */
  Arr_polyhedral_sgm() : m_dirty_center(true) {}
  
  /*! Copy Constructor */
  Arr_polyhedral_sgm(const Arr_polyhedral_sgm & sgm)
  {
    // Not implemented yet!
    CGAL_assertion(0);
  }
  
  /*! Destructor */
  virtual ~Arr_polyhedral_sgm() { clear(); }

  /*! \brief clears the internal representation and auxiliary data structures
   */
  void clear()
  {
    m_dirty_center = true;
    Base::clear();
  }
  
  /*! Compute the minkowski sum of a range of objects of type
   * Arr_polyhedral_sgm
   */
  template <class SgmIterator>  
  void minkowski_sum(SgmIterator begin, SgmIterator end)
  {
    Arr_polyhedral_sgm * sgm1 = *begin++;
    Arr_polyhedral_sgm * sgm2 = *begin;
    minkowski_sum(sgm1, sgm2);
  }

  /*! Compute the Minkowski sum of 2 objects of type Arr_polyhedral_sgm
   * \param sgm1 the first Arr_polyhedral_sgm object
   * \param sgm2 the second Arr_polyhedral_sgm object
   */
  void minkowski_sum(Arr_polyhedral_sgm * sgm1, Arr_polyhedral_sgm * sgm2)
  {
    // Compute the overlays:
    Arr_polyhedral_sgm_overlay sgm_overlay;
    CGAL::overlay(*sgm1, *sgm2, *this, sgm_overlay);
    // print_stat();
  }

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

CGAL_END_NAMESPACE

#endif
