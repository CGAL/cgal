// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s): Ron Wein          <wein@post.tau.ac.il>
//            Efi Fogel         <efif@post.tau.ac.il>
//            Eric Berberich    <ericb@post.tau.ac.il>
//            (based on old version by: Iddo Hanniel,
//                                      Eyal Flato,
//                                      Oren Nechushtan,
//                                      Ester Ezra,
//                                      Shai Hirsch,
//                                      and Eugene Lipovetsky)
#ifndef CGAL_ARRANGEMENT_ON_SURFACE_2_H
#define CGAL_ARRANGEMENT_ON_SURFACE_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * The header file for the Arrangement_on_surface_2<Traits,Dcel> class.
 */

#include <map>
#include <vector>
#include <algorithm>
#include <boost/mpl/assert.hpp>

#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/HalfedgeDS_iterator.h>
#include <CGAL/Arrangement_2/Arrangement_2_iterators.h>
#include <CGAL/In_place_list.h>
#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_observer.h>
#include <CGAL/Arr_accessor.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/Iterator_transform.h>

namespace CGAL {

/*! \class Arrangement_on_surface_2
 * The arrangement class, representing 2-dimensional subdivisions induced on
 * an arbitrary surface by a set of arbitrary planar curves.
 * The GeomTraits parameter corresponds to a geometry-traits class that
 * defines the Point_2 and X_monotone_curve_2 types and implements the
 * geometric predicates and constructions for the family of curves it defines.
 * The TopTraits parameter corresponds to a topology-traits class that defines
 * the topological structure of the surface. Note that the geometry traits
 * class should also be aware of the kind of surface on which its curves and
 * points are defined.
 */
template <typename GeomTraits_, typename TopTraits_>
class Arrangement_on_surface_2 {
public:
  typedef GeomTraits_                                     Geometry_traits_2;
  typedef TopTraits_                                      Topology_traits;

  // first define adaptor ...
  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2>   Traits_adaptor_2;

  // .. as it completes (potentially) missing side tags
  typedef typename Traits_adaptor_2::Left_side_category   Left_side_category;
  typedef typename Traits_adaptor_2::Bottom_side_category Bottom_side_category;
  typedef typename Traits_adaptor_2::Top_side_category    Top_side_category;
  typedef typename Traits_adaptor_2::Right_side_category  Right_side_category;

  BOOST_MPL_ASSERT(
                   (typename
                    Arr_sane_identified_tagging<Left_side_category,
                    Bottom_side_category,
                    Top_side_category,
                    Right_side_category>::result)
                   );

public:
  typedef Arrangement_on_surface_2<Geometry_traits_2, Topology_traits>
  Self;

  typedef typename Geometry_traits_2::Point_2             Point_2;
  typedef typename Geometry_traits_2::X_monotone_curve_2  X_monotone_curve_2;

  // maybe remove this in a future version (that supports complete handling
  // of all sides)
  typedef typename Arr_are_all_sides_oblivious_tag<Left_side_category,
                                                   Bottom_side_category,
                                                   Top_side_category,
                                                   Right_side_category>::result
    Are_all_sides_oblivious_category;

  typedef typename Arr_has_identified_sides<Left_side_category,
                                            Bottom_side_category>::result
    Has_identified_sides_category;

  typedef typename Arr_two_sides_category<Bottom_side_category,
                                          Top_side_category>::result
    Top_or_bottom_sides_category;

public:
  typedef typename Topology_traits::Dcel            Dcel;
  typedef typename Dcel::Size                       Size;

protected:
  friend class Arr_observer<Self>;
  friend class Arr_accessor<Self>;

  // Internal DCEL types:
  typedef typename Dcel::Vertex                     DVertex;
  typedef typename Dcel::Halfedge                   DHalfedge;
  typedef typename Dcel::Face                       DFace;
  typedef typename Dcel::Outer_ccb                  DOuter_ccb;
  typedef typename Dcel::Inner_ccb                  DInner_ccb;
  typedef typename Dcel::Isolated_vertex            DIso_vertex;

  typedef typename Dcel::difference_type            DDifference;
  typedef typename Dcel::iterator_category          DIterator_category;

  typedef typename Dcel::Vertex_iterator            DVertex_iter;
  typedef typename Dcel::Vertex_const_iterator      DVertex_const_iter;

  typedef typename Dcel::Halfedge_iterator          DHalfedge_iter;
  typedef typename Dcel::Halfedge_const_iterator    DHalfedge_const_iter;

  typedef typename Dcel::Edge_iterator              DEdge_iter;
  typedef typename Dcel::Edge_const_iterator        DEdge_const_iter;

  typedef typename Dcel::Face_iterator              DFace_iter;
  typedef typename Dcel::Face_const_iterator        DFace_const_iter;

  typedef typename DFace::Outer_ccb_iterator        DOuter_ccb_iter;
  typedef typename DFace::Outer_ccb_const_iterator  DOuter_ccb_const_iter;

  typedef typename DFace::Inner_ccb_iterator        DInner_ccb_iter;
  typedef typename DFace::Inner_ccb_const_iterator  DInner_ccb_const_iter;

  typedef typename DFace::Isolated_vertex_iterator  DIso_vertex_iter;
  typedef typename DFace::Isolated_vertex_const_iterator
                                                    DIso_vertex_const_iter;

protected:
  /*! \class
   * A functor for filtering DCEL vertices at infinity.
   */
  class _Is_concrete_vertex {
  private:
    const Topology_traits* m_topol_traits;

  public:
    _Is_concrete_vertex() : m_topol_traits(NULL) {}

    _Is_concrete_vertex(const Topology_traits* topol_traits) :
      m_topol_traits(topol_traits)
    {}

    bool operator()(const DVertex& v) const
    {
      if (m_topol_traits == NULL)
        return true;

      return (m_topol_traits->is_concrete_vertex(&v));
    }
  };

  /*! \class
   * A functor for filtering fictitious DCEL vertices.
   */
  class _Is_valid_vertex {
  private:
    const Topology_traits* m_topol_traits;

  public:
    _Is_valid_vertex() : m_topol_traits(NULL) {}

    _Is_valid_vertex(const Topology_traits* topol_traits) :
      m_topol_traits(topol_traits)
    {}

    bool operator()(const DVertex& v) const
    {
      if (m_topol_traits == NULL)
        return true;

      return (m_topol_traits->is_valid_vertex(&v));
    }
  };

  /*! \struct
   * A functor for filtering fictitious DCEL halfedges.
   */
  class _Is_valid_halfedge {
  private:
    const Topology_traits* m_topol_traits;

  public:
    _Is_valid_halfedge() : m_topol_traits(NULL) {}

    _Is_valid_halfedge(const Topology_traits* topol_traits) :
      m_topol_traits(topol_traits)
    {}

    bool operator()(const DHalfedge& he) const
    {
      if (m_topol_traits == NULL)
        return true;

      return (m_topol_traits->is_valid_halfedge(&he));
    }
  };

  /*! \struct
   * A functor for filtering the fictitious faces.
   */
  class _Is_valid_face {
  private:
    const Topology_traits* m_topol_traits;

  public:
    _Is_valid_face() : m_topol_traits(NULL) {}

    _Is_valid_face(const Topology_traits* topol_traits) :
      m_topol_traits(topol_traits)
    {}

    bool operator()(const DFace& f) const
    {
      if (m_topol_traits == NULL)
        return true;

      return (m_topol_traits->is_valid_face(&f));
    }
  };

  /*! \struct
   * A functor for filtering bounded faces.
   */
  class _Is_unbounded_face {
  private:
    const Topology_traits* m_topol_traits;

  public:
    _Is_unbounded_face() : m_topol_traits(NULL) {}

    _Is_unbounded_face(const Topology_traits* topol_traits) :
      m_topol_traits(topol_traits)
    {}

    const Topology_traits* topology_traits() const { return m_topol_traits; }

    bool operator()(const DFace& f) const
    {
      return (m_topol_traits->is_valid_face(&f) &&
              m_topol_traits->is_unbounded(&f));
    }
  };

public:
  // Forward declerations:
  class Vertex;
  class Halfedge;
  class Face;

  // Definition of the halfedge data-structure itereators and circulators:
  typedef I_Filtered_iterator<DVertex_iter, _Is_concrete_vertex,
                              Vertex, DDifference, DIterator_category>
    Vertex_iterator;

  typedef I_Filtered_const_iterator<DVertex_const_iter, _Is_concrete_vertex,
                                    DVertex_iter, Vertex, DDifference,
                                    DIterator_category>
    Vertex_const_iterator;

  typedef I_Filtered_iterator<DHalfedge_iter, _Is_valid_halfedge,
                              Halfedge, DDifference, DIterator_category>
    Halfedge_iterator;

  typedef I_Filtered_const_iterator<DHalfedge_const_iter, _Is_valid_halfedge,
                                    DHalfedge_iter, Halfedge, DDifference,
                                    DIterator_category>
    Halfedge_const_iterator;

  /*! \class
   * Edges iterator - defined as a derived class to make it assignable
   * to the halfedge iterator type.
   */
  class Edge_iterator :
    public I_Filtered_iterator<DEdge_iter, _Is_valid_halfedge,
                               Halfedge, DDifference, DIterator_category>
  {
    typedef I_Filtered_iterator<DEdge_iter, _Is_valid_halfedge,
                                Halfedge, DDifference, DIterator_category>
    Base;

  public:
    Edge_iterator() {}

    Edge_iterator(DEdge_iter iter, DEdge_iter iend,
                  const _Is_valid_halfedge& pred) :
      Base(iter, iend, pred)
    {}

    // Casting to a halfedge iterator.
    operator Halfedge_iterator() const
    {
      return (Halfedge_iterator(DHalfedge_iter(this->current_iterator())));
    }

    operator Halfedge_const_iterator() const
    {
      return (Halfedge_const_iterator
              (DHalfedge_const_iter(this->current_iterator())));
    }
  };

  class Edge_const_iterator :
    public I_Filtered_const_iterator<DEdge_const_iter, _Is_valid_halfedge,
                                     DEdge_iter, Halfedge, DDifference,
                                     DIterator_category>
  {
    typedef I_Filtered_const_iterator<DEdge_const_iter, _Is_valid_halfedge,
                                      DEdge_iter, Halfedge, DDifference,
                                      DIterator_category>    Base;

  public:
    Edge_const_iterator() {}

    Edge_const_iterator(Edge_iterator iter) :
      Base(iter.current_iterator(), iter.past_the_end(), iter.filter())
    {}

    Edge_const_iterator(DEdge_const_iter iter, DEdge_const_iter iend,
                        const _Is_valid_halfedge& pred) :
      Base(iter, iend, pred)
    {}

    // Casting to a halfedge iterator.
    operator Halfedge_const_iterator() const
    {
      return (Halfedge_const_iterator
              (DHalfedge_const_iter(this->current_iterator())));
    }
  };

  typedef I_Filtered_iterator<DFace_iter, _Is_valid_face,
                              Face, DDifference,
                              DIterator_category>     Face_iterator;

  typedef I_Filtered_const_iterator<DFace_const_iter, _Is_valid_face,
                                    DFace_iter, Face,
                                    DDifference, DIterator_category>
    Face_const_iterator;

  typedef _HalfedgeDS_vertex_circ<Halfedge, Halfedge_iterator,
                                  Bidirectional_circulator_tag>
    Halfedge_around_vertex_circulator;

  typedef _HalfedgeDS_vertex_const_circ<Halfedge, Halfedge_const_iterator,
                                        Bidirectional_circulator_tag>
    Halfedge_around_vertex_const_circulator;

  typedef _HalfedgeDS_facet_circ<Halfedge, Halfedge_iterator,
                                 Bidirectional_circulator_tag>
    Ccb_halfedge_circulator;

  typedef _HalfedgeDS_facet_const_circ<Halfedge, Halfedge_const_iterator,
                                       Bidirectional_circulator_tag>
    Ccb_halfedge_const_circulator;

  /*! \class
   * Unbounded faces iterator - defined as a derived class to make it
   * assignable to the face iterator type.
   */
  class Unbounded_face_iterator :
    public I_Filtered_iterator<DFace_iter, _Is_unbounded_face,
                               Face, DDifference, DIterator_category>
  {
    typedef I_Filtered_iterator<DFace_iter, _Is_unbounded_face,
                                Face, DDifference, DIterator_category>
      Base;

  public:
    Unbounded_face_iterator() {}

    Unbounded_face_iterator(DFace_iter iter, DFace_iter iend,
                            const _Is_unbounded_face& is_unbounded) :
      Base(iter, iend, is_unbounded)
    {}

    // Casting to a face iterator.
    operator Face_iterator() const
    {
      return (Face_iterator(DFace_iter(this->current_iterator()),
                            DFace_iter(this->past_the_end()),
                            _Is_valid_face(this->filter().topology_traits())));
    }

    operator Face_const_iterator() const
    {
      return (Face_const_iterator
              (DFace_const_iter(this->current_iterator()),
               DFace_const_iter(this->past_the_end()),
               _Is_valid_face(this->filter().topology_traits())));
    }
  };

  class Unbounded_face_const_iterator :
    public I_Filtered_const_iterator<DFace_const_iter, _Is_unbounded_face,
                                     DFace_iter, Face, DDifference,
                                     DIterator_category>
  {
    typedef I_Filtered_const_iterator<DFace_const_iter, _Is_unbounded_face,
                                      DFace_iter, Face, DDifference,
                                      DIterator_category>   Base;

  public:
    Unbounded_face_const_iterator() {}

    Unbounded_face_const_iterator(Unbounded_face_iterator iter) : Base(iter) {}

    Unbounded_face_const_iterator(DFace_const_iter iter,
                                  DFace_const_iter iend,
                                  const _Is_unbounded_face& is_unbounded) :
      Base(iter, iend, is_unbounded)
    {}

    // Casting to a face iterator.
    operator Face_const_iterator() const
    {
      return (Face_const_iterator(DFace_const_iter(this->current_iterator()),
                                  DFace_const_iter(this->past_the_end())));
    }
  };

protected:
  struct _Halfedge_to_ccb_circulator {
    typedef DHalfedge*               argument_type;
    typedef Ccb_halfedge_circulator  result_type;

    result_type operator()(argument_type s) const
    { return Ccb_halfedge_circulator(Halfedge_iterator(s)); }
  };

  struct _Const_halfedge_to_ccb_circulator {
    typedef const DHalfedge*               argument_type;
    typedef Ccb_halfedge_const_circulator  result_type;

    result_type operator()(argument_type s) const
    { return Ccb_halfedge_const_circulator(Halfedge_const_iterator(s)); }
  };

  typedef Cast_function_object<DVertex, Vertex>   _Vertex_to_vertex;

public:
  typedef Iterator_transform<DOuter_ccb_iter, _Halfedge_to_ccb_circulator>
                                      Outer_ccb_iterator;

  typedef Iterator_transform<DOuter_ccb_const_iter,
                             _Const_halfedge_to_ccb_circulator>
                                      Outer_ccb_const_iterator;

  typedef Iterator_transform<DInner_ccb_iter, _Halfedge_to_ccb_circulator>
                                      Inner_ccb_iterator;

  typedef Iterator_transform<DInner_ccb_const_iter,
                             _Const_halfedge_to_ccb_circulator>
                                      Inner_ccb_const_iterator;

  /*! \class
   * Isolated vertices iterator - defined as a class to make it assignable
   * to the vertex iterator type.
   */
  class Isolated_vertex_iterator :
    public Iterator_project<DIso_vertex_iter, _Vertex_to_vertex>
  {
    typedef Iterator_project<DIso_vertex_iter, _Vertex_to_vertex>       Base;

  public:
    Isolated_vertex_iterator() {}

    Isolated_vertex_iterator(DIso_vertex_iter iter) : Base(iter) {}

    // Casting to a vertex iterator.
    operator Vertex_iterator() const
    { return (Vertex_iterator(DVertex_iter(this->ptr()))); }

    operator Vertex_const_iterator() const
    { return (Vertex_const_iterator(DVertex_const_iter(this->ptr()))); }
  };

  class Isolated_vertex_const_iterator :
    public Iterator_project<DIso_vertex_const_iter, _Vertex_to_vertex>
  {
    typedef Iterator_project<DIso_vertex_const_iter, _Vertex_to_vertex>  Base;

  public:
    Isolated_vertex_const_iterator() {}

    Isolated_vertex_const_iterator(Isolated_vertex_iterator iter) :
      Base(iter)
    {}

    Isolated_vertex_const_iterator(DIso_vertex_const_iter iter) :
      Base(iter)
    {}

    // Casting to a vertex iterator.
    operator Vertex_const_iterator() const
    { return (Vertex_const_iterator(DVertex_const_iter(this->ptr()))); }
  };

protected:
  class _Valid_vertex_iterator :
    public I_Filtered_iterator<DVertex_iter, _Is_valid_vertex, Vertex,
                               DDifference, DIterator_category>
  {
    typedef I_Filtered_iterator<DVertex_iter, _Is_valid_vertex, Vertex,
                                DDifference, DIterator_category> Base;

  public:
    _Valid_vertex_iterator() {}

    _Valid_vertex_iterator(DVertex_iter iter, DVertex_iter iend,
                           const _Is_valid_vertex& pred) :
      Base(iter, iend, pred)
    {}

    // Casting to a vertex iterator.
    operator Vertex_iterator() const
    { return (Vertex_iterator(DVertex_iter(this->current_iterator()))); }

    operator Vertex_const_iterator() const
    {
      return (Vertex_const_iterator(DVertex_const_iter
                                    (this->current_iterator())));
    }
  };

public:
  // Definition of handles (equivalent to iterators):
  typedef Vertex_iterator              Vertex_handle;
  typedef Halfedge_iterator            Halfedge_handle;
  typedef Face_iterator                Face_handle;

  typedef Vertex_const_iterator        Vertex_const_handle;
  typedef Halfedge_const_iterator      Halfedge_const_handle;
  typedef Face_const_iterator          Face_const_handle;

  /*! \class
   * The arrangement vertex class.
   */
  class Vertex : public DVertex {
    typedef DVertex                     Base;

  public:
    /*! Default constrcutor. */
    Vertex() {}

    /*! Check whether the vertex lies on an open boundary. */
    bool is_at_open_boundary() const { return (Base::has_null_point()); }

    /*! Get the vertex degree (number of incident edges). */
    Size degree() const
    {
      if (this->is_isolated())
        return (0);

      // Go around the vertex and count the incident halfedges.
      const DHalfedge* he_first = Base::halfedge();
      const DHalfedge* he_curr = he_first;
      Size n = 0;

      if (he_curr != NULL) {
        do {
          ++n;
          he_curr = he_curr->next()->opposite();
        } while (he_curr != he_first);
      }
      return (n);
    }

    /*!
     * Get the incident halfedges (non-const version).
     * \pre The vertex is not isolated.
     */
    Halfedge_around_vertex_circulator incident_halfedges()
    {
      CGAL_precondition(! this->is_isolated());
      return Halfedge_around_vertex_circulator
        (DHalfedge_iter(Base::halfedge()));
    }

    /*!
     * Get the incident halfedges (const version).
     * \pre The vertex is not isolated.
     */
    Halfedge_around_vertex_const_circulator incident_halfedges() const
    {
      CGAL_precondition(! this->is_isolated());
      return Halfedge_around_vertex_const_circulator
        (DHalfedge_const_iter(Base::halfedge()));
    }

    /*!
     * Get the face that contains the vertex (non-const version).
     * \pre The vertex is isolated.
     */
    Face_handle face()
    {
      CGAL_precondition(this->is_isolated());
      return (DFace_iter(Base::isolated_vertex()->face()));
    }

    /*!
     * Get the face that contains the vertex (const version).
     * \pre The vertex is isolated.
     */
    Face_const_handle face() const
    {
      CGAL_precondition(this->is_isolated());
      return (DFace_const_iter(Base::isolated_vertex()->face()));
    }


  private:
    // Blocking access to inherited functions from the Dcel::Vertex.
    bool has_null_point() const;
    void set_point(Point_2* );
    void set_boundary(Arr_parameter_space , Arr_parameter_space );
    const DHalfedge* halfedge() const;
    DHalfedge* halfedge();
    void set_halfedge(DHalfedge* );
    const DIso_vertex* isolated_vertex() const;
    DIso_vertex* isolated_vertex();
    void set_isolated_vertex(DIso_vertex* );
  };

  /*!
   * \class The arrangement halfedge class.
   */
  class Halfedge : public DHalfedge {
    typedef DHalfedge             Base;

  public:
    /*! Default constrcutor. */
    Halfedge() {}

    /*! Check whether the halfedge is fictitious. */
    bool is_fictitious() const
    { return (Base::has_null_curve()); }

    /*! Get the source vertex (non-const version). */
    Vertex_handle source()
    { return (DVertex_iter(Base::opposite()->vertex())); }

    /*! Get the source vertex (const version). */
    Vertex_const_handle source() const
    { return (DVertex_const_iter(Base::opposite()->vertex())); }

    /*! Get the target vertex (non-const version). */
    Vertex_handle target()
    { return (DVertex_iter(Base::vertex())); }

    /*! Get the target vertex (const version). */
    Vertex_const_handle target() const
    { return (DVertex_const_iter(Base::vertex())); }

    /*! Get the incident face (non-const version). */
    Face_handle face()
    {
      return (! Base::is_on_inner_ccb()) ?
        DFace_iter(Base::outer_ccb()->face()) :
        DFace_iter(Base::inner_ccb()->face());
    }

    /*! Get the incident face (const version). */
    Face_const_handle face() const
    {
      return (! Base::is_on_inner_ccb()) ?
        DFace_const_iter(Base::outer_ccb()->face()) :
        DFace_const_iter(Base::inner_ccb()->face());
    }

    /*! Get the twin halfedge (non-const version). */
    Halfedge_handle twin()
    { return (DHalfedge_iter(Base::opposite())); }

    /*! Get the twin halfedge (const version). */
    Halfedge_const_handle twin() const
    { return (DHalfedge_const_iter(Base::opposite())); }

    /*! Get the previous halfegde in the chain (non-const version). */
    Halfedge_handle prev()
    { return (DHalfedge_iter(Base::prev())); }

    /*! Get the previous halfegde in the chain (const version). */
    Halfedge_const_handle prev() const
    { return (DHalfedge_const_iter(Base::prev())); }

    /*! Get the next halfegde in the chain (non-const version). */
    Halfedge_handle next()
    { return (DHalfedge_iter(Base::next())); }

    /*! Get the next halfegde in the chain (const version). */
    Halfedge_const_handle next() const
    { return (DHalfedge_const_iter(Base::next())); }

    /*! Get the connected component of the halfedge (non-const version). */
    Ccb_halfedge_circulator ccb()
    { return Ccb_halfedge_circulator(DHalfedge_iter(this)); }

    /*! Get the connected component of the halfedge (const version). */
    Ccb_halfedge_const_circulator ccb() const
    { return Ccb_halfedge_const_circulator(DHalfedge_const_iter(this)); }

  private:

    // Blocking access to inherited functions from the Dcel::Halfedge.
    bool has_null_curve() const;
    void set_curve(X_monotone_curve_2* );
    const DHalfedge* opposite() const;
    DHalfedge* opposite();
    void set_opposite(DHalfedge* );
    void set_direction(Arr_halfedge_direction );
    void set_prev(DHalfedge* );
    void set_next(DHalfedge* );
    const DVertex* vertex() const ;
    DVertex* vertex();
    void set_vertex(DVertex* );
    const DOuter_ccb* outer_ccb() const;
    DOuter_ccb* outer_ccb();
    void set_outer_ccb(DOuter_ccb* );
    const DInner_ccb* inner_ccb() const;
    DInner_ccb* inner_ccb();
    void set_inner_ccb(DInner_ccb* );
  };

  /*!
   * \class The arrangement face class.
   */
  class Face : public DFace {
    typedef DFace                 Base;

  public:
    /*! Default constrcutor. */
    Face() {}

    /*! Get an iterator for the outer CCBs of the face (non-const version). */
    Outer_ccb_iterator outer_ccbs_begin()
    { return (DOuter_ccb_iter(Base::outer_ccbs_begin())); }

    /*! Get an iterator for the outer CCBs the face (const version). */
    Outer_ccb_const_iterator outer_ccbs_begin() const
    { return (DOuter_ccb_const_iter(Base::outer_ccbs_begin())); }

    /*! Get a past-the-end iterator for the outer CCBs (non-const version). */
    Outer_ccb_iterator outer_ccbs_end()
    { return (DOuter_ccb_iter(Base::outer_ccbs_end())); }

    /*! Get a past-the-end iterator for the outer CCBs (const version). */
    Outer_ccb_const_iterator outer_ccbs_end() const
    { return (DOuter_ccb_const_iter(Base::outer_ccbs_end())); }

    /*! Get an iterator for the inner CCBs of the face (non-const version). */
    Inner_ccb_iterator inner_ccbs_begin()
    { return (DInner_ccb_iter(Base::inner_ccbs_begin())); }

    /*! Get an iterator for the inner CCBs the face (const version). */
    Inner_ccb_const_iterator inner_ccbs_begin() const
    { return (DInner_ccb_const_iter(Base::inner_ccbs_begin())); }

    /*! Get a past-the-end iterator for the inner CCBs (non-const version). */
    Inner_ccb_iterator inner_ccbs_end()
    { return (DInner_ccb_iter(Base::inner_ccbs_end())); }

    /*! Get a past-the-end iterator for the inner CCBs (const version). */
    Inner_ccb_const_iterator inner_ccbs_end() const
    { return (DInner_ccb_const_iter(Base::inner_ccbs_end())); }

    /*! Get an iterator for the isolated_vertices inside the face
     * (non-const version).
     */
    Isolated_vertex_iterator isolated_vertices_begin()
    { return (DIso_vertex_iter(Base::isolated_vertices_begin())); }

    /*! Get an iterator for the isolated_vertices inside the face
     * (const version).
     */
    Isolated_vertex_const_iterator isolated_vertices_begin() const
    { return (DIso_vertex_const_iter(Base::isolated_vertices_begin())); }

    /*! Get a past-the-end iterator for the isolated_vertices
     * (non-const version).
     */
    Isolated_vertex_iterator isolated_vertices_end()
    { return (DIso_vertex_iter(Base::isolated_vertices_end())); }

    /*! Get a past-the-end iterator for the isolated_vertices
     * (const version).
     */
    Isolated_vertex_const_iterator isolated_vertices_end() const
    { return (DIso_vertex_const_iter(Base::isolated_vertices_end())); }

    /// \name These functions are kept for Arrangement_2 compatibility:
    //@{

    /*!
     * Check whether the face has an outer CCB.
     */
    bool has_outer_ccb() const
    { return (Base::number_of_outer_ccbs() > 0); }

    /*!
     * Get a circulator for the outer boundary (non-const version).
     * \pre The face has a single outer CCB.
     */
    Ccb_halfedge_circulator outer_ccb()
    {
      CGAL_precondition(Base::number_of_outer_ccbs() == 1);

      DOuter_ccb_iter iter = Base::outer_ccbs_begin();
      DHalfedge* he = *iter;
      return Ccb_halfedge_circulator(DHalfedge_iter(he));
    }

    /*!
     * Get a circulator for the outer boundary (const version).
     * \pre The face has a single outer CCB.
     */
    Ccb_halfedge_const_circulator outer_ccb() const
    {
      CGAL_precondition(Base::number_of_outer_ccbs() == 1);

      DOuter_ccb_const_iter iter = Base::outer_ccbs_begin();
      const DHalfedge* he = *iter;
      return Ccb_halfedge_const_circulator(DHalfedge_const_iter(he));
    }

    /*! Get the number of holes (inner CCBs) inside the face. */
    Size number_of_holes() const
    { return (Base::number_of_inner_ccbs()); }

    /*! Get an iterator for the holes inside the face (non-const version). */
    Inner_ccb_iterator holes_begin()
    { return (this->inner_ccbs_begin()); }

    /*! Get an iterator for the holes inside the face (const version). */
    Inner_ccb_const_iterator holes_begin() const
    { return (this->inner_ccbs_begin()); }

    /*! Get a past-the-end iterator for the holes (non-const version). */
    Inner_ccb_iterator holes_end()
    { return (this->inner_ccbs_end()); }

    /*! Get a past-the-end iterator for the holes (const version). */
    Inner_ccb_const_iterator holes_end() const
    { return (this->inner_ccbs_end()); }
    //@}

  private:
    // Blocking access to inherited functions from the Dcel::Face.
    void set_unbounded(bool);
    void set_fictitious(bool);
    void add_outer_ccb(DOuter_ccb*, Halfedge*);
    void erase_outer_ccb(DOuter_ccb*);
    void add_inner_ccb(DInner_ccb*, Halfedge*);
    void erase_inner_ccb(DInner_ccb*);
    void add_isolated_vertex(DIso_vertex*, DVertex*);
    void erase_isolated_vertex(DIso_vertex*);
  };

protected:
  typedef CGAL_ALLOCATOR(Point_2)                 Points_alloc;
  typedef CGAL_ALLOCATOR(X_monotone_curve_2)      Curves_alloc;

  typedef Arr_observer<Self>                      Observer;
  typedef std::list<Observer*>                    Observers_container;
  typedef typename Observers_container::iterator  Observers_iterator;

  typedef typename Observers_container::reverse_iterator
                                                  Observers_rev_iterator;

  // Data members:
  Topology_traits         m_topol_traits;  // the topology traits.
  Points_alloc            m_points_alloc;  // allocator for the points.
  Curves_alloc            m_curves_alloc;  // allocator for the curves.
  Observers_container     m_observers;     // pointers to existing observers.
  const Traits_adaptor_2* m_geom_traits;   // the geometry-traits adaptor.
  bool                    m_own_traits;    // inidicates whether the geometry
                                           // traits should be freed up.

public:
  /// \name Constructors.
  //@{

  /*! Default constructor. */
  Arrangement_on_surface_2();

  /*! Copy constructor. */
  Arrangement_on_surface_2(const Self & arr);

  /*! Constructor given a traits object. */
  Arrangement_on_surface_2(const Geometry_traits_2* geom_traits);
  //@}

  /// \name Assignment functions.
  //@{

  /*! Assignment operator. */
  Self& operator=(const Self& arr);

  /*! Assign an arrangement. */
  void assign(const Self& arr);
  //@}

  /// \name Destruction functions.
  //@{

  /*! Destructor. */
  virtual ~Arrangement_on_surface_2();

  /*! Clear the arrangement. */
  virtual void clear();
  //@}

  /// \name Access the traits-class objects.
  //@{

  /*! Access the geometry-traits object (const version). */
  inline const Traits_adaptor_2* traits_adaptor() const
  { return (m_geom_traits); }

  /*! Access the geometry-traits object (const version). */
  inline const Geometry_traits_2* geometry_traits() const
  { return (m_geom_traits); }

  /*! Access the topology-traits object (non-const version). */
  inline Topology_traits* topology_traits()
  { return (&m_topol_traits); }

  /*! Access the topology-traits object (const version). */
  inline const Topology_traits* topology_traits() const
  { return (&m_topol_traits); }
  //@}

  /// \name Access the arrangement dimensions.
  //@{

  /*! Check whether the arrangement is empty. */
  bool is_empty() const
  { return (m_topol_traits.is_empty_dcel()); }

  /*!
   * Check whether the arrangement is valid. In particular, check the
   * validity of each vertex, halfedge and face, their incidence relations
   * and the geometric properties of the arrangement.
   */
  bool is_valid() const;

  /*! Get the number of arrangement vertices. */
  Size number_of_vertices() const
  { return (m_topol_traits.number_of_concrete_vertices()); }

  /*! Get the number of isolated arrangement vertices. */
  Size number_of_isolated_vertices() const
  { return (_dcel().size_of_isolated_vertices()); }

  /*! Get the number of arrangement halfedges (the result is always even). */
  Size number_of_halfedges() const
  { return (m_topol_traits.number_of_valid_halfedges()); }

  /*! Get the number of arrangement edges. */
  Size number_of_edges() const
  { return (m_topol_traits.number_of_valid_halfedges() / 2); }

  /*! Get the number of arrangement faces. */
  Size number_of_faces() const
  { return (m_topol_traits.number_of_valid_faces()); }

  /*! Get the number of unbounded faces in the arrangement. */
  Size number_of_unbounded_faces() const
  {
    Unbounded_face_const_iterator iter = unbounded_faces_begin();
    Unbounded_face_const_iterator end = unbounded_faces_end();
    Size n_unb = 0;

    while (iter != end) {
      ++n_unb;
      ++iter;
    }

    return (n_unb);
  }
  //@}

  /// \name Traversal functions for the arrangement vertices.
  //@{

  /*! Get an iterator for the first vertex in the arrangement. */
  Vertex_iterator vertices_begin()
  {
    return (Vertex_iterator(_dcel().vertices_begin(), _dcel().vertices_end(),
                            _Is_concrete_vertex(&m_topol_traits)));
  }

  /*! Get a past-the-end iterator for the arrangement vertices. */
  Vertex_iterator vertices_end()
  {
    return (Vertex_iterator(_dcel().vertices_end(), _dcel().vertices_end(),
                            _Is_concrete_vertex(&m_topol_traits)));
  }

  /*!
  returns a range over handles of the arrangement vertices .
  */
  Iterator_range<Prevent_deref<Vertex_iterator> >
  vertex_handles()
  {
    return make_prevent_deref_range(vertices_begin(), vertices_end());
  }

  /*! Get a const iterator for the first vertex in the arrangement. */
  Vertex_const_iterator vertices_begin() const
  {
    return (Vertex_const_iterator(_dcel().vertices_begin(),
                                  _dcel().vertices_end(),
                                  _Is_concrete_vertex(&m_topol_traits)));
  }

  /*! Get a past-the-end const iterator for the arrangement vertices. */
  Vertex_const_iterator vertices_end() const
  {
    return (Vertex_const_iterator(_dcel().vertices_end(),
                                  _dcel().vertices_end(),
                                  _Is_concrete_vertex(&m_topol_traits)));
  }

  /*!
  returns a const range (model of `ConstRange`) over handles of the arrangement vertices .
  */
  Iterator_range<Prevent_deref<Vertex_iterator> >
  vertex_handles() const
  {
    return make_prevent_deref_range(vertices_begin(), vertices_end());
  }

  //@}

  /// \name Traversal functions for the arrangement halfedges.
  //@{

  /*! Get an iterator for the first halfedge in the arrangement. */
  Halfedge_iterator halfedges_begin()
  {
    return (Halfedge_iterator(_dcel().halfedges_begin(),
                              _dcel().halfedges_end(),
                              _Is_valid_halfedge(&m_topol_traits)));
  }

  /*! Get a past-the-end iterator for the arrangement halfedges. */
  Halfedge_iterator halfedges_end()
  {
    return (Halfedge_iterator(_dcel().halfedges_end(),
                              _dcel().halfedges_end(),
                              _Is_valid_halfedge(&m_topol_traits)));
  }

  /*!
  returns a range over handles of the arrangement halfedges .
  */
  Iterator_range<Prevent_deref<Halfedge_iterator> >
  halfedge_handles()
  {
    return make_prevent_deref_range(halfedges_begin(), halfedges_end());
  }

  /*! Get a const iterator for the first halfedge in the arrangement. */
  Halfedge_const_iterator halfedges_begin() const
  {
    return (Halfedge_const_iterator(_dcel().halfedges_begin(),
                                    _dcel().halfedges_end(),
                                    _Is_valid_halfedge(&m_topol_traits)));
  }

  /*! Get a past-the-end const iterator for the arrangement halfedges. */
  Halfedge_const_iterator halfedges_end() const
  {
    return (Halfedge_const_iterator(_dcel().halfedges_end(),
                                    _dcel().halfedges_end(),
                                    _Is_valid_halfedge(&m_topol_traits)));
  }
  /*!
  returns a const range (model of `ConstRange`) over handles of the arrangement halfedges .
  */
  Iterator_range<Prevent_deref<Halfedge_iterator> >
  halfedge_handles() const
  {
    return make_prevent_deref_range(halfedges_begin(), halfedges_end());
  }
  //@}

  /// \name Traversal functions for the arrangement edges.
  //@{

  /*! Get an iterator for the first edge in the arrangement. */
  Edge_iterator edges_begin()
  {
    return (Edge_iterator(_dcel().edges_begin(), _dcel().edges_end(),
                          _Is_valid_halfedge(&m_topol_traits)));
  }

  /*! Get a past-the-end iterator for the arrangement edges. */
  Edge_iterator edges_end()
  {
    return (Edge_iterator(_dcel().edges_end(), _dcel().edges_end(),
                          _Is_valid_halfedge(&m_topol_traits)));
  }

  /*!
  returns a range over handles of the arrangement edges .
  */
  Iterator_range<Prevent_deref<Edge_iterator> >
  edge_handles()
  {
    return make_prevent_deref_range(edges_begin(), edges_end());
  }

  /*! Get a const iterator for the first edge in the arrangement. */
  Edge_const_iterator edges_begin() const
  {
    return (Edge_const_iterator(_dcel().edges_begin(), _dcel().edges_end(),
                                _Is_valid_halfedge(&m_topol_traits)));
  }

  /*! Get a past-the-end const iterator for the arrangement edges. */
  Edge_const_iterator edges_end() const
  {
    return (Edge_const_iterator(_dcel().edges_end(), _dcel().edges_end(),
                                _Is_valid_halfedge(&m_topol_traits)));
  }

  /*!
  returns a const range (model of `ConstRange`) over handles of the arrangement edges .
  */
  Iterator_range<Prevent_deref<Edge_iterator> >
  edge_handles() const
  {
    return make_prevent_deref_range(edges_begin(), edges_end());
  }
  //@}

  /// \name Traversal functions for the arrangement faces.
  //@{

  /*! Get an iterator for the first face in the arrangement. */
  Face_iterator faces_begin()
  {
    return (Face_iterator(_dcel().faces_begin(), _dcel().faces_end(),
                          _Is_valid_face(&m_topol_traits)));
  }

  /*! Get a past-the-end iterator for the arrangement faces. */
  Face_iterator faces_end()
  {
    return (Face_iterator(_dcel().faces_end(), _dcel().faces_end(),
                          _Is_valid_face(&m_topol_traits)));
  }

  /*!
  returns a range over handles of the arrangement faces .
  */
  Iterator_range<Prevent_deref<Face_iterator> >
  face_handles()
  {
    return make_prevent_deref_range(faces_begin(), faces_end());
  }
  /*! Get a const iterator for the first face in the arrangement. */
  Face_const_iterator faces_begin() const
  {
    return (Face_const_iterator(_dcel().faces_begin(), _dcel().faces_end(),
                                _Is_valid_face(&m_topol_traits)));
  }

  /*! Get a past-the-end const iterator for the arrangement faces. */
  Face_const_iterator faces_end() const
  {
    return (Face_const_iterator(_dcel().faces_end(), _dcel().faces_end(),
                                _Is_valid_face(&m_topol_traits)));
  }

  /*!
  returns a const range (model of `ConstRange`) over handles of the arrangement faces .
  */
  Iterator_range<Prevent_deref<Face_iterator> >
  face_handles() const
  {
    return make_prevent_deref_range(faces_begin(), faces_end());
  }
  //! reference_face (const version).
  /*! The function returns a reference face of the arrangement.
   * All reference faces of arrangements of the same type have a common
   * point.
   * \return A const handle to the reference face.
   */
  Face_const_handle reference_face() const
  {
    return _const_handle_for(this->topology_traits()->reference_face());
  }

  //! reference_face (non-const version).
  /*! The function returns a reference face of the arrangement.
    All reference faces of arrangements of the same type have a common
    point.
    \return A handle to the reference face.
  */
  Face_handle reference_face()
  { return _handle_for(this->topology_traits()->reference_face()); }

  //@}

  /// \name Traversal functions for the unbounded faces of the arrangement.
  //@{

  /*! Get an iterator for the first unbounded face in the arrangement. */
  Unbounded_face_iterator unbounded_faces_begin()
  {
    return Unbounded_face_iterator(_dcel().faces_begin(), _dcel().faces_end(),
                                   _Is_unbounded_face(&m_topol_traits));
  }

  /*! Get a past-the-end iterator for the unbounded arrangement faces. */
  Unbounded_face_iterator unbounded_faces_end()
  {
    return Unbounded_face_iterator(_dcel().faces_end(), _dcel().faces_end(),
                                   _Is_unbounded_face(&m_topol_traits));
  }

  /*! Get a const iterator for the first unbounded face in the arrangement. */
  Unbounded_face_const_iterator unbounded_faces_begin() const
  {
    return Unbounded_face_const_iterator(_dcel().faces_begin(),
                                         _dcel().faces_end(),
                                         _Is_unbounded_face(&m_topol_traits));
  }

  /*! Get a past-the-end const iterator for the unbounded arrangement faces. */
  Unbounded_face_const_iterator unbounded_faces_end() const
  {
    return Unbounded_face_const_iterator(_dcel().faces_end(),
                                         _dcel().faces_end(),
                                         _Is_unbounded_face(&m_topol_traits));
  }

  /*! Get the fictitious face (non-const version). */
  Face_handle fictitious_face()
  {
    // The fictitious contains all other faces in a single hole inside it.
    return
      Face_handle(const_cast<DFace*>(this->topology_traits()->initial_face()));
  }

  /*!
   * Get the unbounded face (const version).
   * The fictitious contains all other faces in a single hole inside it.
   */
  Face_const_handle fictitious_face() const
  { return DFace_const_iter(this->topology_traits()->initial_face()); }
  //@}

  /// \name Casting away constness for handle types.
  //@{
  Vertex_handle non_const_handle(Vertex_const_handle vh)
  {
    DVertex* p_v = (DVertex*)&(*vh);
    return (Vertex_handle(p_v));
  }

  Halfedge_handle non_const_handle(Halfedge_const_handle hh)
  {
    DHalfedge* p_he = (DHalfedge*)&(*hh);
    return (Halfedge_handle(p_he));
  }

  Face_handle non_const_handle(Face_const_handle fh)
  {
    DFace* p_f = (DFace*) &(*fh);
    return (Face_handle(p_f));
  }
  //@}

  /// \name Specilaized insertion functions.
  //@{

  /*!
   * Insert a point that forms an isolated vertex in the interior of a given
   * face.
   * \param p The given point.
   * \param f The face into which we insert the new isolated vertex.
   * \return A handle for the isolated vertex that has been created.
   */
  Vertex_handle insert_in_face_interior(const Point_2& p, Face_handle f);

  /*!
   * Insert an x-monotone curve into the arrangement as a new hole (inner
   * component) inside the given face.
   * \param cv The given x-monotone curve.
   * \param f The face into which we insert the new hole.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, directed (lexicographically) from left to right.
   */
  Halfedge_handle insert_in_face_interior(const X_monotone_curve_2& cv,
                                          Face_handle f);

  /*!
   * Insert an x-monotone curve into the arrangement, such that its left
   * endpoint corresponds to a given arrangement vertex.
   * \param cv The given x-monotone curve.
   * \param v The given vertex.
   * \param f The face that contains v (in case it has no incident edges).
   * \pre The left endpoint of cv is incident to the vertex v.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target is the new vertex.
   */
  Halfedge_handle insert_from_left_vertex(const X_monotone_curve_2& cv,
                                          Vertex_handle v,
                                          Face_handle f = Face_handle());

  /*!
   * Insert an x-monotone curve into the arrangement, such that its left
   * endpoints corresponds to a given arrangement vertex, given the exact
   * place for the curve in the circular list around this vertex.
   * \param cv The given x-monotone curve.
   * \param prev The reference halfedge. We should represent cv as a pair
   *             of edges, one of them should become prev's successor.
   * \pre The target vertex of prev is cv's left endpoint.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target is the new vertex that was created.
   */
  Halfedge_handle insert_from_left_vertex(const X_monotone_curve_2& cv,
                                          Halfedge_handle prev);

  /*!
   * Insert an x-monotone curve into the arrangement, such that its right
   * endpoint corresponds to a given arrangement vertex.
   * \param cv The given x-monotone curve.
   * \param v The given vertex.
   * \param f The face that contains v (in case it has no incident edges).
   * \pre The right endpoint of cv is incident to the vertex v.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target is the new vertex.
   */
  Halfedge_handle insert_from_right_vertex(const X_monotone_curve_2& cv,
                                           Vertex_handle v,
                                           Face_handle f = Face_handle());

  /*!
   * Insert an x-monotone curve into the arrangement, such that its right
   * endpoints corresponds to a given arrangement vertex, given the exact
   * place for the curve in the circular list around this vertex.

   * \param cv The given x-monotone curve.
   * \param prev The reference halfedge. We should represent cv as a pair
   *             of edges, one of them should become prev's successor.
   * \pre The target vertex of prev is cv's right endpoint.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve, whose target is the new vertex that was created.
   */
  Halfedge_handle insert_from_right_vertex(const X_monotone_curve_2& cv,
                                           Halfedge_handle prev);

  /*!
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to given arrangement vertices.
   * \param cv The given x-monotone curve.
   * \param v1 The first vertex.
   * \param v2 The second vertex.
   * \param f The face that contains v1 and v2
   *          (in case both have no incident edges).
   * \pre v1 and v2 corresponds to cv's endpoints.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve directed from v1 to v2.
   */
  Halfedge_handle insert_at_vertices(const X_monotone_curve_2& cv,
                                     Vertex_handle v1,
                                     Vertex_handle v2,
                                     Face_handle f = Face_handle());

  /*!
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to given arrangement vertices, given the exact
   * place for the curve in one of the circular lists around a vertex.
   * \param cv The given x-monotone curve.
   * \param prev1 The reference halfedge for the first vertex.
   * \param v2 The second vertex.
   * \pre The target vertex of prev1 and v2 corresponds to cv's endpoints.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve directed from prev1 to v2.
   */
  Halfedge_handle insert_at_vertices(const X_monotone_curve_2& cv,
                                     Halfedge_handle prev1,
                                     Vertex_handle v2);

  /*!
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to given arrangement vertices, given the exact
   * place for the curve in both circular lists around these two vertices.
   * \param cv the given curve.
   * \param prev1 The reference halfedge for the first vertex.
   * \param prev2 The reference halfedge for the second vertex.
   * \pre The target vertices of prev1 and prev2 are cv's endpoints.
   * \return A handle for one of the halfedges corresponding to the inserted
   *         curve directed from prev1's target to prev2's target.
   */
  Halfedge_handle insert_at_vertices(const X_monotone_curve_2 & cv,
                                     Halfedge_handle prev1,
                                     Halfedge_handle prev2);

  //@}

  /// \name Vertex manipulation functions.
  //@{

  /*!
   * Replace the point associated with the given vertex.
   * \param v The vertex to modify.
   * \param p The point that should be associated with the edge.
   * \pre p is geometrically equivalent to the current point
   *      associated with v.
   * \return A handle for a the modified vertex (same as v).
   */
  Vertex_handle modify_vertex(Vertex_handle v, const Point_2& p);

  /*!
   * Remove an isolated vertex from the interior of a given face.
   * \param v The vertex to remove.
   * \pre v is an isolated vertex (it has no incident halfedges).
   * \return A handle for the face containing v.
   */
  Face_handle remove_isolated_vertex(Vertex_handle v);

  ///@}

  /// \name Halfedge manipulation functions.
  //@{

  /*!
   * Replace the x-monotone curve associated with the given edge.
   * \param e The edge to modify.
   * \param cv The curve that should be associated with the edge.
   * \pre cv is geometrically equivalent to the current curve
   *      associated with e.
   * \return A handle for a the modified halfedge (same as e).
   */
  Halfedge_handle modify_edge(Halfedge_handle e, const X_monotone_curve_2& cv);

  /*!
   * Split a given edge into two, and associate the given x-monotone
   * curves with the split edges.
   * \param e The edge to split (one of the pair of twin halfegdes).
   * \param cv1 The curve that should be associated with the first split edge.
   * \param cv2 The curve that should be associated with the second split edge.

   * \pre cv1's source and cv2's target equal the endpoints of the curve
   *      currently assoicated with e (respectively), and cv1's target equals
   *      cv2's target, and this is the split point (ot vice versa).
   * \return A handle for the halfedge whose source is the source of the the
   *         original halfedge e, and whose target is the split point.
   */
  Halfedge_handle split_edge(Halfedge_handle e,
                             const X_monotone_curve_2& cv1,
                             const X_monotone_curve_2& cv2);

  /*!
   * Merge two edges to form a single edge, and associate the given x-monotone
   * curve with the merged edge.
   * \param e1 The first edge to merge (one of the pair of twin halfegdes).
   * \param e2 The second edge to merge (one of the pair of twin halfegdes).
   * \param cv The curve that should be associated with merged edge.
   * \return A handle for the merged halfedge.
   */
  Halfedge_handle merge_edge(Halfedge_handle e1, Halfedge_handle e2,
                             const X_monotone_curve_2& cv);

  /*!
   * Remove an edge from the arrangement.
   * \param e The edge to remove (one of the pair of twin halfegdes).
   * \param remove_source Should the source vertex of e be removed if it
   *                      becomes isolated (true by default).
   * \param remove_target Should the target vertex of e be removed if it
   *                      becomes isolated (true by default).
   * \return A handle for the remaining face.
   */
  Face_handle remove_edge(Halfedge_handle e,
                          bool remove_source = true,
                          bool remove_target = true);

  //@}

protected:
  /// \name Determining the boundary-side conditions.
  //@{

  /*! Determines whether a boundary-side categoty indicates an open side.
   */
  inline bool is_open(Arr_boundary_side_tag) const { return false; }
  inline bool is_open(Arr_open_side_tag) const { return true; }

  /*! Determines whether the given x and y parameter spaces are open.
   * These parameter spaces are typically associated with a particular curve
   * end.
   * \param ps_x The parameter space in x.
   * \param ps_y The parameter space in y.
   */
  inline bool is_open(Arr_parameter_space ps_x, Arr_parameter_space ps_y) const
  {
    return
      (((ps_x == ARR_LEFT_BOUNDARY) && is_open(Left_side_category())) ||
       ((ps_x == ARR_RIGHT_BOUNDARY) && is_open(Right_side_category())) ||
       ((ps_y == ARR_BOTTOM_BOUNDARY) && is_open(Bottom_side_category())) ||
       ((ps_y == ARR_TOP_BOUNDARY) && is_open(Top_side_category())));

  }

  /*! Determines whether a boundary-side categoty indicates a constracted side.
   */
  inline bool is_contracted(Arr_boundary_side_tag) const { return false; }
  inline bool is_contracted(Arr_contracted_side_tag) const { return true; }

  /*! Determines whether a boundary-side categoty indicates a constracted side.
   */
  inline bool is_identified(Arr_boundary_side_tag) const { return false; }
  inline bool is_identified(Arr_identified_side_tag) const { return true; }
  //@}

  /// \name Allocating and de-allocating points and curves.
  //@{

  /*! Allocate a new point. */
  Point_2*_new_point(const Point_2& pt)
  {
    Point_2* p_pt = m_points_alloc.allocate(1);

    m_points_alloc.construct(p_pt, pt);
    return (p_pt);
  }

  /*! De-allocate a point. */
  void _delete_point(Point_2& pt)
  {
    Point_2* p_pt = &pt;

    m_points_alloc.destroy(p_pt);
    m_points_alloc.deallocate(p_pt, 1);
  }

  /*! Allocate a new curve. */
  X_monotone_curve_2* _new_curve(const X_monotone_curve_2& cv)
  {
    X_monotone_curve_2* p_cv = m_curves_alloc.allocate(1);
    m_curves_alloc.construct(p_cv, cv);
    return (p_cv);
  }

  /*! De-allocate a curve. */
  void _delete_curve(X_monotone_curve_2& cv)
  {
    X_monotone_curve_2* p_cv = &cv;

    m_curves_alloc.destroy(p_cv);
    m_curves_alloc.deallocate(p_cv, 1);
  }
  //@}

  /// \name Converting handles to pointers (for the arrangement accessor).
  //@{
  /*! Access the DCEL (non-const version). */
  inline Dcel& _dcel() { return (m_topol_traits.dcel()); }
  /*! Access the DCEL (const version). */
  inline const Dcel& _dcel() const
  { return (m_topol_traits.dcel()); }

  /*! Convert a vertex handle to a pointer to a DCEL vertex. */
  inline DVertex* _vertex(Vertex_handle vh) const
  { return (&(*vh)); }

  /*! Convert a constant vertex handle to a pointer to a DCEL vertex. */
  inline const DVertex* _vertex(Vertex_const_handle vh) const
  { return (&(*vh)); }

  /*! Convert a halfedge handle to a pointer to a DCEL halfedge. */
  inline DHalfedge* _halfedge(Halfedge_handle hh) const
  { return (&(*hh)); }

  /*! Convert a constant halfedge handle to a pointer to a DCEL halfedge. */
  inline const DHalfedge* _halfedge(Halfedge_const_handle hh) const
  { return (&(*hh)); }

  /*! Convert a face handle to a pointer to a DCEL face. */
  inline DFace* _face(Face_handle fh) const
  { return (&(*fh)); }

  /*! Convert a constant face handle to a pointer to a DCEL face. */
  inline const DFace* _face(Face_const_handle fh) const
  { return (&(*fh)); }
  //@}

  /// \name Converting pointers to handles (for the arrangement accessor).
  //@{

  /*! Convert a pointer to a DCEL vertex to a vertex handle. */
  Vertex_handle _handle_for(DVertex* v)
  { return (Vertex_handle(v)); }

  /*! Convert a pointer to a DCEL vertex to a constant vertex handle. */
  Vertex_const_handle _const_handle_for(const DVertex* v) const
  { return (Vertex_const_handle(v)); }

  /*! Convert a pointer to a DCEL halfedge to a halfedge handle. */
  Halfedge_handle _handle_for(DHalfedge* he)
  { return (Halfedge_handle(he)); }


  /*! Convert a pointer to a DCEL halfedge to a constant halfedge handle. */
  Halfedge_const_handle _const_handle_for(const DHalfedge* he) const
  { return (Halfedge_const_handle(he)); }

  /*! Convert a pointer to a DCEL face to a face handle. */
  Face_handle _handle_for(DFace* f)
  { return (Face_handle(f)); }

  /*! Convert a pointer to a DCEL face to a constant face handle. */
  Face_const_handle _const_handle_for(const DFace* f) const
  { return (Face_const_handle(f)); }
  //@}

  /// \name Auxiliary (protected) functions.
  //@{

  /*! Is the vertex incident to a given halfedge lexicographically smaller than
   * the vertex incident to another given halfedge. Recall that the incident
   * vertex is the target vertex. This function is used, for example, in the
   * search for lexicographically smallest vertex in a CCB, when an edge is
   * about to be removed from the DCEL.
   *
   * This is the implementation for the case where all 4 boundary sides are
   * oblivious.
   *
   * \param he1 the given first halfedge
   * \param ps_x1 the parameter space in x of the vertex incident to he1
   * \param ps_y1 the parameter space in y of the vertex incident to he1
   * \param he2 the given second halfedge
   * \param ps_x2 the parameter space in x of the vertex incident to he2
   * \param ps_y2 the parameter space in y of the vertex incident to he2
   * \precondition he1 is directed from right to left
   * \precondition he2 is directed from right to left
   * \precondition the vertex incident to he1 (he1->vertex()) is different
   *        than the vertex incident to he1 (he2->vertex()), and thus their
   *        geometric mappings (he1->vertex()->point() and
   *        he2->vertex()->point()) are not equal.
   */
  bool _is_smaller(const DHalfedge* he1,
                   Arr_parameter_space ps_x1, Arr_parameter_space ps_y1,
                   const DHalfedge* he2,
                   Arr_parameter_space ps_x2, Arr_parameter_space ps_y2,
                   Arr_all_sides_oblivious_tag) const;

  /*! This is a wrapper for the case where any boundary side is not
   * necessarily oblivious.
   */
  bool _is_smaller(const DHalfedge* he1,
                   Arr_parameter_space ps_x1, Arr_parameter_space ps_y1,
                   const DHalfedge* he2,
                   Arr_parameter_space ps_x2, Arr_parameter_space ps_y2,
                   Arr_not_all_sides_oblivious_tag) const;

  /*! Is the lexicographically minimal vertex of a given x-monotone curve
   * lexicographically smaller than the lexicographically minimal vertex of
   * another given x-monotone curve. This function is used, for example, when
   * a new curve is to be inserted into the arrangement. In this case the
   * search is conducted over the curves that will comprise a new CCB.
   *
   * This is the implementation for the case where all 4 boundary sides are
   * oblivious.
   *
   * \param cv1 the given first x-monotone curve
   * \param ps_x1 the parameter space in x of the minimal point of cv1
   * \param ps_y1 the parameter space in y of the minimal point of cv1
   * \param cv2 the given second x-monotone curve
   * \param ps_x2 the parameter space in x of the minimal point of cv2
   * \param ps_y2 the parameter space in y of the minimal point of cv2
   * \precondition the minimal points of cv1 and cv2 are not equal.
   */
  bool _is_smaller(const X_monotone_curve_2& cv1, const Point_2& p1,
                   Arr_parameter_space ps_x1, Arr_parameter_space ps_y1,
                   const X_monotone_curve_2& cv2, const Point_2& p2,
                   Arr_parameter_space ps_x2, Arr_parameter_space ps_y2,
                   Arr_all_sides_oblivious_tag) const;

  /*! This is the implementation for the case where any one of the 4 boundary
   * sides can be of any type.
   */
  bool _is_smaller(const X_monotone_curve_2& cv1, const Point_2& p1,
                   Arr_parameter_space ps_x1, Arr_parameter_space ps_y1,
                   const X_monotone_curve_2& cv2, const Point_2& p2,
                   Arr_parameter_space ps_x2, Arr_parameter_space ps_y2,
                   Arr_not_all_sides_oblivious_tag) const;

  /*! Given two x-monotone curves that share their minimal end point.
   * The function return true if the y-coordinate of the first curve curve
   * near its minimal end smaller than the y-coordinate of the second curve
   * (near its minimal end). This function is used, for example, when
   * a new curve is to be inserted into the arrangement. In this case the
   * search is conducted over the curves that will comprise a new CCB.
   *
   * This is the implementation for the case where all 4 boundary sides are
   * oblivious.
   *
   * \param cv1 the given first x-monotone curve
   * \param cv2 the given second x-monotone curve
   * \param p the shared minimal point of cv1 and cv2
   * \param ps_x the parameter space in x of the minimal point of cv1
   * \param ps_y the parameter space in y of the minimal point of cv1
   * \precondition the minimal points of cv1 and cv2 are equal.
   */
  bool _is_smaller_near_right(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              const Point_2& p,
                              Arr_parameter_space ps_x,
                              Arr_parameter_space ps_y,
                              Arr_all_sides_oblivious_tag) const;

  /*! This is the implementation for the case where any one of the 4 boundary
   * sides can be of any type.
   */
  bool _is_smaller_near_right(const X_monotone_curve_2& cv1,
                              const X_monotone_curve_2& cv2,
                              const Point_2& p,
                              Arr_parameter_space ps_x,
                              Arr_parameter_space ps_y,
                              Arr_not_all_sides_oblivious_tag) const;

  /*!
   * Locate the place for the given curve around the given vertex.
   * \param v The given arrangement vertex.
   * \param cv The given x-monotone curve.
   * \param ind Whether we refer to the minimal or maximal end of cv.
   * \return A pointer to a halfedge whose target is v, where cv should be
   *         inserted between this halfedge and the next halfedge around this
   *         vertex (in a clockwise order).
   *         A NULL return value indicates a precondition violation.
   */
  DHalfedge* _locate_around_vertex(DVertex* v, const X_monotone_curve_2& cv,
                                   Arr_curve_end ind) const;

  /*!
   * Compute the distance (in halfedges) between two halfedges.
   * \param e1 The source halfedge.
   * \param e2 The destination halfedge.
   * \pre e1 and e2 belong to the same connected component
   * \return The number of halfedges along the component boundary between the
   *         two halfedges.
   */
  unsigned int _halfedge_distance(const DHalfedge* e1,
                                  const DHalfedge* e2) const;

  /*!
   * Compare the length of the induced paths from e1 to e2 and
   *  from e2 to e1.
   * \pre e1 and e2 belong to the same connected component
   * \return The comparison result
   */
  Comparison_result _compare_induced_path_length(const DHalfedge* e1,
                                                 const DHalfedge* e2) const;

  /*!
   * Update the indices according to boundary locations
   */
  void
  _compute_indices(Arr_parameter_space ps_x_curr, Arr_parameter_space ps_y_curr,
                   Arr_parameter_space ps_x_next, Arr_parameter_space ps_y_next,
                   int& x_index, int& y_index,  boost::mpl::bool_<true>) const;

  /*!
   * Update the indices according to boundary locations (i.e. does nothing)
   */
  void
  _compute_indices(Arr_parameter_space ps_x_curr, Arr_parameter_space ps_y_curr,
                   Arr_parameter_space ps_x_next, Arr_parameter_space ps_y_next,
                   int& x_index, int& y_index,  boost::mpl::bool_<false>) const;

  /*!
   * Is the first given x-monotone curve above the second given?
   * \param xcv1 the first given curve
   * \param ps_y1 the parameter space in y of xcv1
   * \param xcv2 the second given curve
   * \param Arr_identified_side_tag used for dispatching to ensure that this
   *        function is invoked when the bottom and top boundaries are
   *        identified
   */
  bool _is_above(const X_monotone_curve_2& xcv1,
                 const X_monotone_curve_2& xcv2,
                 const Point_2& point,
                 Arr_parameter_space ps_y1,
                 Arr_has_identified_side_tag) const;

  /*!
   * Is the first given x-monotone curve above the second given?
   * \param xcv1 the first given curve
   * \param ps_y1 the parameter space in y of xcv1
   * \param xcv2 the second given curve
   * \param Arr_contracted_side_tag used for dispatching to ensure that this
   *        function is invoked when the bottom or top boundaries are
   *        contracted
   */
  bool _is_above(const X_monotone_curve_2& xcv1,
                 const X_monotone_curve_2& xcv2,
                 const Point_2& point,
                 Arr_parameter_space ps_y1,
                 Arr_has_contracted_side_tag) const;

  /*!
   * Is the first given x-monotone curve above the second given?
   * \param xcv1 the first given curve
   * \param ps_y1 the parameter space in y of xcv1
   * \param xcv2 the second given curve
   * \param Arr_oblivious_side_tag used for dispatching to ensure that this
   *        function is invoked when the bottom and top boundaries are neither
   *        identified nor contracted
   */
  bool _is_above(const X_monotone_curve_2& xcv1,
                 const X_monotone_curve_2& xcv2,
                 const Point_2& point,
                 Arr_parameter_space ps_y1,
                 Arr_boundary_cond_tag) const;

  /*!
   * Compute the signs (in left/right and bottom/top) of a path
   * induced by the sequence he_to=>cv,cv_dir=>he_away, and reports
   * as side-effect the halfedges pointing to local minima copied
   * to an outputiterator.
   * \param he_to The predecessor halfedge.
   * \param cv The x-monotone curve we use to connect he_to's target and
   *           he_away's source vertex.
   * \param cv_dir the direction of the curve between he_to and he_away
   * \param he_away The succcessor halfedge.
   * \param local_mins_it the outputiterator
   * (value_type = std::pair< DHalfedge*, int >, where the int denotes the
   * index) to report the halfedges pointing to local minima (<-shaped
   * situation)
   * \return A pair of signs for the induced path (ZERO if non-perimetric,
   * POSITIVE if perimetric ccb is oriented in positive direction,
   * NEGATIVE if perimetric ccb is oriented in negative direction).
   */
  template <typename OutputIterator>
  std::pair<Sign, Sign>
  _compute_signs_and_local_minima(const DHalfedge* he_to,
                                  const X_monotone_curve_2& cv,
                                  Arr_halfedge_direction cv_dir,
                                  const DHalfedge* he_away,
                                  OutputIterator local_mins_it) const;

  /*!
   * Compute the signs (in left/right and bottom/top) of a closed ccb (loop)
   * represented by a given halfedge, and the halfedge pointing to the smallest
   * vertex on the ccb.
   * \param he The representative halfedge on the ccb.
   * \param ps_x_min The parameter space in x of the smallest vertex.
   * \param ps_y_min The parameter space in y of the smallest vertex.
   * \param index_min The index of the smallest vertex.
   * \return A pair of, a pair of signs for the induced path, and the halfedge
   *     pointing to the smallest vertex.
   *     A sign ZERO is if the ccb is non-perimetric,
   *     POSITIVE if the ccb is perimetric and oriented in positive direction,
   *     NEGATIVE if the ccb is perimetric and oriented in negative direction).
   */
  std::pair<std::pair<Sign, Sign>,  const DHalfedge*>
  _compute_signs_and_min(const DHalfedge* he,
                         Arr_parameter_space& ps_x_min,
                         Arr_parameter_space& ps_y_min,
                         int& index_min) const;

  /*!
   * Compute the signs (in left/right and bottom/top) of a closed ccb (loop)
   * represented by a given halfedge.
   * \param he The representative halfedge on the ccb.
   * \return A pair of signs for the induced path.
   *     A sign ZERO is if the ccb is non-perimetric,
   *     POSITIVE if the ccb is perimetric and oriented in positive direction,
   *     NEGATIVE if the ccb is perimetric and oriented in negative direction).
   */
  std::pair<Sign, Sign> _compute_signs(const DHalfedge* he,
                                       boost::mpl::bool_<true>) const;

  /*! Compute the signs (in left/right and bottom/top) of a closed ccb (loop)
   * represented by a given halfedge for the case where non of the boundaries
   * is identified.
   * \return the pair (ZERO, ZERO)
   */
  std::pair<Sign, Sign> _compute_signs(const DHalfedge* he,
                                       boost::mpl::bool_<false>) const;

  /*!
   * Given two predecessor halfedges that will be used for inserting a
   * new halfedge pair (he_to is the predecessor of the directed curve
   * cv, cv_dir and he_away will be the successor), such that the
   * insertion will create a new face that forms a hole inside an existing
   * face, determine whether he_to=>cv,cv_dir=>he_away will be part
   * of the new outer ccb of the new face.
   * \param he_to The predecessor halfedge.
   * \param cv The x-monotone curve we use to connect he_to's target and
   *           he_away's source vertex.
   * \param cv_dir the direction of the curve between he_to and he_away
   * \param he_away The succcessor halfedge.
   * \pre he_to and he_away belong to the same inner CCB.
   * \return true if he_to=>cv,cv_dir=>he_away lie in the interior of the face we
   *         are about to create (i.e.~are part of the new outer ccb),
   *         false otherwise - in which case the subsequence
   *         he_away->next()=>cv,opposite(cv_dir)=>he_to->next()
   *         must be incident to this new face (i.e.~are part
   *         of the new outer ccb).
   */
  template <typename InputIterator>
  bool _defines_outer_ccb_of_new_face(const DHalfedge* he_to,
                                      const X_monotone_curve_2& cv,
                                      const DHalfedge* he_away,
                                      InputIterator lm_begin,
                                      InputIterator lm_end) const;

  /*!
   * Move a given outer CCB from one face to another.
   * \param from_face The face currently containing the component.
   * \param to_face The face into which we should move the component.
   * \param he A halfedge lying on the outer component.
   */
  void _move_outer_ccb(DFace* from_face, DFace* to_face, DHalfedge* he);

  /*!
   * Move a given inner CCB (hole) from one face to another.
   * \param from_face The face currently containing the component.
   * \param to_face The face into which we should move the component.
   * \param he A halfedge lying on the inner component.
   */
  void _move_inner_ccb(DFace* from_face, DFace* to_face, DHalfedge* he);

  /*!
   * Move all inner CCBs (holes) from one face to another.
   * \param from_face The face currently containing the components.
   * \param to_face The face into which we should move the components.
   */
  void _move_all_inner_ccb(DFace* from_face, DFace* to_face);

  /*!
   * Insert the given vertex as an isolated vertex inside the given face.
   * \param f The face that should contain the isolated vertex.
   * \param v The isolated vertex.
   */
  void _insert_isolated_vertex(DFace* f, DVertex* v);

  /*!
   * Move a given isolated vertex from one face to another.
   * \param from_face The face currently containing the isolated vertex.
   * \param to_face The face into which we should move the isolated vertex.
   * \param v The isolated vertex.
   */
  void _move_isolated_vertex(DFace* from_face, DFace* to_face, DVertex* v);

  /*!
   * Move all isolated vertices from one face to another.
   * \param from_face The face currently containing the isolated vertices.
   * \param to_face The face into which we should move the isolated vertices.
   */
  void _move_all_isolated_vertices(DFace* from_face, DFace* to_face);

  /*!
   * Create a new vertex and associate it with the given point.
   * \param p The point.
   * \return A pointer to the newly created vertex.
   */
  DVertex* _create_vertex(const Point_2& p);

  /*!
   * Create a new boundary vertex.
   * \param cv The curve incident to the boundary.
   * \param ind The relevant curve-end.
   * \param bx The boundary condition in x.
   * \param by The boundary condition in y.
   * \pre Either bx or by does not equal ARR_INTERIOR.
   * \return A pointer to the newly created vertex.
   */
  DVertex* _create_boundary_vertex(const X_monotone_curve_2& cv,
                                   Arr_curve_end ind,
                                   Arr_parameter_space bx,
                                   Arr_parameter_space by);

  /*!
   * Locate the DCEL features that will be used for inserting the given curve
   * end, which has a boundary condition, and set a proper vertex there.
   * \param f The face that contains the curve end.
   * \param cv The x-monotone curve.
   * \param ind The curve end.
   * \param bx The boundary condition at the x-coordinate.
   * \param by The boundary condition at the y-coordinate.
   * \param p_pred Output: The predecessor halfedge around this vertex
   *                       (may be NULL, if no such halfedge exists).
   * \return The vertex that corresponds to the curve end.
   */
  DVertex* _place_and_set_curve_end(DFace* f,
                                    const X_monotone_curve_2& cv,
                                    Arr_curve_end ind,
                                    Arr_parameter_space bx,
                                    Arr_parameter_space by,
                                    DHalfedge** p_pred);

  /*!
   * Insert an x-monotone curve into the arrangement, such that both its
   * endpoints correspond to free arrangement vertices (newly created vertices
   * or existing isolated vertices), so a new inner CCB is formed in the face
   * that contains the two vertices.
   * \param f The face containing the two end vertices.
   * \param cv The given x-monotone curve.
   * \param cv_dir The direction of the curve
   * \param v1 The free vertex that corresponds to the left endpoint of cv.
   * \param v2 The free vertex that corresponds to the right endpoint of cv.
   * \return A pointer to one of the halfedges corresponding to the inserted
   *         curve, directed from v1 to v2.
   */
  DHalfedge* _insert_in_face_interior(DFace* f,
                                      const X_monotone_curve_2& cv,
                                      Arr_halfedge_direction cv_dir,
                                      DVertex* v1, DVertex* v2);

  /*!
   * Insert an x-monotone curve into the arrangement, such that one of its
   * endpoints corresponds to a given arrangement vertex, given the exact
   * place for the curve in the circular list around this vertex. The other
   * endpoint corrsponds to a free vertex (a newly created vertex or an
   * isolated vertex).
   * \param he_to The reference halfedge. We should represent cv as a pair
   *              of edges, one of them should become he_to's successor.
   * \param cv The given x-monotone curve.
   * \param cv_dir The direction of cv.
   * \param v The free vertex that corresponds to the other endpoint.
   * \return A pointer to one of the halfedges corresponding to the inserted
   *         curve, whose target is the vertex v.
   */
  DHalfedge* _insert_from_vertex(DHalfedge* he_to, const X_monotone_curve_2& cv,
                                 Arr_halfedge_direction cv_dir,
                                 DVertex* v);

  /*!
   * Insert an x-monotone curve into the arrangement, where the end vertices
   * are given by the target points of two given halfedges.
   * The two halfedges should be given such that in case a new face is formed,
   * it will be the incident face of the halfedge directed from the first
   * vertex to the second vertex.
   * \param he_to The reference halfedge pointing to the insertion vertex
   * \param cv the given curve.
   * \param cv_dir the direction of the curve
   * \param he_away the reference halfedge for the second vertex.
   * \param res the comparison result of the points associated with prev1's
   *            target vertex and prev2's target vertex.
   * \param new_face (Output) indicates whether a new face has been created.
   * \param swapped_predecessors (Output) indicates whether roles of prev1 and
   *                                      prev2 have been switched
   * \param allow_swap_of_predecessors set to false if no swapping should
   *                                   take place at all
   * \return A pointer to one of the halfedges corresponding to the inserted
   *         curve directed from prev1's target to prev2's target.
   *         In case a new face has been created, it is given as the incident
   *         face of this halfedge.
   */
  DHalfedge* _insert_at_vertices(DHalfedge* he_to,
                                 const X_monotone_curve_2& cv,
                                 Arr_halfedge_direction cv_dir,
                                 DHalfedge* he_away,
                                 bool& new_face,
                                 bool& swapped_predecessors,
                                 bool allow_swap_of_predecessors = true);

  /*!
   * Relocate all inner CCBs and isolated vertices to their proper position,
   * immediately after a face has split due to the insertion of a new halfedge.
   * \param new_he The new halfedge that caused the split, such that the new
   *               face lies to its left and the old face to its right.
   */
  void _relocate_in_new_face(DHalfedge* new_he);

  /*!
   * Relocate all inner CCBs to their proper position,
   * immediately after a face has split due to the insertion of a new halfedge.
   * \param new_he The new halfedge that caused the split, such that the new
   *               face lies to its left and the old face to its right.
   */
  void _relocate_inner_ccbs_in_new_face(DHalfedge* new_he);

  /*!
   * Relocate all vertices to their proper position,
   * immediately after a face has split due to the insertion of a new halfedge.
   * \param new_he The new halfedge that caused the split, such that the new
   *               face lies to its left and the old face to its right.
   */
  void _relocate_isolated_vertices_in_new_face(DHalfedge* new_he);

  /*!
   * Replace the point associated with the given vertex.
   * \param v The vertex to modify.
   * \param p The point that should be associated with the edge.
   */
  void _modify_vertex(DVertex* v, const Point_2& p);

  /*!
   * Replace the x-monotone curve associated with the given edge.
   * \param e The edge to modify.
   * \param cv The curve that should be associated with the edge.
   */
  void _modify_edge(DHalfedge* he, const X_monotone_curve_2& cv);

  /*!
   * Check if the given vertex represents one of the ends of a given curve.
   * \param v The vertex.
   * \param cv The curve.
   * \param ind Indicates whether the minimal or the maximal end of cv is
   *            refereed to.
   * \return Whether v represents the left (or right) end of cv.
   */
  bool _are_equal(const DVertex* v,
                  const X_monotone_curve_2& cv, Arr_curve_end ind) const;

  /*!
   * Split a given edge into two at a given point, and associate the given
   * x-monotone curves with the split edges.
   * \param e The edge to split (one of the pair of twin halfegdes).
   * \param p The split point.
   * \param cv1 The curve that should be associated with the first split edge,
   *            whose source equals e's source and its target is p.
   * \param cv2 The curve that should be associated with the second split edge,
   *            whose source is p and its target equals e's target.
   * \return A pointer to the first split halfedge, whose source equals the
   *         source of e, and whose target is the split point.
   */
  DHalfedge* _split_edge(DHalfedge* e, const Point_2& p,
                         const X_monotone_curve_2& cv1,
                         const X_monotone_curve_2& cv2);

  /*!
   * Split a given edge into two at a given vertex, and associate the given
   * x-monotone curves with the split edges.
   * \param e The edge to split (one of the pair of twin halfegdes).
   * \param v The split vertex.
   * \param cv1 The curve that should be associated with the first split edge,
   *            whose source equals e's source and its target is v.
   * \param cv2 The curve that should be associated with the second split edge,
   *            whose source is v and its target equals e's target.
   * \return A pointer to the first split halfedge, whose source equals the
   *         source of e, and whose target is v.
   */
  DHalfedge* _split_edge(DHalfedge* e, DVertex* v,
                         const X_monotone_curve_2& cv1,
                         const X_monotone_curve_2& cv2);

  /*!
   * Remove a pair of twin halfedges from the arrangement.
   * \param e One of the halfedges to be removed.
   * \param remove_source Should the source vertex of e be removed if it
   *                      becomes isolated.
   * \param remove_target Should the target vertex of e be removed if it
   *                      becomes isolated.
   * \pre In case the removal causes the creation of a new inner CCB (hole),
   *      e should point at this hole.
   * \return A pointer to the remaining face.
   */
  DFace* _remove_edge(DHalfedge* e, bool remove_source, bool remove_target);

  /*!
   * Decide whether a hole is created when an edge is removed.
   *
   * \param signs1 signs of future ccb1
   * \param signs2 signs of future ccb2
   * \param same_face to he and he->opposite() belong to same face
   * return true, in case a new hole is created, false otherwise
   */
  bool _hole_creation_on_edge_removal(std::pair< CGAL::Sign, CGAL::Sign > signs1,
                                      std::pair< CGAL::Sign, CGAL::Sign > signs2,
                                      bool same_face);

  /*!
   * Remove a vertex in case it becomes redundant after the deletion of an
   * incident edge.
   * \param v The vertex.
   * \param f The face that contains v (in case it becomes isolated).
   */
  void _remove_vertex_if_redundant(DVertex* v, DFace* f);

  /*!
   * Remove an isolated vertex from the interior of its face (but not from
   * the DCEL).
   * \param v The isolated vertex to remove.
   */
  void _remove_isolated_vertex(DVertex* v);
  //@}

  /// \name Auxiliary (protected) functions for validity checking.
  //@{

  /*! Check the validity of a given vertex. */
  bool _is_valid(Vertex_const_handle v) const;

  /*! Check the validity of a given halfedge. */
  bool _is_valid(Halfedge_const_handle he) const;

  /*! Check the validity of a given face. */
  bool _is_valid(Face_const_handle f) const;

  /*! Check the validity of an outer CCB. */
  bool _is_outer_ccb_valid(const DOuter_ccb* oc, const DHalfedge* first) const;

  /*! Check the validity of an inner CCB. */
  bool _is_inner_ccb_valid(const DInner_ccb* ic, const DHalfedge* first) const;

  /*!
   * Check that all vertices are unique (no two vertices with the same
   * geometric point.
   */
  bool _are_vertices_unique() const;

  /*! Check that the curves around a given vertex are ordered clockwise. */
  bool _are_curves_ordered_cw_around_vertrex(Vertex_const_handle v) const;

  //@}

protected:
  /// \name Managing and notifying the arrangement observers.
  //@{

  /*!
   * Register a new observer (so it starts receiving notifications).
   * \param p_obs A pointer to the observer object.
   */
  void _register_observer(Observer* p_obs) { m_observers.push_back(p_obs); }

  /*!
   * Unregister a new observer (so it stops receiving notifications).
   * \param p_obs A pointer to the observer object.
   * \return Whether the observer was successfully unregistered.
   */
  bool _unregister_observer(Observer* p_obs)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();

    for (iter = m_observers.begin(); iter != end; ++iter) {
      if ((*iter) == p_obs) {
        // Remove the p_ob pointer from the list of observers.
        m_observers.erase (iter);
        return true;
      }
    }

    // If we reached here, the observer was not registered.
    return false;
  }

protected:
  /* Notify the observers on global arrangement operations: */

  void _notify_before_assign(const Self& arr)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_assign(arr);
  }

  void _notify_after_assign()
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_assign();
  }

  void _notify_before_clear()
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_clear();
  }

  void _notify_after_clear()
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();

    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_clear();
  }

  void _notify_before_global_change()
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_global_change();
  }

  void _notify_after_global_change()
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_global_change();
  }

  /* Notify the observers on local changes in the arrangement: */

  void _notify_before_create_vertex(const Point_2& p)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_create_vertex(p);
  }

  void _notify_after_create_vertex(Vertex_handle v)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_create_vertex(v);
  }

  void _notify_before_create_boundary_vertex(const X_monotone_curve_2& cv,
                                             Arr_curve_end ind,
                                             Arr_parameter_space bx,
                                             Arr_parameter_space by)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();

    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_create_boundary_vertex(cv, ind, bx, by);
  }

  void _notify_after_create_boundary_vertex(Vertex_handle v)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_create_boundary_vertex(v);
  }

  void _notify_before_create_edge(const X_monotone_curve_2& c,
                                  Vertex_handle v1, Vertex_handle v2)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_create_edge(c, v1, v2);
  }

  void _notify_after_create_edge(Halfedge_handle e)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_create_edge(e);
  }

  void _notify_before_modify_vertex(Vertex_handle v, const Point_2& p)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_modify_vertex(v, p);
  }

  void _notify_after_modify_vertex(Vertex_handle v)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_modify_vertex(v);
  }

  void _notify_before_modify_edge(Halfedge_handle e,
                                  const X_monotone_curve_2& c)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_modify_edge(e, c);
  }

  void _notify_after_modify_edge(Halfedge_handle e)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_modify_edge(e);
  }

  void _notify_before_split_edge(Halfedge_handle e, Vertex_handle v,
                                 const X_monotone_curve_2& c1,
                                 const X_monotone_curve_2& c2)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_split_edge(e, v, c1, c2);
  }

  void _notify_after_split_edge(Halfedge_handle e1, Halfedge_handle e2)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_split_edge(e1, e2);
  }

  void _notify_before_split_fictitious_edge(Halfedge_handle e, Vertex_handle v)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_split_fictitious_edge(e, v);
  }

  void _notify_after_split_fictitious_edge(Halfedge_handle e1,
                                           Halfedge_handle e2)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_split_fictitious_edge(e1, e2);
  }

  void _notify_before_split_face(Face_handle f, Halfedge_handle e)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_split_face(f, e);
  }

  void _notify_after_split_face(Face_handle f, Face_handle new_f, bool is_hole)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_split_face(f, new_f, is_hole);
  }

  void _notify_before_split_outer_ccb(Face_handle f, Ccb_halfedge_circulator h,
                                      Halfedge_handle e)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_split_outer_ccb(f, h, e);
  }

  void _notify_after_split_outer_ccb(Face_handle f, Ccb_halfedge_circulator h1,
                                     Ccb_halfedge_circulator h2)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_split_outer_ccb(f, h1, h2);
  }

  void _notify_before_split_inner_ccb(Face_handle f, Ccb_halfedge_circulator h,
                                      Halfedge_handle e)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_split_inner_ccb(f, h, e);
  }

  void _notify_after_split_inner_ccb(Face_handle f,
                                     Ccb_halfedge_circulator h1,
                                     Ccb_halfedge_circulator h2)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_split_inner_ccb(f, h1, h2);
  }

  void _notify_before_add_outer_ccb(Face_handle f, Halfedge_handle e)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_add_outer_ccb(f, e);
  }

  void _notify_after_add_outer_ccb(Ccb_halfedge_circulator h)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_add_outer_ccb(h);
  }

  void _notify_before_add_inner_ccb(Face_handle f, Halfedge_handle e)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_add_inner_ccb(f, e);
  }

  void _notify_after_add_inner_ccb(Ccb_halfedge_circulator h)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_add_inner_ccb(h);
  }

  void _notify_before_add_isolated_vertex(Face_handle f, Vertex_handle v)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_add_isolated_vertex(f, v);
  }

  void _notify_after_add_isolated_vertex(Vertex_handle v)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_add_isolated_vertex(v);
  }

  void _notify_before_merge_edge(Halfedge_handle e1, Halfedge_handle e2,
                                 const X_monotone_curve_2& c)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_merge_edge(e1, e2, c);
  }

  void _notify_after_merge_edge(Halfedge_handle e)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_merge_edge(e);
  }

  void _notify_before_merge_fictitious_edge(Halfedge_handle e1,
                                            Halfedge_handle e2)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_merge_fictitious_edge(e1, e2);
  }

  void _notify_after_merge_fictitious_edge(Halfedge_handle e)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_merge_fictitious_edge(e);
  }

  void _notify_before_merge_face(Face_handle f1, Face_handle f2,
                                 Halfedge_handle e)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_merge_face(f1, f2, e);
  }

  void _notify_after_merge_face(Face_handle f)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_merge_face(f);
  }

  void _notify_before_merge_outer_ccb(Face_handle f,
                                      Ccb_halfedge_circulator h1,
                                      Ccb_halfedge_circulator h2,
                                      Halfedge_handle e)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_merge_outer_ccb(f, h1, h2, e);
  }

  void _notify_after_merge_outer_ccb(Face_handle f,
                                     Ccb_halfedge_circulator h)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_merge_outer_ccb(f, h);
  }

  void _notify_before_merge_inner_ccb(Face_handle f,
                                      Ccb_halfedge_circulator h1,
                                      Ccb_halfedge_circulator h2,
                                      Halfedge_handle e)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_merge_inner_ccb(f, h1, h2, e);
  }

  void _notify_after_merge_inner_ccb(Face_handle f,
                                     Ccb_halfedge_circulator h)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_merge_inner_ccb(f, h);
  }

  void _notify_before_move_outer_ccb(Face_handle from_f,
                                     Face_handle to_f,
                                     Ccb_halfedge_circulator h)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_move_outer_ccb(from_f, to_f, h);
  }

  void _notify_after_move_outer_ccb(Ccb_halfedge_circulator h)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_move_outer_ccb(h);
  }

  void _notify_before_move_inner_ccb(Face_handle from_f,
                                     Face_handle to_f,
                                     Ccb_halfedge_circulator h)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_move_inner_ccb(from_f, to_f, h);
  }

  void _notify_after_move_inner_ccb(Ccb_halfedge_circulator h)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_move_inner_ccb(h);
  }

  void _notify_before_move_isolated_vertex(Face_handle from_f,
                                           Face_handle to_f,
                                           Vertex_handle v)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_move_isolated_vertex(from_f, to_f, v);
  }


  void _notify_after_move_isolated_vertex(Vertex_handle v)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_move_isolated_vertex(v);
  }

  void _notify_before_remove_vertex(Vertex_handle v)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_remove_vertex(v);
  }

  void _notify_after_remove_vertex()
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_remove_vertex();
  }

  void _notify_before_remove_edge(Halfedge_handle e)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_remove_edge(e);
  }

  void _notify_after_remove_edge()
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_remove_edge();
  }

  void _notify_before_remove_outer_ccb(Face_handle f, Ccb_halfedge_circulator h)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_remove_outer_ccb(f, h);
  }

  void _notify_after_remove_outer_ccb(Face_handle f)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_remove_outer_ccb(f);
  }

  void _notify_before_remove_inner_ccb(Face_handle f, Ccb_halfedge_circulator h)
  {
    Observers_iterator iter;
    Observers_iterator end = m_observers.end();
    for (iter = m_observers.begin(); iter != end; ++iter)
      (*iter)->before_remove_inner_ccb(f, h);
  }

  void _notify_after_remove_inner_ccb(Face_handle f)
  {
    Observers_rev_iterator iter;
    Observers_rev_iterator end = m_observers.rend();
    for (iter = m_observers.rbegin(); iter != end; ++iter)
      (*iter)->after_remove_inner_ccb(f);
  }
  //@}
};

//-----------------------------------------------------------------------------
// Declarations of the various global insertion and removal functions.
//-----------------------------------------------------------------------------

// In some compilers there is a template deduction disambiguity between this
// function and the following function receiving two InputIterator.
// For now the solution is to add a dummy variable at the end (referring
// to point-location). Maybe the proper solution is to use boost::enable_if
// together with appropriate tag.
/*!
 * Insert a curve or x-monotone curve into the arrangement (incremental
 * insertion).
 * The inserted curve can be x-monotone (or not) and may intersect the
 * existing arrangement.
 * \param arr The arrangement.
 * \param cv The curve to be inserted.
 * \param pl A point-location object associated with the arrangement.
 */
template <typename GeomTraits, typename TopTraits, typename Curve,
          typename PointLocation>
void insert(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
            const Curve& c, const PointLocation& pl,
            typename PointLocation::Point_2* = 0);

/*!
 * Insert a curve or x-monotone curve into the arrangement (incremental
 * insertion).
 * The inserted curve can be x-monotone (or not) and may intersect the
 * existing arrangement. The default "walk" point-location strategy is used
 * for the curve insertion.
 * \param arr The arrangement.
 * \param cv The curve to be inserted.
 */
template <typename GeomTraits, typename TopTraits, typename Curve>
void insert(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
            const Curve& c);

/*!
 * Insert a range of curves or x-monotone curves into the arrangement
 * (aggregated insertion).
 * The inserted curves may intersect one another and may also intersect the
 * existing arrangement.
 * \param arr The arrangement.
 * \param begin An iterator for the first curve in the range.
 * \param end A past-the-end iterator for the curve range.
 * \pre The value type of the iterators must be Curve_2.
 */
template <typename GeomTraits, typename TopTraits, typename InputIterator>
void insert(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
            InputIterator begin, InputIterator end);

/*!
 * Insert an x-monotone curve into the arrangement (incremental insertion)
 * when the location of the left endpoint of the curve is known and is
 * given as an isertion hint.
 * The inserted x-monotone curve may intersect the existing arrangement.
 * \param arr The arrangement.
 * \param cv The x-monotone curve to be inserted.
 * \param obj An object that represents the location of cv's left endpoint
 *            in the arrangement.
 */

template <typename GeomTraits, typename TopTraits>
void insert(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
            const typename GeomTraits::X_monotone_curve_2& c,
            const Object& obj);

/*!
 * Insert an x-monotone curve into the arrangement, such that the curve
 * interior does not intersect with any existing edge or vertex in the
 * arragement (incremental insertion).
 * \param arr The arrangement.
 * \param c The x-monotone curve to be inserted.
 * \param pl A point-location object associated with the arrangement.
 * \pre The interior of c does not intersect any existing edge or vertex.
 * \return A handle for one of the new halfedges corresponding to the
 *         inserted curve, directed (lexicographically) from left to right.
 */
template <typename GeomTraits, typename TopTraits, typename PointLocation>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
insert_non_intersecting_curve
(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
 const typename GeomTraits::X_monotone_curve_2& c,
 const PointLocation& pl);

/*!
 * Insert an x-monotone curve into the arrangement, such that the curve
 * interior does not intersect with any existing edge or vertex in the
 * arragement (incremental insertion). The default point-location strategy
 * is used for the curve insertion.
 * \param arr The arrangement.
 * \param c The x-monotone curve to be inserted.
 * \pre The interior of c does not intersect any existing edge or vertex.
 * \return A handle for one of the new halfedges corresponding to the inserted
 *         curve, directed (lexicographically) from left to right.
 */
template <typename GeomTraits, typename TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
insert_non_intersecting_curve
(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
 const typename GeomTraits::X_monotone_curve_2& c);

/*!
 * Insert a range of pairwise interior-disjoint x-monotone curves into
 * the arrangement, such that the curve interiors do not intersect with
 * any existing edge or vertex in the arragement (aggregated insertion).
 * \param arr The arrangement.
 * \param begin An iterator for the first x-monotone curve in the range.
 * \param end A past-the-end iterator for the x-monotone curve range.
 * \pre The value type of the iterators must be X_monotone_curve_2.
 *      The curves in the range are pairwise interior-disjoint, and their
 *      interiors do not intersect any existing edge or vertex.
 */
template <typename GeomTraits, typename TopTraits, typename InputIterator>
void insert_non_intersecting_curves
(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
 InputIterator begin, InputIterator end);

/*!
 * Remove an edge from the arrangement. In case it is possible to merge
 * the edges incident to the end-vertices of the removed edge after its
 * deletion, the function performs these merges as well.
 * \param arr The arrangement.
 * \param e The edge to remove (one of the pair of twin halfegdes).
 * \return A handle for the remaining face.
 */
template <typename GeomTraits, typename TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Face_handle
remove_edge(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
            typename Arrangement_on_surface_2<GeomTraits,
                                              TopTraits>::Halfedge_handle e);

/*!
 * Insert a vertex that corresponds to a given point into the arrangement.
 * The inserted point may lie on any existing arrangement feature.
 * \param arr The arrangement.
 * \param p The point to be inserted.
 * \param pl A point-location object associated with the arrangement.
 * \return A handle to the vertex that corresponds to the given point.
 */
template <typename GeomTraits, typename TopTraits, typename PointLocation>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Vertex_handle
insert_point(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
             const typename GeomTraits::Point_2& p,
             const PointLocation& pl);

/*!
 * Insert a vertex that corresponds to a given point into the arrangement.
 * The inserted point may lie on any existing arrangement feature.
 * \param arr The arrangement.
 * \param p The point to be inserted.
 * \return A handle to the vertex that corresponds to the given point.
 */
template <typename GeomTraits, typename TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Vertex_handle
insert_point(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
             const typename GeomTraits::Point_2& p);

/*!
 * Remove a vertex from the arrangement.
 * \param arr The arrangement.
 * \param v The vertex to remove.
 * \return Whether the vertex has been removed or not.
 */
template <typename GeomTraits, typename TopTraits>
bool
remove_vertex(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
              typename Arrangement_on_surface_2<GeomTraits,
                                                TopTraits>::Vertex_handle v);


/*!
 * Check the validity of the arrangement. In particular, check that the
 * edegs are disjoint-interior, and the holes are located in their proper
 * position.
 * \param arr The arrangement.
 * \return Whether the arrangement is valid.
 */
template <typename GeomTraits, typename TopTraits>
bool is_valid(const Arrangement_on_surface_2<GeomTraits, TopTraits>& arr);

/*!
 * Compute the zone of the given x-monotone curve in the existing arrangement.
 * Meaning, it output the arrangment's vertices, edges and faces that the
 * x-monotone curve intersects.
 * \param arr The arrangement.
 * \param c The x-monotone curve that its zone was computed.
 * \param oi Output iterator of CGAL::Object to insert the zone elements to.
 * \param pi The point location strategy that is used to locate the starting
 * point.
 * \return The output iterator that the curves were inserted to.
 */
template <typename GeomTraits, typename TopTraits,
          typename OutputIterator, typename PointLocation>
OutputIterator zone(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                    const typename GeomTraits::X_monotone_curve_2& c,
                    OutputIterator oi,
                    const PointLocation& pl);

/*!
 * Compute the zone of the given x-monotone curve in the existing arrangement.
 * Overloaded version with no point location object - the walk point-location
 * strategy is used as default.
 * \param arr The arrangement.
 * \param c The x-monotone curve that its zone was computed.
 * \param oi Output iterator of CGAL::Object to insert the zone elements to.
 * \return The output iterator that the curves were inserted to.
 */
template <typename GeomTraits, typename TopTraits, typename OutputIterator>
OutputIterator zone(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                    const typename GeomTraits::X_monotone_curve_2& c,
                    OutputIterator oi);

/*!
 * Checks if the given curve/x-monotone curve intersects the existing
 * arrangement.
 * \param arr The arrangement.
 * \param c The curve/x-monotone curve.
 * \param pi The point location strategy that is used to locate the starting
 * point.
 * \return True if the curve intersect the arrangement, false otherwise.
 */
template <typename GeomTraits, typename TopTraits, typename Curve,
          typename PointLocation>
bool do_intersect(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                  const Curve& c, const PointLocation& pl);

/*!
 * Checks if the given curve/x-monotone curve intersects the existing
 * arrangement.
 * Overloaded version with no point location object - the walk point-location
 * strategy is used as default.
 * \param arr The arrangement.
 * \param c The x-monotone curve/curve.
 * \return True if the curve intersect the arrangement, false otherwise.
 */
template <typename GeomTraits, typename TopTraits, typename Curve>
bool do_intersect(Arrangement_on_surface_2<GeomTraits, TopTraits>& arr,
                  const Curve& c);

} //namespace CGAL

// The function definitions can be found under:
#include <CGAL/Arrangement_2/Arrangement_on_surface_2_impl.h>
#include <CGAL/Arrangement_2/Arrangement_on_surface_2_global.h>

#endif
