// Copyright (c) 2015  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Frank Da, David Cohen-Steiner, Andreas Fabri

#ifndef CGAL_ADVANCING_FRONT_SURFACE_RECONSTRUCTION_VERTEX_BASE_WITH_ID_3_H
#define CGAL_ADVANCING_FRONT_SURFACE_RECONSTRUCTION_VERTEX_BASE_WITH_ID_3_H

#include <CGAL/license/Advancing_front_surface_reconstruction.h>


#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/utility.h>

#include <utility>

#include <list>



namespace CGAL {

  template <class B, class C> class Advancing_front_surface_reconstruction;

  /*!
  \ingroup PkgAdvancingFrontSurfaceReconstruction

  The class `Advancing_front_surface_reconstruction_vertex_base_3` is the default
  vertex type for the class  `Advancing_front_surface_reconstruction`.

  \tparam Traits has to be a model of `DelaunayTriangulationTraits_3`.

  \tparam Vb has to be a model of `TriangulationVertexBase_3`.
  */
  template <typename Traits, typename Vb = Triangulation_vertex_base_3<Traits> >
  class Advancing_front_surface_reconstruction_vertex_base_3 : public Vb
  {
  public:

    template < typename TDS2 >
    struct Rebind_TDS {
      typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
      typedef Advancing_front_surface_reconstruction_vertex_base_3<Traits,Vb2>  Other;
    };


    template <class B,class C> friend class Advancing_front_surface_reconstruction;


    typedef Vb Base;
    typedef typename Base::Vertex_handle Vertex_handle;
    typedef typename Base::Cell_handle Cell_handle;
    typedef typename Vb::Point Point;
    typedef double coord_type;

    typedef Triple< Cell_handle, int, int > Edge;
    typedef std::pair< Edge, int > Edge_incident_facet;
    typedef std::pair< Edge_incident_facet, Edge_incident_facet > IO_edge_type;

    typedef coord_type criteria;

    typedef std::pair< criteria, IO_edge_type > Radius_edge_type;
    typedef std::pair< Radius_edge_type, int > Border_elt;
    typedef std::pair< Vertex_handle, Border_elt > Next_border_elt;


    //par convention je remplis d'abord first et si necessaire second...
    typedef std::pair< Next_border_elt*,  Next_border_elt*> Intern_successors_type;



  public:

    typedef std::pair< criteria, IO_edge_type* > Radius_ptr_type;
    typedef std::pair< Vertex_handle, Vertex_handle > Edge_like;
    typedef std::pair< criteria,  Edge_like > Incidence_request_elt;
    typedef std::list< Incidence_request_elt > Incidence_request_type;
    typedef typename Incidence_request_type::iterator Incidence_request_iterator;


    //-------------------- DATA MEMBERS ---------------------------------

    typedef std::size_t Info;  // so that we are a model of TriangulationVertexBaseWithInfo_3

  private:
    std::size_t m_id;
    int m_mark;
    int m_post_mark;
    Intern_successors_type* m_incident_border;

    // Instead of having a set per vertex, there is a global list
    // in the surface reconstruction class
    // and two iterators per vertex in this list
    // Note that m_ie_last is not past the end
    // m_ie_first == m_ie_last == interior_edge.end()  iff  the set is empty
    typename std::list<Vertex_handle>::iterator m_ie_first, m_ie_last;


    // We do the same for the incidence requests
    typename std::list< Incidence_request_elt >::iterator m_ir_first, m_ir_last;


    //-------------------- CONSTRUCTORS ---------------------------------

  public:

    Advancing_front_surface_reconstruction_vertex_base_3()
      : Vb(), m_mark(-1),
        m_post_mark(-1)
    {}

    Advancing_front_surface_reconstruction_vertex_base_3(const Point & p)
      : Vb(p), m_mark(-1),
        m_post_mark(-1)
    {}

    Advancing_front_surface_reconstruction_vertex_base_3(const Point & p, Cell_handle f)
      : Vb(p, f), m_mark(-1),
        m_post_mark(-1)
    {}

    Advancing_front_surface_reconstruction_vertex_base_3(Cell_handle f)
      : Vb(f), m_mark(-1),
        m_post_mark(-1)
    {}

    Advancing_front_surface_reconstruction_vertex_base_3(const Advancing_front_surface_reconstruction_vertex_base_3& other)
      : Vb(other), m_mark(-1),
        m_post_mark(-1)
    {}

    //-------------------- MEMBER FUNCTIONS -----------------------------

  public:

    std::size_t& id()
    {
      return m_id;
    }

    const std::size_t& id() const
    {
      return m_id;
    }

    std::size_t& info()
    {
      return m_id;
    }

    const std::size_t& info() const
    {
      return m_id;
    }


    //-------------------------------------------------------------------
  private:

    void delete_border()
    {
      m_incident_border = NULL;
    }


    inline Next_border_elt* next_on_border(const int& i) const
    {
      if (m_incident_border == NULL) return NULL; //vh is interior
      if (m_incident_border->first->first != NULL)
	if (m_incident_border->first->second.second == i)
	  return m_incident_border->first;
      if (m_incident_border->second->first != NULL)
	if (m_incident_border->second->second.second == i)
	  return m_incident_border->second;
      return NULL;
    }




    inline bool is_border_edge(Vertex_handle v) const
    {
      if (m_incident_border == NULL) return false;
      return ((m_incident_border->first->first == v)||
	      (m_incident_border->second->first == v));
    }

    inline Next_border_elt* border_elt(Vertex_handle v) const
    {
      if (m_incident_border == NULL) return NULL;
      if (m_incident_border->first->first == v) return m_incident_border->first;
      if (m_incident_border->second->first == v) return m_incident_border->second;
      return NULL;
    }

  public:
    inline Next_border_elt* first_incident() const
    {
      if (m_incident_border == NULL) return NULL;
      return m_incident_border->first;
    }
  private:
    inline Next_border_elt* second_incident() const
    {
      if (m_incident_border == NULL) return NULL;
      return m_incident_border->second;
    }


    inline  void set_next_border_elt(const Next_border_elt& elt)
    {
      if (m_incident_border->first->first == NULL)
	*m_incident_border->first = elt;
      else
	{
	  if (m_incident_border->second->first != NULL)
	    std::cerr << "+++probleme de MAJ du bord <Vertex_base>" << std::endl;
	  *m_incident_border->second = elt;
	}
    }



    //-------------------------------------------------------------------

  public:

    inline bool is_on_border() const
    {
      return (m_mark > 0);
    }

    inline bool not_interior() const
    {
      return (m_mark != 0);
    }

    inline int number_of_incident_border() const
    {
      return m_mark;
    }

    inline bool is_exterior() const
    {
      return (m_mark < 0);
    }

    //-------------------------------------------------------------------
  private:

    inline void inc_mark()
    {
      if (m_mark==-1)
	m_mark=1;
      else
	m_mark++;
    }

    //-------------------------------------------------------------------
  public:
    inline void set_post_mark(const int& i)
    {
      m_post_mark = i;
    }

    inline bool is_post_marked(const int& i)
    {
      return (m_post_mark == i);
    }
  };

} // namespace CGAL

#endif // CGAL_ADVANCING_FRONT_SURFACE_RECONSTRUCTION_VERTEX_BASE_3_H

