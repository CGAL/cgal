// Copyright (c) 2004-2006  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_2_REFINE_FACES_H
#define CGAL_MESH_2_REFINE_FACES_H

#include <CGAL/license/Mesh_2.h>


#include <CGAL/Mesh_2/Face_badness.h>
#include <CGAL/Double_map.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>

#include <string>
#include <sstream>

namespace CGAL {

namespace Mesh_2 {

// Previous is the whole previous edges_level.
template <typename Tr, typename Criteria, typename Previous>
class Refine_faces_base :
    public Triangulation_mesher_level_traits_2<Tr>,
    public No_test_point_conflict,
    public No_after_no_insertion
{
  /** \name Types from Tr. */

  typedef Tr Triangulation;

  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Geom_traits::FT FT;
  typedef FT      Squared_length;

  typedef typename Tr::Vertex_handle        Vertex_handle;
  typedef typename Tr::Face_handle          Face_handle;

  typedef typename Tr::Face_circulator        Face_circulator;
  typedef typename Tr::Finite_faces_iterator  Finite_faces_iterator;
  typedef typename Tr::All_faces_iterator     All_faces_iterator;
  typedef typename Tr::Point                  Point;

  typedef Triangulation_mesher_level_traits_2<Tr> Triangulation_traits;
  typedef typename Triangulation_traits::Zone Zone;

public:
  using Triangulation_mesher_level_traits_2<Tr>::triangulation_ref_impl;

protected: // --- PROTECTED TYPES ---
  /** Meshing criteria. */
  typedef typename Criteria::Is_bad s_bad;
  typedef typename Criteria::Quality Quality;

  /** \name typedefs for private members types */

  struct Face_compare {
    bool operator()(const Face_handle& fh1, const Face_handle& fh2) const {
      if(fh1->vertex(0)->point() < fh2->vertex(0)->point())
        return true;
      else if(fh1->vertex(0)->point() == fh2->vertex(0)->point()) {
        if(fh1->vertex(1)->point() < fh2->vertex(1)->point())
          return true;
        else if(fh1->vertex(1)->point() == fh2->vertex(1)->point() &&
                fh1->vertex(2)->point() < fh2->vertex(2)->point())
          return true;
      }
      return false;
    }
  };
  
  typedef CGAL::Double_map<Face_handle, Quality, Face_compare> Bad_faces;

protected:
  // --- PROTECTED MEMBER DATAS ---

  Criteria& criteria; /**<The meshing criteria */
  Previous& previous;

  /** List of bad finite faces */
  Bad_faces bad_faces;

  Mesh_2::Face_badness current_badness;

public:
  /** \name CONSTRUCTORS */

  Refine_faces_base(Tr& t, Criteria& criteria_, Previous& prev) 
    : Triangulation_traits(t), criteria(criteria_), previous(prev)
  {
  }

  /** \Name MESHER_LEVEL FUNCTIONS */

  /** Scans all faces in domain and put them in the map if they are
      bad. */
  void scan_triangulation_impl()
  {
    bad_faces.clear();
#ifdef CGAL_MESH_2_DEBUG_BAD_FACES
    std::cerr << "bad_faces.clear()\n";
#endif // CGAL_MESH_2_DEBUG_BAD_FACES

    for(typename Tr::Finite_faces_iterator fit =
	  triangulation_ref_impl().finite_faces_begin();
        fit != triangulation_ref_impl().finite_faces_end();
        ++fit)
    {
      if( fit->is_in_domain() )
	{
	  Quality q;
	  Mesh_2::Face_badness badness = is_bad(fit, q);
	  if( badness != Mesh_2::NOT_BAD )
	    push_in_bad_faces(fit, q);
	}
    }
  }

  Zone conflicts_zone_impl(const Point& p, Face_handle fh)
  {
    Zone zone;

    zone.fh = triangulation_ref_impl().locate(p, zone.locate_type, zone.i, fh);

    zone.parent_face = fh;

    triangulation_ref_impl().
      get_conflicts_and_boundary(p,
                                 std::back_inserter(zone.faces),
                                 std::back_inserter(zone.boundary_edges),
				 fh
                                 );
#ifdef CGAL_MESH_2_DEBUG_CONFLICTS_ZONE
    std::cerr << "get_conflicts_and_boundary(" << p << "):" << std::endl
              << "faces: " << zone.faces.size() << std::endl
              << "edges: " << zone.boundary_edges.size() << std::endl;
#endif // CGAL_MESH_2_DEBUG_CONFLICTS_ZONE
    return zone;
  }

  /** Tells if the map of faces to be conformed is empty or not. */
  bool no_longer_element_to_refine_impl() const
  {
    return bad_faces.empty();
  }

  /** Get the next face to conform. */
  Face_handle get_next_element_impl()
  {
    Face_handle fh = bad_faces.front()->second;
    current_badness = is_bad(bad_faces.front()->first);

    CGAL_assertion_code
      (typename Geom_traits::Orientation_2 orientation =
       triangulation_ref_impl().geom_traits().orientation_2_object()
       );
    CGAL_assertion(orientation(fh->vertex(0)->point(),
                               fh->vertex(1)->point(),
                               fh->vertex(2)->point()) != COLLINEAR );
    return fh;
  }

  /** Pop the first face of the map. */
  void pop_next_element_impl()
  {
    bad_faces.pop_front();
  }

  /** Returns the circumcenter of the face. */
  Point refinement_point_impl(const Face_handle& f) const
  {
    return triangulation_ref_impl().circumcenter(f);
  }

  /** \todo ?? */
  void before_conflicts_impl(const Face_handle&, const Point&)
  { /// @todo modularize
    previous.set_imperative_refinement(current_badness == 
				       Mesh_2::IMPERATIVELY_BAD);
  }

  /** Remove the conflicting faces from the bad faces map. */
  void before_insertion_impl(const Face_handle&, const Point&,
			     Zone& zone)
  {
    /** @todo Perhaps this function is useless. */
    for(typename Zone::Faces_iterator fh_it = zone.faces.begin();
        fh_it != zone.faces.end();
        ++fh_it)
      {
        if((*fh_it)->is_in_domain() )
          remove_bad_face(*fh_it);
        (*fh_it)->set_in_domain(false);
      }
  }

  /** Restore markers in the star of \c v. */
  void after_insertion_impl(const Vertex_handle& v)
  {
#ifdef CGAL_MESH_2_VERBOSE
    std::cerr << "*";
#endif
    typename Tr::Face_circulator fc = 
      triangulation_ref_impl().incident_faces(v), fcbegin(fc);
    do {
      fc->set_in_domain(true);
    } while (++fc != fcbegin);
    compute_new_bad_faces(v);
  }

private:
  /** \name AUXILIARY FUNCTIONS */

  /** Auxiliary function called to put a new face in the map. */
  void push_in_bad_faces(Face_handle fh, const Quality& q);


public:
  /** \name Functions that maintain the map of bad faces. */

  /**
   * Updates the map with faces incident to the vertex \a v.
   * @todo The visitor should be made friend, instead of this function to
   * be public.
   */
  void compute_new_bad_faces(Vertex_handle v);

  /** Auxiliary function called to erase a face handle from the map. */
  void remove_bad_face(Face_handle fh);

public:
  /** \name ACCESS FUNCTION */

  Mesh_2::Face_badness is_bad(const Face_handle fh, Quality& q) const;
  Mesh_2::Face_badness is_bad(const Face_handle fh) const;
  Mesh_2::Face_badness is_bad(Quality q) const;

  /**
   * Adds the sequence [\c begin, \c end[ to the list
   * of bad faces.
   * Use this overriden function if the list of bad faces can be
   * computed easily without testing all faces.
   * \param Fh_it is an iterator of \c Face_Handle.
   */
  template <class Fh_it>
  void set_bad_faces(Fh_it begin, Fh_it end)
  {
    bad_faces.clear();
#ifdef CGAL_MESH_2_DEBUG_BAD_FACES
    std::cerr << "bad_faces.clear()\n";
#endif // CGAL_MESH_2_DEBUG_BAD_FACES
    for(Fh_it pfit=begin; pfit!=end; ++pfit)
      push_in_bad_faces(*pfit, Quality());
  }

}; // end class Refine_faces_base
  
// --- PRIVATE MEMBER FUNCTIONS ---

template <typename Tr, typename Criteria, typename Previous>
inline
void Refine_faces_base<Tr, Criteria, Previous>::
push_in_bad_faces(Face_handle fh, const Quality& q)
{
#ifdef CGAL_MESH_2_DEBUG_BAD_FACES
  std::cerr << "push_in_bad_faces("
            << fh->vertex(0)->point() << ","
            << fh->vertex(1)->point() << ","
            << fh->vertex(2)->point() << ")\n";
#endif // CGAL_MESH_2_DEBUG_BAD_FACES
  CGAL_assertion_code
    (typename Geom_traits::Orientation_2 orientation =
     triangulation_ref_impl().geom_traits().orientation_2_object()
     );
  CGAL_assertion( orientation(fh->vertex(0)->point(),
                              fh->vertex(1)->point(),
                              fh->vertex(2)->point()) != COLLINEAR );
  CGAL_assertion(fh->is_in_domain());
  bad_faces.insert(fh, q);
}

template <typename Tr, typename Criteria, typename Previous>
inline
void Refine_faces_base<Tr, Criteria, Previous>::
remove_bad_face(Face_handle fh)
{
#ifdef CGAL_MESH_2_DEBUG_BAD_FACES
  std::cerr << "bad_faces.erase("
            << fh->vertex(0)->point() << ","
            << fh->vertex(1)->point() << ","
            << fh->vertex(2)->point() << ")\n";
#endif // CGAL_MESH_2_DEBUG_BAD_FACES
  bad_faces.erase(fh);
}

template <typename Tr, typename Criteria, typename Previous>
void Refine_faces_base<Tr, Criteria, Previous>::
compute_new_bad_faces(Vertex_handle v)
{
  typename Tr::Face_circulator fc = triangulation_ref_impl().incident_faces(v), fcbegin(fc);
  do {
    if(!triangulation_ref_impl().is_infinite(fc))
      if( fc->is_in_domain() )
	{
	  Quality q;
	  Mesh_2::Face_badness badness = is_bad(fc, q);
	  if( badness != Mesh_2::NOT_BAD )
	    push_in_bad_faces(fc, q);
	}
    fc++;
  } while(fc!=fcbegin);
}

template <typename Tr, typename Criteria, typename Previous>
inline
Mesh_2::Face_badness
Refine_faces_base<Tr, Criteria, Previous>::
is_bad(const Face_handle f, Quality& q) const
{
  return criteria.is_bad_object()(f, q);
}

template <typename Tr, typename Criteria, typename Previous>
inline
Mesh_2::Face_badness
Refine_faces_base<Tr, Criteria, Previous>::
is_bad(const Face_handle f) const
{
  Quality q;
  return criteria.is_bad_object()(f, q);
}

template <typename Tr, typename Criteria, typename Previous>
inline
Mesh_2::Face_badness
Refine_faces_base<Tr, Criteria, Previous>::
is_bad(Quality q) const
{
  return criteria.is_bad_object()(q);
}

  namespace details {
    template <typename Tr, typename Self, typename Previous>
    struct Refine_faces_types
    {
      typedef Mesher_level <
	Tr,
        Self,
        typename Tr::Face_handle,
        Previous,
	Triangulation_mesher_level_traits_2<Tr>
	>
      Faces_mesher_level;
    }; // end Refine_faces_types
  } // end namespace details

template <typename Tr,
          typename Criteria,
          typename Previous,
          typename Base = Refine_faces_base<Tr, Criteria, Previous> >
class Refine_faces : 
  public Base, 
  public details::Refine_faces_types<Tr, 
    Refine_faces<Tr, Criteria, Previous, Base>,
    Previous>::Faces_mesher_level
{
  typedef typename Tr::Geom_traits Geom_traits;

  template <class Pair>
  struct Pair_get_first: public CGAL::unary_function<Pair,
                                                    typename Pair::first_type>
  {
    typedef typename Pair::first_type result;
    const result& operator()(const Pair& p) const
    {
      return p.first;
    }
  };

public:
  typedef Refine_faces<Tr, Criteria, Previous, Base> Self;
  typedef typename details::Refine_faces_types<Tr, Self, Previous>
    ::Faces_mesher_level Mesher;

  typedef Tr Triangulation;

  typedef typename Base::Bad_faces Bad_faces;

  typedef typename boost::transform_iterator<
    Pair_get_first<typename Bad_faces::Direct_entry>,
    typename Bad_faces::const_iterator>
   Bad_faces_const_iterator;

public:
  Refine_faces(Tr& t, Criteria& criteria, Previous& previous)
    : Base(t, criteria, previous), Mesher(previous)
  {
  }

  /** \name DEBUGGING FUNCTIONS */

  Bad_faces_const_iterator begin() const
  {
    return Bad_faces_const_iterator(this->bad_faces.begin());
  }

  Bad_faces_const_iterator end() const
  {
    return Bad_faces_const_iterator(this->bad_faces.end());
  }

  bool check_bad_faces()
  {
    struct Output_bad_face {
      std::string operator()(typename Tr::Face_handle fh)
      {
	std::stringstream str;
	
	str << "("
	    << fh->vertex(0)->point() << ", "
            << fh->vertex(1)->point() << ", "
            << fh->vertex(2)->point()
	    << ")";
	return str.str();
      }
    };

    CGAL_assertion_code
      (typename Geom_traits::Orientation_2 orientation =
       this->triangulation_ref_impl().geom_traits().orientation_2_object()
       );
    for(Bad_faces_const_iterator fit = begin();
        fit != end();
        ++fit)
      if( orientation((*fit)->vertex(0)->point(),
                      (*fit)->vertex(1)->point(),
                      (*fit)->vertex(2)->point()) == COLLINEAR )
        {
          std::cerr << "collinear("
                    << (*fit)->vertex(0)->point() << ", "
                    << (*fit)->vertex(1)->point() << ", "
                    << (*fit)->vertex(2)->point()<< ") == true"
                    << std::endl;
          std::cerr << "Dump of bad_faces:" << std::endl;
          this->bad_faces.dump_direct_func(std::cerr, Output_bad_face());

          return false;
        }
    return true;  
  }

}; // end Refine_faces

} // end namespace Mesh_2

} // end namespace CGAL

#endif // CGAL_MESH_2_REFINE_FACES_H
