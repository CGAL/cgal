// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
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
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_MESH_2_REFINE_FACES_H
#define CGAL_MESH_2_REFINE_FACES_H

#include <CGAL/Double_map.h>
#include <boost/iterator/transform_iterator.hpp>

namespace CGAL {

namespace Mesh_2 {

template <typename Tr, typename Criteria>
class Refine_faces_base
{
  /** \name Types from Tr. */

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

protected: // --- PROTECTED TYPES ---
  /** Meshing criteria. */
  typedef typename Criteria::Is_bad Is_bad;
  typedef typename Criteria::Quality Quality;

  /** \name typedefs for private members types */

  typedef CGAL::Double_map<Face_handle, Quality> Bad_faces;

protected:
  // --- PROTECTED MEMBER DATAS ---

  Tr& tr; /**< The triangulation itself. */
  Triangulation_mesher_level_traits_2<Tr> traits;
  Criteria& criteria; /**<The meshing criteria */

  /** List of bad finite faces */
  Bad_faces bad_faces;
public:
  /** \name CONSTRUCTORS */

  Refine_faces_base(Tr& t, Criteria& criteria_) 
    : tr(t), criteria(criteria_) 
  {
  }

  /** \Name MESHER_LEVEL FUNCTIONS */

  Tr& get_triangulation_ref()
  {
    return tr;
  }

  const Tr& get_triangulation_ref() const
  {
    return tr;
  }

  Triangulation_traits& get_triangulation_traits()
  {
    return traits;
  }

  const Triangulation_traits& get_triangulation_traits() const
  {
    return traits;
  }

  /** Scans all marked faces and put them in the map if they are
      bad. */
  void do_scan_triangulation()
  {
    for(typename Tr::Finite_faces_iterator fit = tr.finite_faces_begin();
        fit != tr.finite_faces_end();
        ++fit)
    {
      Quality q;
      if( fit->is_marked() && is_bad(fit, q) )
        push_in_bad_faces(fit, q);
    }
  }

  /** Tells if the map of faces to be conformed is empty or not. */
  bool is_no_longer_element_to_refine() const
  {
    return bad_faces.empty();
  }

  /** Get the next face to conform. */
  Face_handle do_get_next_element()
  {
    Face_handle fh = bad_faces.front()->second;

    CGAL_assertion_code(typename Geom_traits::Orientation_2 orientation =
                        tr.geom_traits().orientation_2_object());
    CGAL_assertion(orientation(fh->vertex(0)->point(),
                               fh->vertex(1)->point(),
                               fh->vertex(2)->point()) != COLLINEAR );
    return fh;
  }

  /** Pop the first face of the map. */
  void do_pop_next_element()
  {
    bad_faces.pop_front();
  }

  /** Returns the circumcenter of the face. */
  Point get_refinement_point(const Face_handle& f) const
  {
    return tr.circumcenter(f);
  }

  /** Do nothing */
  void do_before_conflicts(const Face_handle&, const Point&)
  {
  }

  /** Do nothing */
  std::pair<bool, bool>
  do_test_point_conflict_from_superior(const Point&,
                                       Zone&)
  {
    return std::make_pair(true, true);
  }

  /** Do nothing */
  std::pair<bool, bool>
  do_private_test_point_conflict(const Point&, Zone& ) const
  {
    return std::make_pair(true, true);
  }

  /** Remove the conflicting faces from the bad faces map. */
  void do_before_insertion(const Face_handle& fh, const Point&,
                           Zone& zone)
  {
    /** @todo Perhaps this function is useless. */
    for(typename Zone::Faces_iterator fh_it = zone.faces.begin();
        fh_it != zone.faces.end();
        ++fh_it)
      {
        if(*fh_it != fh && (*fh_it)->is_marked() )
          remove_bad_face(*fh_it);
        (*fh_it)->set_marked(false);
      }
  }

  /** Do nothing. */
  void do_after_no_insertion(const Face_handle&, const Point&,
                             Zone&)
  {
  }

  /** Restore markers in the star of \c v. */
  void do_after_insertion(const Vertex_handle& v)
  {
    typename Tr::Face_circulator fc = tr.incident_faces(v), fcbegin(fc);
    do {
      fc->set_marked(true);
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

  bool is_bad(const Face_handle fh, Quality& q) const;
  bool is_bad(const Face_handle fh) const ;

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
    for(Fh_it pfit=begin; pfit!=end; ++pfit)
      push_in_bad_faces(*pfit, Quality());
  }

}; // end class Refine_faces_base
  
// --- PRIVATE MEMBER FUNCTIONS ---

template <typename Tr, typename Criteria>
inline
void Refine_faces_base<Tr, Criteria>::
push_in_bad_faces(Face_handle fh, const Quality& q)
{
#ifdef DEBUG
  std::cerr << "push_in_bad_faces("
            << fh->vertex(0)->point() << ","
            << fh->vertex(1)->point() << ","
            << fh->vertex(2)->point() << ")\n";
#endif // DEBUG
  CGAL_assertion_code(typename Geom_traits::Orientation_2 orientation =
                      tr.geom_traits().orientation_2_object());
  CGAL_assertion( orientation(fh->vertex(0)->point(),
                              fh->vertex(1)->point(),
                              fh->vertex(2)->point()) != COLLINEAR );
  CGAL_assertion(fh->is_marked());
  bad_faces.insert(fh, q);
}

template <typename Tr, typename Criteria>
inline
void Refine_faces_base<Tr, Criteria>::
remove_bad_face(Face_handle fh)
{
#ifdef DEBUG
  std::cerr << "bad_faces.erase("
            << fh->vertex(0)->point() << ","
            << fh->vertex(1)->point() << ","
            << fh->vertex(2)->point() << ")\n";
#endif // DEBUG
  bad_faces.erase(fh);
}

template <typename Tr, typename Criteria>
void Refine_faces_base<Tr, Criteria>::
compute_new_bad_faces(Vertex_handle v)
{
  typename Tr::Face_circulator fc = v->incident_faces(), fcbegin(fc);
  do {
    Quality q;
    if(!tr.is_infinite(fc))
      if( fc->is_marked() && is_bad(fc, q) )
        push_in_bad_faces(fc, q);
    fc++;
  } while(fc!=fcbegin);
}

template <typename Tr, typename Criteria>
inline
bool Refine_faces_base<Tr, Criteria>::
is_bad(const Face_handle f, Quality& q) const
{
  return criteria.is_bad_object()(f, q);
}

template <typename Tr, typename Criteria>
inline
bool Refine_faces_base<Tr, Criteria>::
is_bad(const Face_handle f) const
{
  Quality q;
  return criteria.is_bad(f, q);
}

  namespace details {
    template <typename Tr, typename Self, typename Previous>
    struct Refine_faces_types
    {
      typedef Mesher_level <
        Triangulation_mesher_level_traits_2<Tr>,
        Self,
        typename Tr::Face_handle,
        Previous > Faces_mesher_level;
    }; // end Refine_faces_types
  } // end namespace details

template <typename Tr,
          typename Criteria,
          typename Previous,
          typename Base = Refine_faces_base<Tr, Criteria> >
class Refine_faces : 
  public Base, 
  public details::Refine_faces_types<Tr, 
    Refine_faces<Tr, Criteria, Previous, Base>,
    Previous>::Faces_mesher_level
{
  typedef typename Tr::Geom_traits Geom_traits;

  template <class Pair>
  struct Pair_get_first: public std::unary_function<Pair,
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

  typedef typename Base::Bad_faces Bad_faces;

  typedef typename boost::transform_iterator<
    Pair_get_first<typename Bad_faces::Direct_entry>,
    typename Bad_faces::const_iterator>
   Bad_faces_const_iterator;

public:
  Refine_faces(Tr& t, Criteria& criteria, Previous& previous)
    : Base(t, criteria), Mesher(previous)
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
    CGAL_assertion_code(typename Geom_traits::Orientation_2 orientation =
                        this->tr.geom_traits().orientation_2_object());
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
          this->bad_faces.dump_direct_func(std::cerr);

          return false;
        }
    return true;  
  }

}; // end Refine_faces

} // end namespace Mesh_2

} // end namespace CGAL

#endif // CGAL_MESH_2_REFINE_FACES_H
