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

#ifndef CGAL_MESH_2_REFINE_EDGES_H
#define CGAL_MESH_2_REFINE_EDGES_H

#include <CGAL/Mesher_level.h>
#include <CGAL/Mesh_2/Triangulation_mesher_level_traits_2.h>
#include <CGAL/Mesh_2/Filtered_queue_container.h>

#include <utility>

namespace CGAL {

/**
 * \namespace Mesh_2
 *   Defines classes that are not yet documented.
 * 
 * \namespace Mesh_2::details
 *   Namespace for internal use.
 */

namespace Mesh_2 {

  namespace details {

    /** This class defines several auxiliary types for \c Refine_edges. */
    template <typename Tr>
    struct Refine_edges_base_types
    {
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Face_handle Face_handle;

      typedef std::pair<Vertex_handle,
                        Vertex_handle> Constrained_edge;

      /** Object predicate that tests if a given \c Constrained_Edge is
          really an edge of the triangulation and is constrained.
      */
      class Is_really_a_constrained_edge {
        const Tr& tr;
      public:
        /** \param tr_ points to the triangulation. */
        explicit Is_really_a_constrained_edge(const Tr& tr_) : tr(tr_) {}

        bool operator()(const Constrained_edge& ce) const
        {
          Face_handle fh;
          int i;
          return tr.is_edge(ce.first, ce.second, fh,i) &&
            fh->is_constrained(i);
        }
      };

      typedef ::CGAL::Mesh_2::Filtered_queue_container<Constrained_edge,
                                               Is_really_a_constrained_edge>
                           Default_container;
    };

  }; // end namespace details


  /**
   * Predicate class that verifies that an edge is locally conforming
   * Gabriel. Moreover, This classes defines a predicate that test if an
   * edge is encroached by a given point.
   * \param Tr The type of the trianglation.
   */
  template <typename Tr>
  struct Is_locally_conforming_Gabriel
  {
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Face_handle Face_handle;
    typedef typename Tr::Point Point;
    typedef typename Tr::Geom_traits Geom_traits;

    /** Operator that takes an edge (\c fh, \c index). */
    bool operator()(Tr& ct,
                    const Face_handle& fh,
                    const int i) const
    {
//         if( ct.extras().is_bad( static_cast<const Tr&>(ct),
//                                 fh, i ) ) return false;

      typedef typename Geom_traits::Angle_2 Angle_2;
      
      const Angle_2 angle = ct.geom_traits().angle_2_object();

      const Vertex_handle& va = fh->vertex(ct. cw(i));
      const Vertex_handle& vb = fh->vertex(ct.ccw(i));      

      const Point& a = va->point();
      const Point& b = vb->point();

      const Vertex_handle& vi = fh->vertex(i);
      const Vertex_handle& mvi = fh->mirror_vertex(i);

      return( ( ct.is_infinite(vi) || 
                angle(a, vi->point(), b) != OBTUSE)
              &&
              ( ct.is_infinite(mvi) || 
                angle(a, mvi->point(), b) != OBTUSE)
              );
    }

    /** Operator that takes an edge (\c va, \c vb). */
    bool operator()(Tr& ct,
                    const Vertex_handle& va,
                    const Vertex_handle& vb) const
    {
      Face_handle fh;
      int i;
      CGAL_assertion_code( bool should_be_true = )
      ct.is_edge(va, vb, fh, i);
      CGAL_assertion( should_be_true );
      
      return this->operator()(ct, fh, i);
    }

    /**
     * Operator that takes an edge (\c fh, \c index) and a point \c p.
     * Tests if the point encroached the edge.
     */
    bool operator()(Tr& ct,
                    const Face_handle& fh,
                    const int i,
                    const Point& p) const
    {
      return this->operator()(ct,
                              fh->vertex(ct. cw(i)),
                              fh->vertex(ct.ccw(i)),
                              p);
    }

    /**
     * Operator that takes an edge (\c va, \c vb) and a point \c p.
     * Tests if the point encroached the edge.
     */
    bool operator()(Tr& ct,
                    const Vertex_handle& va,
                    const Vertex_handle& vb,
                    const Point& p) const
      {
//         if( ct.extras().is_bad( static_cast<const Tr&>(ct),
//                                  fh, i ) ) return false;

        typedef typename Geom_traits::Angle_2 Angle_2;

        const Angle_2 angle = ct.geom_traits().angle_2_object();

        const Point& a = va->point();
        const Point& b = vb->point();

        return( angle(a, p, b) != OBTUSE );
      }
  };

  /**
   * Predicate class that verifies that an edge is locally conforming
   * Delaunay.
   * \param Tr The type of the trianglation.
   */
  template <typename Tr>
  struct Is_locally_conforming_Delaunay
  {
    typedef typename Tr::Vertex_handle Vertex_handle;
    typedef typename Tr::Face_handle Face_handle;
    typedef typename Tr::Point Point;
    typedef typename Tr::Geom_traits Geom_traits;

    /** Operator that takes an edge (\c fh, \c index). */
    bool operator()(Tr& ct,
                    const Face_handle& fh,
                    const int i) const
    {
//       if( ct.extras().is_bad( static_cast<const Tr&>(ct),
//                                fh, i ) ) return false;

      typedef typename Geom_traits::Side_of_oriented_circle_2
        Side_of_oriented_circle_2;

      Side_of_oriented_circle_2 in_circle =
        ct.geom_traits().side_of_oriented_circle_2_object();
      
      const Vertex_handle& vi = fh->vertex(i);
      const Vertex_handle& mvi = fh->mirror_vertex(i);

      if(ct.is_infinite(vi) || ct.is_infinite(mvi)){
        return true;
      }

      const Point& a = fh->vertex(ct. cw(i))->point();
      const Point& b = fh->vertex(ct.ccw(i))->point();
      const Point& c = vi->point();
      const Point& d = mvi->point();

      return( in_circle(c, b, a, d) == ON_NEGATIVE_SIDE );
    }

    /** Operator that takes an edge (\c va, \c vb). */
    bool operator()(Tr& ct,
                    const Vertex_handle& va,
                    const Vertex_handle& vb) const
    {
//       if( ct.extras().is_bad( static_cast<const Tr&>(ct),
//                                fh, i ) ) return false;

      typedef typename Geom_traits::Side_of_oriented_circle_2
        Side_of_oriented_circle_2;

      Side_of_oriented_circle_2 in_circle =
        ct.geom_traits().side_of_oriented_circle_2_object();

      Face_handle fh;
      int i;
      CGAL_assertion_code( bool test = )
        ct.is_edge(va, vb, fh, i);
      CGAL_assertion( test == true );

      const Vertex_handle& vi = fh->vertex(i);
      const Vertex_handle& mvi = fh->mirror_vertex(i);

      if(ct.is_infinite(vi) || ct.is_infinite(mvi)){
        return true;
      }

      const Point& a = va->point();
      const Point& b = vb->point();
      const Point& c = vi->point();
      const Point& d = mvi->point();

      return( in_circle(c, b, a, d) == ON_NEGATIVE_SIDE );
    }
  };

/**
 * This class is the base for the first level of Mesh_2: the edge
 * conforming level. It does not handle clusters.
 *
 * \param Tr is the type of triangulation on which the level acts.
 * \param Is_locally_conform defines the locally conform criterion: Gabriel
 *        or Delaunay. It defaults to the Garbriel criterion.
 * \param Container is the type of container. It defaults to a filtered
 *        queue of \c Vertex_handle pair (see \c Filtered_queue_container).
 */
template <
  class Tr,
  class Is_locally_conform = Is_locally_conforming_Gabriel<Tr>,
  class Container = 
    typename details::Refine_edges_base_types<Tr>::Default_container
>
class Refine_edges_base
{
  typedef typename Tr::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Tr::Face_circulator Face_circulator;
  
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Face_handle Face_handle;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Point Point;
  typedef typename Tr::Geom_traits Geom_traits;

  typedef typename Triangulation_mesher_level_traits_2<Tr>::Zone Zone;

  typedef typename details::Refine_edges_base_types<Tr>::Constrained_edge
                               Constrained_edge;

protected:
  /* --- protected datas --- */

  Tr& tr; /**< The triangulation itself. */

  /** Predicates to filter edges. */
  const typename details::Refine_edges_base_types<Tr>::
    Is_really_a_constrained_edge is_really_a_constrained_edge;

  Container edges_to_be_conformed; /**< Edge queue */

  /** The object predicate that defines the locally conform criteria. */
  Is_locally_conform is_locally_conform;

  Vertex_handle va, vb;

public:
  /** \name CONSTRUCTORS */

  Refine_edges_base(Tr& tr_) :
    tr(tr_), is_really_a_constrained_edge(tr_),
    edges_to_be_conformed(is_really_a_constrained_edge),
    is_locally_conform()
  {
  }

  /** \name HELPING FUNCTIONS */

  void clear()
  {
    edges_to_be_conformed.clear();
  }

  /** \name Functions that this level must declare. */

  Tr& get_triangulation_ref()
  {
    return tr;
  }

  const Tr& get_triangulation_ref() const
  {
    return tr;
  }

  /** Scans all constrained edges and put them in the queue if they are
      encroached. */
  void do_scan_triangulation()
  {
    clear();
#ifndef CGAL_IT_IS_A_CONSTRAINED_TRIANGULATION_PLUS
    for(Finite_edges_iterator ei = tr.finite_edges_begin();
        ei != tr.finite_edges_end();
        ++ei)
      if(ei->first->is_constrained(ei->second) &&
         !is_locally_conform(tr, ei->first, ei->second) )
        add_constrained_edge_to_be_conformed(*ei);
    
#else
  for(typename Tr::Subconstraint_iterator it = tr.subconstraints_begin();
      it != tr.subconstraints_end(); ++it) 
    {
      const Vertex_handle& va = it->first.first;
      const Vertex_handle& vb = it->first.second;
      
      if(fh->is_constrained(i) &&
         !is_locally_conform(tr, va, vb) )
        add_constrained_edge_to_be_conformed(va, vb);
    }
#endif
  } // end do_scan_triangulation()

  /** Tells if the queue of edges to be conformed is empty or not. */
  bool is_no_longer_element_to_refine()
  {
    return edges_to_be_conformed.empty();
  }

  /** Get the next edge to conform. */
  Constrained_edge do_get_next_element()
  {
    return edges_to_be_conformed.get_next_element();
  }

  /** Pop the first edge of the queue. */
  void do_pop_next_element()
  {
    edges_to_be_conformed.remove_next_element();
  }

  /** This version computes the refinement point without handling
      clusters. The refinement point of an edge is just the middle point of
      the segment. */
  Point get_refinement_point(const Constrained_edge& edge) const
  {
    typename Geom_traits::Construct_midpoint_2
      midpoint = tr.geom_traits().construct_midpoint_2_object();

    const Vertex_handle& va = edge.first;
    const Vertex_handle& vb = edge.second;

    return midpoint(va->point(), vb->point());
  }

  /** Unmark as constrained. */
  void do_before_conflicts(const Constrained_edge& e, const Point&)
  {
    /*const Vertex_handle& */va = e.first;
    /*const Vertex_handle& */vb = e.second;
    
    Face_handle fh;
    int i;
    CGAL_assertion_code( bool should_be_true= )
    tr.is_edge(va, vb, fh, i);
    CGAL_assertion( should_be_true );
    std::cerr <<"ok\n";
    
    tr.remove_constrained_edge(fh, i);
    tr.is_edge(va, vb);
  }

  void do_after_no_insertion(const Constrained_edge& e, const Point&,
                             const Zone& )
  {
    const Vertex_handle& va = e.first;
    const Vertex_handle& vb = e.second;
    
    tr.insert_constraint(va, vb);
  }

  /**
   * Test if the edges of the boundary are locally conforming.
   * Push which that are not in the list of edges to be conformed.
   */
  std::pair<bool, bool>
  do_test_point_conflict_from_superior(const Point& p,
                                       Zone& z)
  {
    bool no_edge_is_encroached = true;

    for(typename Zone::Edges_iterator eit = z.boundary_edges.begin();
        eit != z.boundary_edges.end(); ++eit)
      { 
        const Face_handle& fh = it->first;
        const int& i = it->second;

        if(fh->is_constrained(i) && !is_locally_conform(tr, fh, i, p))
          {
            add_constrained_edge_to_be_conformed(*eit);
            no_edge_is_encroached = false;
          }
      }
    return std::make_pair(no_edge_is_encroached, true);
  }

  /** Do nothing */
  std::pair<bool, bool>
  do_private_test_point_conflict(const Point&, Zone& ) const
  {
    CGAL_assertion(tr.is_edge(va, vb));
    return std::make_pair(true, true);
  }

  /** Saves the handles of the edge that will be splitted. */
  void do_before_insertion(const Constrained_edge& e, const Point&,
                           const Zone&)
  {
    CGAL_assertion(tr.is_edge(va, vb));
    va = e.first;
    vb = e.second;
    CGAL_assertion(tr.is_edge(va, vb));
    std::cerr << "OK ok\n";
  }

  /**
   * Scans the edges of the star boundary, to test if they are both
   * locally conforming. If not, push them in the list of edges to be
   * conformed.
   */
  void do_after_insertion(const Vertex_handle& v)
  {
    // @todo Perhaps we should remove destroyed edges too.
    // @warning This code has been rewroten!

    Face_circulator fc = tr.incident_faces(v), fcbegin(fc);
    if( fc == 0 ) return;

    do {
      const int i = fc->index(v);
      CGAL_assertion( i>=0 && i < 4);
      if( fc->is_constrained(i) &&
          !is_locally_conform(tr, fc, i) )
        add_constrained_edge_to_be_conformed(Edge(fc, i));
      ++fc;
    } while( fc != fcbegin );

  if(!is_locally_conform(tr, va, v))
    add_constrained_edge_to_be_conformed(va, v);
  
  if(!is_locally_conform(tr, vb, v))
    add_constrained_edge_to_be_conformed(vb, v);
  } // end do_after_insertion

protected:
  /** \name Auxiliary functions */

  /** Add an \c Edge \c e in the queue. */
  void add_constrained_edge_to_be_conformed(const Edge& e)
  {
    const Vertex_handle& va = e.first->vertex(tr. cw(e.second));
    const Vertex_handle& vb = e.first->vertex(tr.ccw(e.second));
    edges_to_be_conformed.add_element(std::make_pair(va, vb));
  }

  /** Add an edge (\c va,\c  vb) in the queue. */
  void add_constrained_edge_to_be_conformed(const Vertex_handle& va,
                                            const Vertex_handle& vb)
  {
    edges_to_be_conformed.add_element(std::make_pair(va, vb));
  }

}; // end class Refine_edges_base

  namespace details {
    template <typename Tr, typename Self>
    struct Refine_edges_types
    {
      typedef Mesher_level <
        Triangulation_mesher_level_traits_2<Tr>,
        Self,
        typename ::CGAL::Mesh_2::details::Refine_edges_base_types<Tr>
        ::Constrained_edge,
        Null_mesher_level > Edges_mesher_level;
    }; // end Refine_edges_types
  } // end namespace details

template <typename Tr,
          typename Is_locally_conform = Is_locally_conforming_Gabriel<Tr>,
          typename Base = Refine_edges_base<Tr, Is_locally_conform> >
struct Refine_edges : 
  public Base, 
  public details::Refine_edges_types<Tr, 
    Refine_edges<Tr, Is_locally_conform, Base> >::Edges_mesher_level
{
  typedef Refine_edges<Tr, Is_locally_conform, Base> Self;
  typedef typename details::Refine_edges_types<Tr, Self>::Edges_mesher_level
                                 Mesher;
public:
  Refine_edges(Tr& t, Null_mesher_level& null_level)
    : Base(t), Mesher(null_level)
  {
  }
}; // end Refine_edges


} // end namespace Mesh_2

} // end namespace CGAL

#endif // CGAL_MESH_2_REFINE_EDGES_H
