#ifndef CGAL_MESHER_LEVEL_H
#define CGAL_MESHER_LEVEL_H

#include <boost/tuple/tuple.hpp>
#include <utility>

namespace CGAL {

  struct Null_mesh_visitor {
    Null_mesh_visitor previous_level() const { return *this; }

    template <typename E, typename P, typename Z>
    void before_insertion(E, P, Z) const {}

    template <typename V>
    void after_insertion(V) const {}
  };

struct Null_mesher_level {

  bool no_longer_element_to_refine() { return true; }

  template <typename Visitor>
  void refine(Visitor) {}
  
  template <typename P, typename Z>
  std::pair<bool, bool> test_point_conflict_from_superior(P, Z)
  { 
    return std::make_pair(true, false);
  }

  template <typename E, typename P>
  void before_conflicts(E, P) {}

  template <typename E, typename P, typename Z, typename Visitor>
  void before_insertion(E, P, Z, Visitor) {}

  template <typename V, typename Visitor>
  void after_insertion(V, Visitor) {}

}; // end Null_mesher_level

template <
  class Triangulation_traits, /**< Traits class that defines operations
                                    from the trianguilation. */
  class Derived, /**< Base class, that implement methods */
  class Element, /**< Type of elements that this level refines. */
  class Previous = Null_mesher_level /**< Previous level type, 
                                        defaults to \c Null_level
                                               */
  > 
class Mesher_level
{
public:
  /** Type of triangulation that is meshed. */
  typedef typename Triangulation_traits::Triangulation Triangulation;
  /** Type of point that are inserted into the triangulation. */
  typedef typename Triangulation_traits::Point Point;
  /** Type of vertex handles that are returns by insertions into the
      triangulation. */
  typedef typename Triangulation_traits::Vertex_handle Vertex_handle;
  /** Type of the conflict zone for a point that can be inserted. */
  typedef typename Triangulation_traits::Zone Zone;

private:
  /** \name Private member functions */
  
  /** Curiously recurring template pattern. */
  //@{
  Derived& derived()
  {
    return static_cast<Derived&>(*this);
  }

  const Derived& derived() const
  {
    return static_cast<const Derived&>(*this);
  }
  //@}

  /** \name Private member datas */

  Previous* previous_level; /**< The previous level of the refinement
                                    process. */
public:
  typedef Previous Previous_level;

  typedef Mesher_level<Triangulation_traits,
                       Derived,
                       Element,
                       Previous_level> Self;

  /** \name CONSTRUCTORS */
  Mesher_level(Previous_level* previous): previous_level(previous) {};

  /** \name FUNCTIONS IMPLEMENTED IN THE CLASS \c Derived */

  /**  Access to the triangulation */
  //@{
  Triangulation& triangulation()
  {
    return derived().get_triangulation_ref();
  }

  const Triangulation& triangulation() const
  {
    return derived().get_triangulation_ref();
  }
  //@}


  /** Called before the first refinement, to initialized the queue of
      elements that should be refined. */
  void scan_triangulation()
  {
    derived().do_scan_triangulation();
  }

  /** Tells if, as regards the elements of type \c Element, the refinement is
      done. */
  bool no_longer_element_to_refine()
  {
    return derived().is_no_longer_element_to_refine();
  }

  /** Retrieves the next element that could be refined. */
  Element get_next_element()
  {
    return derived().do_get_next_element();
  }

  /** Remove from the list the next element that could be refined. */
  void pop_next_element()
  {
    derived().do_pop_next_element();
  }

  /** Gives the point that should be inserted to refine the element \c e */
  Point refinement_point(const Element& e)
  {
    return derived().get_refinement_point(e);
  }

  /** Actions before testing conflicts for point \c p and element \c e */
  void before_conflicts(const Element& e, const Point& p)
  {
    previous_level->before_conflicts(e, p);
    derived().do_before_conflicts(e, p);
  }

  /** Tells if, as regards this level of the refinement process, if the
      point conflicts with something, and do what is needed. The return
      type is made of two booleans:
        - the first one tells if the point can be inserted,
        - in case of, the first one is \c false, the second one tells if
        the tested element should be reconsidered latter.
  */
  std::pair<bool, bool> private_test_point_conflict(const Point& p,
                                                    Zone& zone)
  {
    return derived().do_private_test_point_conflict(p, zone);
  }

  /** Tells if, as regards this level of the refinement process, if the
      point conflicts with something, and do what is needed. The return
      type is made of two booleans:
        - the first one tells if the point can be inserted,
        - in case of, the first one is \c false, the second one tells if
        the tested element should be reconsidered latter.
      This function is called by the superior level, if any.
  */
  std::pair<bool, bool> test_point_conflict_from_superior(const Point& p,
                                                          Zone& zone)
  {
    return derived().do_test_point_conflict_from_superior(p, zone);
  }

  /** Actions before inserting the point \c p in order to refine the
      element \c e. The zone of conflicts is \c zone. */  
  template <class Mesh_visitor>
  void before_insertion(Element& e, const Point& p, Zone& zone, 
                        Mesh_visitor visitor)
  {
    derived().do_before_insertion(e, p, zone);
    visitor.before_insertion(e, p, zone);
  }

  /** Actions after having inserted the point.
      \param vh is the vertex handle of the inserted point,
      \param visitor is the visitor.*/
  template <class Mesh_visitor>
  void after_insertion(Vertex_handle vh, Mesh_visitor visitor)
  {
    visitor.after_insertion(vh);
    derived().do_after_insertion(vh);
    previous_level->after_insertion(vh, visitor.previous_level());
  }

  /** \name MESHING PROCESS */
  template <class Mesh_visitor>
  void refine(Mesh_visitor visitor)
  {
    while(! previous_level->no_longer_element_to_refine() ||
          ! no_longer_element_to_refine() )
    {
      previous_level->refine(visitor.previous_level());
      if(! no_longer_element_to_refine() )
        process_one_element(visitor);
    }
  }

  template <class Mesh_visitor>
  void process_one_element(Mesh_visitor visitor)
  {
    Element e = get_next_element();
    if( try_to_refine_element(e, visitor) )
      pop_next_element();
  }

  template <class Mesh_visitor>
  bool try_to_refine_element(Element e, Mesh_visitor visitor)
  {
    const Point& p = refinement_point(e);

    before_conflicts(e, p);

    Zone zone = Triangulation_traits().get_conflicts_zone(triangulation(), p);

    bool can_refine_element;
    bool drop_element;
    ::boost::tie(can_refine_element, drop_element) =
      test_point_conflict(p, zone);

    if(can_refine_element)
    {
      before_insertion(e, p, zone, visitor);

      Vertex_handle v = Triangulation_traits().insert(triangulation(),
                                                      p, zone);

      after_insertion(v, visitor);
    }

    return drop_element;
  }

  std::pair<bool, bool> test_point_conflict(const Point& p, Zone& zone)
  {
    bool can_refine_element;
    bool drop_element;
    ::boost::tie(can_refine_element, drop_element) = 
      previous_level->test_point_conflict_from_superior(p, zone);

    if( !can_refine_element )
      return std::make_pair(can_refine_element, drop_element);
    return private_test_point_conflict(p, zone);
  }

}; // end Mesher_level

}; // end namespace CGAL

#endif // CGAL_MESHER_LEVEL_H
