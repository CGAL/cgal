#ifndef CGAL_MESH_3_REFINE_TETS_H
#define CGAL_MESH_3_REFINE_TETS_H

#include <CGAL/Mesher_level.h>
#include <CGAL/Mesh_3/Triangulation_mesher_level_traits_3.h>

#include <CGAL/Mesh_3/Refine_edges.h>
#include <CGAL/Mesh_3/Refine_facets.h>
#include <CGAL/Mesh_3/Double_map_container.h>

#include <list>
#include <utility>

namespace CGAL {
namespace Mesh_3 {

template <class Tr,
          class Container = 
          Double_map_container<typename Tr::Cell_handle,
                               typename Tr::typename Geom_traits::Quality> >
>
class Refine_tets_base
{
  typedef typename Tr::Bare_point Bare_point;
  typedef typename Tr::Point Point;
  typedef typename Tr::Weighted_point Weighted_point;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;

  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Triangulation_mesher_level_traits_3<Tr>::Zone Zone;

  typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Tr::Facet_circulator Facet_circulator;

  typedef typename Tr::Facet Facet;

  typedef typename Geom_traits::Quality Quality;

public:
  /** \name CONSTRUCTORS */

  Refine_tets_base(Tr& t) : tr(t) {}

private:
  /* --- private datas --- */
  Tr& tr; /**< The triangulation itself. */
  Container tets_to_refine/**< The set of tets to refine. */


  bool should_be_refined(const Cell_handle c, Quality& q) const
  {
    const Bare_point
      & a = c->vertex(0)->point().point(),
      & b = c->vertex(1)->point().point(),
      & c = c->vertex(2)->point().point();
      & d = c->vertex(3)->point().point();

    return geom_traits().is_bad_object()(a,b,c,d,q);

  }

  bool should_be_refined(const Cell_handle c) const
  {
    Quality q;
    return should_be_refined(c, q);
  }
}

public:
  /** \name Functions that this level must declare. */

  Tr& get_triangulation_ref()
  {
    return tr;
  }

  const Tr& get_triangulation_ref() const
  {
    return tr;
  }

  void do_scan_triangulation()
  {
    for(Cell_iterator cit = tr.finite_cells_begin();
        cit != tr.finite_cells_end();
        ++cit)
      test_cell(cit);
  }

  bool test_cell(const Cell_handle c)
  {
    Quality q;
    if( cit->is_is_domain() && should_be_refined(cit, q) )
      tets_to_refine.add_element(cit, q);
  }

  bool is_no_longer_element_to_refine() const 
  {
    return tets_to_refine.empty();
  }

  Facet do_get_next_element()
  {
    return tets_to_refine.get_next_element();
  }

  void do_pop_next_element()
  {
    tets_to_refine.remove_next_element();
  }

  /** @todo Handle point().point() and Weighted_point. */
  Weighted_point get_refinement_point(const Cell_handle& c) const
  {
    typename Geom_traits::Construct_circumcenter_3 circumcenter = 
      tr.geom_traits().construct_circumcenter_3_object();

    const Bare_point
      & a = c->vertex(0)->point().point(),
      & b = c->vertex(1)->point().point(),
      & c = c->vertex(2)->point().point();
      & d = c->vertex(3)->point().point();

    return Weighted_point(circumcenter(a, b, c, d);
  }

  void do_before_conflicts(const Facet&, const Point&)
  {}

  std::pair<bool, bool> do_private_test_point_conflict(const Point&,
                                               Zone& )
  {
    return std::make_pair(true, true);
  }

  std::pair<bool, bool> do_test_point_conflict_from_superior(const Point&,
                                                             Zone& )
  {
    return std::make_pair(true, true);
  }

  void do_before_insertion(const Facet&, const Point&,
                           const Zone&)
  {
  }

  void do_after_insertion(const Vertex_handle& v)
  {
    // scan tets
    typedef std::list<Cell_handle> Cells;
    typedef typename Cells::iterator Cell_iterator;
    Cells incident_cells;

    tr.incident_cells(v, std::back_inserter(incident_cells));

    for(Cell_iterator it = incident_cells.begin();
        it != incident_cells.end();
        ++it)
      test_cell(cit);
  }

}; // end Refine_tets_base  

  namespace tets {

    template <typename Tr, typename Refine_tets>
    class Refine_edges_visitor {
      Refine_tets* refine_tets;

    public:
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Cell_handle Cell_handle;
      typedef typename Tr::Facet Facet;
      typedef typename Tr::Cell_iterator Cell_iterator;
      typedef ::CGAL::Mesh_3::details::Refine_edges_base_types<Tr> Types;
      typedef typename Types::Constrained_edge Constrained_edge;
      typedef ::CGAL::Triangulation_mesher_level_traits_3<Tr> Traits;
      typedef typename Traits::Zone Zone;
      typedef typename Traits::Point Point;

      Refine_edges_visitor(Refine_tets* refine_tets_)
        : refine_tets(refine_tets_) {}

      void before_insertion(const Constrained_edge&,
                            const Point&,
                            Zone& zone) const 
      {
        for(typename Zone::Facets_iterator fit = zone.internal_facets.begin();
            fit != zone.internal_facets.end();
            ++fit)
          refine_tets->try_to_removed_element(*fit);
      }

      void after_insertion(const Vertex_handle& v) const
      {
        // scan edges 
        typedef std::list<Cell_handle> Cells;
        typedef typename Cells::iterator Cell_iterator;
        Cells incident_cells;
        
        refine_tets->triangulation().
          incident_cells(v, std::back_inserter(incident_cells));

        for(Cell_iterator it = incident_cells.begin();
            it != incident_cells.end();
            ++it)
          refine_tets->test_for_facet(Facet(*it, (*it)->index(v)));
      }

      Null_mesh_visitor previous_level() const 
      {
        return Null_mesh_visitor();
      }
      
    }; // end Refine_edges_visitor

    template <class Tr, typename Refine_tets>
    class Refine_tets_visitor {
      Refine_tets* refine_tets;

    public:
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Facet Facet;
      typedef typename Tr::Cell_iterator Cell_iterator;
      typedef ::CGAL::Triangulation_mesher_level_traits_3<Tr> Traits;
      typedef typename Traits::Zone Zone;
      typedef typename Traits::Point Point;

      typedef Refine_edges_visitor<Tr, Refine_tets> Previous_visitor;

      Refine_tets_visitor(Refine_tets* refine_tets_)
        : refine_tets(refine_tets_) {}

      void before_insertion(const Facet&,
                            const Point&,
                            Zone&) const 
      {
      }

      void after_insertion(const Vertex_handle&) const
      {
      }

      Previous_visitor previous_level() const 
      {
        return Previous_visitor(refine_tets);
      }
    }; // end Refine_tets_visitor

  }; // end namespace Mesh_3::tets
  
template <typename Tr,
          typename Base = Refine_tets_base<Tr>,
          typename Edges_level = Refine_edges<Tr>
 >
class Refine_tets : 
  public Base, 
  public Mesher_level <
    Triangulation_mesher_level_traits_3<Tr>,
    Refine_tets<Tr, Base, Edges_level>,
    typename Tr::Facet,
    Edges_level
  >
{
public:
  typedef Refine_tets<Tr, Base> Self;
  typedef Mesher_level <
    Triangulation_mesher_level_traits_3<Tr>,
    Refine_tets<Tr, Base, Edges_level>,
    typename Tr::Facet,
    Edges_level
  > Mesher;
  Refine_tets(Tr& t, Edges_level* edges_level)
    : Base(t), Mesher(edges_level) {}
};

}; // end namespace Mesh_3
}; // end namespace CGAL

#endif // CGAL_MESH_3_REFINE_TETS_H
