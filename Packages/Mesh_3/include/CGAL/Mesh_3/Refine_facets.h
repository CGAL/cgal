#ifndef CGAL_MESH_3_REFINE_FACETS_H
#define CGAL_MESH_3_REFINE_FACETS_H

#include <CGAL/Mesher_level.h>
#include <CGAL/Mesh_3/Triangulation_mesher_level_traits_3.h>
#include <CGAL/Mesh_3/Complex_2_in_triangulation_3.h>

#include <CGAL/Mesh_3/Refine_edges.h>
#include <CGAL/Mesh_3/Simple_set_container.h>

#include <list>
#include <utility>

namespace CGAL {
namespace Mesh_3 {

template <class Tr,
          class Container = Simple_set_container<typename Tr::Facet>
>
class Refine_facets_base
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

  typedef Complex_2_in_triangulation_3<Tr, void, void> Complex_2;

  typedef Container Facets_to_be_conformed;
public:
  /** \name CONSTRUCTORS */

  Refine_facets_base(Tr& t, Complex_2& comp) 
    : tr(t), complex_2(comp)
  {
  }

private:
  /* --- private datas --- */
  Tr& tr; /**< The triangulation itself. */
  Complex_2& complex_2; /**< The Complex_2_in_triangulation_3 */
  Facets_to_be_conformed facets_to_be_conformed; /**< The set of facets to
                                                    refine. */
  void fill_facets_to_be_conformed()
  {
    for(Finite_facets_iterator it = tr.finite_facets_begin();
	it != tr.finite_facets_end();
	++it)
      if(complex_2.is_in_complex(*it) && is_not_locally_Delaunay(*it))
	facets_to_be_conformed.add_element(*it);
  }

  bool is_not_locally_Delaunay(Facet f)
  {
    typename Geom_traits::Power_test_3 power_test = 
      tr.geom_traits().power_test_3_object();

    const Cell_handle& c = f.first;
    const int i = f.second;
    const Cell_handle& n = c->neighbor(i);
    const Vertex_handle& va = c->vertex( (i+1) % 4 );
    const Vertex_handle& vb = c->vertex( (i+2) % 4 );
    const Vertex_handle& vc = c->vertex( (i+3) % 4 );

    const Vertex_handle& v1 = c->vertex(i);
    const Vertex_handle& v2 = n->vertex( n->index(c) );

    return( power_test(va->point(),
		       vb->point(),
		       vc->point(),
		       v1->point(),
		       v2->point()) != ON_NEGATIVE_SIDE ); // À VÉRIFIER !
  }

public:
  bool test_for_facet(const Facet& facet)
  {
    if(complex_2.is_in_complex(facet) &&
       is_not_locally_Delaunay(facet))
      {
        facets_to_be_conformed.add_element(facet);
        return true;
      }
    return false;
  }

  void try_to_removed_element(Facet& facet)
  {
    facets_to_be_conformed.remove_element(facet);
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

  void do_scan_triangulation()
  {
    fill_facets_to_be_conformed();
  }

  bool is_no_longer_element_to_refine() const 
  {
    return facets_to_be_conformed.empty();
  }

  Facet do_get_next_element()
  {
    return facets_to_be_conformed.get_next_element();
  }

  void do_pop_next_element()
  {
    facets_to_be_conformed.remove_next_element();
  }

  /** @todo Handle point().point() and Weighted_point. */
  Weighted_point get_refinement_point(const Facet& facet) const
  {
    typename Geom_traits::Construct_circumcenter_3 circumcenter = 
      tr.geom_traits().construct_circumcenter_3_object();

    const int& index = facet.second; // aliases index

    const int i = (index+1) % 4;
    const int j = (index+2) % 4;
    const int k = (index+3) % 4;

    const Vertex_handle& vi = facet.first->vertex(i);
    const Vertex_handle& vj = facet.first->vertex(j);
    const Vertex_handle& vk = facet.first->vertex(k);
    return Weighted_point(circumcenter(vi->point().point(),
                                       vj->point().point(),
                                       vk->point().point()));
  }

  void do_before_conflicts(const Facet&, const Point&)
  {}

  std::pair<bool, bool> do_private_test_point_conflict(const Point&,
                                               Zone& )
  {
    return std::make_pair(true, true);
  }

  std::pair<bool, bool> do_test_point_conflict_from_superior(const Point& p,
                                                             Zone& zone)
  {
    typename Geom_traits::Power_test_3 power_test = 
      tr.geom_traits().power_test_3_object();
    
    bool one_facet_is_encroached = false;

    for(typename Zone::Facets_iterator fit = zone.boundary_facets.begin();
        fit != zone.boundary_facets.end();
        ++fit)
      {
        const Cell_handle& c = (*fit).first;
        const int& i = (*fit).second; // aliases i

        const Vertex_handle& va = c->vertex( (i+1) % 4 );
        const Vertex_handle& vb = c->vertex( (i+2) % 4 );
        const Vertex_handle& vc = c->vertex( (i+3) % 4 );
        
        if( power_test(va->point(),
		       vb->point(),
		       vc->point(),
                       p) != ON_NEGATIVE_SIDE )
          {
            one_facet_is_encroached = true;
            facets_to_be_conformed.add_element(*fit);
          }
      }
     return std::make_pair(!one_facet_is_encroached, true);
  }

  void do_before_insertion(const Facet&, const Point&,
                           const Zone&)
  {
  }

  void do_after_insertion(const Vertex_handle& v)
  {
    // scan facets
    typedef std::list<Cell_handle> Cells;
    typedef typename Cells::iterator Cell_iterator;
    Cells incident_cells;

    tr.incident_cells(v, std::back_inserter(incident_cells));

    for(Cell_iterator it = incident_cells.begin();
        it != incident_cells.end();
        ++it)
      test_for_facet(Facet(*it, (*it)->index(v)));


  }

}; // end Refine_edges_base  

  namespace facets {

    template <typename Tr, typename Refine_facets>
    class Refine_edges_visitor {
      Refine_facets* refine_facets;

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

      Refine_edges_visitor(Refine_facets* refine_facets_)
        : refine_facets(refine_facets_) {}

      void before_insertion(const Constrained_edge&,
                            const Point&,
                            Zone& zone) const 
      {
        for(typename Zone::Facets_iterator fit = zone.internal_facets.begin();
            fit != zone.internal_facets.end();
            ++fit)
          refine_facets->try_to_removed_element(*fit);
      }

      void after_insertion(const Vertex_handle& v) const
      {
        // scan edges 
        typedef std::list<Cell_handle> Cells;
        typedef typename Cells::iterator Cell_iterator;
        Cells incident_cells;
        
        refine_facets->triangulation().
          incident_cells(v, std::back_inserter(incident_cells));

        for(Cell_iterator it = incident_cells.begin();
            it != incident_cells.end();
            ++it)
          refine_facets->test_for_facet(Facet(*it, (*it)->index(v)));
      }

      Null_mesh_visitor previous_level() const 
      {
        return Null_mesh_visitor();
      }
      
    }; // end Refine_edges_visitor

    template <class Tr, typename Refine_facets>
    class Refine_facets_visitor {
      Refine_facets* refine_facets;

    public:
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Tr::Facet Facet;
      typedef typename Tr::Cell_iterator Cell_iterator;
      typedef ::CGAL::Triangulation_mesher_level_traits_3<Tr> Traits;
      typedef typename Traits::Zone Zone;
      typedef typename Traits::Point Point;

      typedef Refine_edges_visitor<Tr, Refine_facets> Previous_visitor;

      Refine_facets_visitor(Refine_facets* refine_facets_)
        : refine_facets(refine_facets_) {}

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
        return Previous_visitor(refine_facets);
      }
    }; // end Refine_facets_visitor

  }; // end namespace Mesh_3::facets
  
template <typename Tr,
          typename Base = Refine_facets_base<Tr>,
          typename Edges_level = Refine_edges<Tr>
 >
class Refine_facets : 
  public Base, 
  public Mesher_level <
    Triangulation_mesher_level_traits_3<Tr>,
    Refine_facets<Tr, Base, Edges_level>,
    typename Tr::Facet,
    Edges_level
  >
{
public:
  typedef Refine_facets<Tr, Base> Self;
  typedef Mesher_level <
    Triangulation_mesher_level_traits_3<Tr>,
    Refine_facets<Tr, Base, Edges_level>,
    typename Tr::Facet,
    Edges_level
  > Mesher;
  Refine_facets(Tr& t, Edges_level* edges_level)
    : Base(t), Mesher(edges_level) {}
};

}; // end namespace Mesh_3
}; // end namespace CGAL

#endif // CGAL_MESH_3_REFINE_FACETS_H
