#ifndef CGAL_MESH_3_REFINE_EDGES_H
#define CGAL_MESH_3_REFINE_EDGES_H

#include <CGAL/Mesher_level.h>
#include <CGAL/Mesh_3/Triangulation_mesher_level_traits_3.h>
#include <CGAL/Mesh_3/Simple_queue_container.h>

#include <list>
#include <utility>

namespace CGAL {
  
namespace Mesh_3 {

  namespace details {
    template <typename Tr>
    struct Refine_edges_base_types
    {
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef std::pair<Vertex_handle,
                        Vertex_handle>   Constrained_edge;

      typedef ::CGAL::Mesh_3::Simple_queue_container<Constrained_edge>
                                         Default_container;
    };
  }; // end namespace details

template <
  class Tr,
  class Container = 
    typename details::Refine_edges_base_types<Tr>::Default_container
>
class Refine_edges_base
{
protected:
  typedef typename Tr::Bare_point Bare_point;
  typedef typename Tr::Point Point;
  typedef typename Tr::Weighted_point Weighted_point;
  typedef typename Tr::Edge Edge;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Cell_handle Cell_handle;

  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Triangulation_mesher_level_traits_3<Tr>::Zone Zone;

  typedef typename Tr::Facet_circulator Facet_circulator;
  typedef typename Tr::Finite_edges_iterator Finite_edges_iterator;

  typedef details::Refine_edges_base_types<Tr> Types;

  typedef typename Types::Constrained_edge Constrained_edge;

public:
  /** \name CONSTRUCTORS */

  Refine_edges_base(Tr& t) : tr(t) {}

protected:
  /* --- protected datas --- */
  Tr& tr; /**< The triangulation itself. */
  Container edges_to_be_conformed; /**< Edge queue */
  Vertex_handle va, vb; /**< (\c va, \c vb) is the edge that will be
                           refined. */

  bool is_encroached(Edge e) const
  {
    const Vertex_handle& va = e.first->vertex(e.second);
    const Vertex_handle& vb = e.first->vertex(e.third);
    
    Facet_circulator f_circ = tr.incident_facets(e), end(f_circ);

    do {

      // Small algorithm to retrieve the third vertex of a facet
      // (ch, index), knowing two vertices (Vertex_handle) of the
      // facets.
      // Another implementation would be:
      //   int ia = c->index(va);
      //   int ib = c->index(vb);
      //   int k = 6-index-ia-ib;
      // 
      Vertex_handle v;
      for(int i = 0; i < 4; ++i)
	{
	  if( i == (*f_circ).second ) continue; 
	  v = (*f_circ).first->vertex(i);
	  if( v != va && v != vb) break;
	}
      // here v is the third point of the facet

      typename Geom_traits::Angle_3 angle = 
	tr.geom_traits().angle_3_object();
      
      /** @todo S'occuper de ces point().point() */
      if( angle(va->point().point(), v->point().point(),
		vb->point().point()) == OBTUSE )
	return true;

      ++f_circ;
    }
    while( f_circ != end );
    return false;
  }

  bool is_encroached(Edge e, const Bare_point& p) const 
  {
    const Vertex_handle& va = e.first->vertex(e.second);
    const Vertex_handle& vb = e.first->vertex(e.third);
    return is_encroached(va, vb, p);
  }
  
  bool is_encroached(const Vertex_handle& va, const Vertex_handle& vb, 
		     const Bare_point& p) const 
  {
    typename Geom_traits::Angle_3 angle = 
      tr.geom_traits().angle_3_object();
    
    return( angle(va->point().point(), p,
		  vb->point().point()) == OBTUSE );
  }

  bool test_if_encroached(const Vertex_handle va, const Vertex_handle vb)
  {
    Cell_handle c;
    int i, j;

    CGAL_assertion_code(bool b =) tr.is_edge(va, vb, c, i, j);
    CGAL_assertion(b);

    if( is_encroached(Edge(c, i, j)) )
      {
	edges_to_be_conformed.add_element(std::make_pair(va, vb));
	return true;
      }
    else
      return false;
  }

  void fill_edges_to_be_conformed()
  {
    for(Finite_edges_iterator it = tr.finite_edges_begin();
	it != tr.finite_edges_end();
	++it)
      if(complex_2.is_in_complex(*it) && is_encroached(*it))
	edges_to_be_conformed.add_element(std::make_pair(it->first->
						  vertex(it->second),
						  it->first->
						  vertex(it->third)));
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
    fill_edges_to_be_conformed();
  }

  bool is_no_longer_element_to_refine() const 
  {
    return edges_to_be_conformed.empty();
  }

  Constrained_edge do_get_next_element()
  {
    return edges_to_be_conformed.get_next_element();
  }

  void do_pop_next_element()
  {
    edges_to_be_conformed.remove_next_element();
  }

  /** @todo Handle point().point() and Weighted_point. */
  Weighted_point get_refinement_point(const Constrained_edge& edge) const
  {
    typename Geom_traits::Construct_midpoint_3 midpoint = 
      tr.geom_traits().construct_midpoint_3_object();

    return Weighted_point(midpoint(edge.first->point().point(),
                                   edge.second->point().point()));
  }

  void do_before_conflicts(const Constrained_edge&, const Point&)
  {}

  std::pair<bool, bool>
  do_test_point_conflict_from_superior(const Point& p,
                                       Zone& z)
  {
    bool one_edge_is_encroached = false;

    for(typename Zone::Facets_iterator fit = zone.boundary_facets.begin();
        fit != zone.boundary_facets.end();
        ++fit)
      {
        const Cell_handle& c = (*fit).first;
        const int index = (*fit).second;

	// @todo Utiliser la table vertex_triple_index et boucler
	// entre 0 et 2 au lieu de faire la boucle suivante.

        // 0 < i_off < j_off <4
        for(int i_off = 1; i_off < 3; ++i_off)
          for(int j_off = i_off+1; j_off < 4; ++j_off)
            {
              const int i = (index+i_off) % 4;
              const int j = (index+j_off) % 4;
              // (i,j) is one of the three edges of the facet *fit
              const Vertex_handle& vi = c->vertex(i);
              const Vertex_handle& vj = c->vertex(j);
              if( vi < vj ) // (vi, vj) will be visited twice.
		// BUG!!!
		if( is_not_locally_Delaunay(vi, vj) )
                  one_edge_is_encroached = true;
            }
      }
    // The first boolean tells if there is a conflict.
    // The second one tells if the element can be dropped.
    return std::make_pair(!one_edge_is_encroached, false);
    // Here, can_drop=false, because the facet (from which the point
    // is the circumcenter) should be refined latter.
  }

  std::pair<bool, bool> do_private_test_point_conflict(const Point&,
                                                       Zone& ) const
  {
    // nothing
    return std::make_pair(true, true);
  }

  void do_before_insertion(const Constrained_edge& e, const Point&,
                           const Zone&)
  {
    va = e.first;
    vb = e.first;
  }

  void do_after_insertion(const Vertex_handle& v)
  {
    test_if_encroached(va, v);
    test_if_encroached(vb, v);

    // scan edges 
    typedef std::list<Cell_handle> Cells;
    typedef typename Cells::iterator Cell_iterator;
    Cells incident_cells;

    tr.incident_cells(v, std::back_inserter(incident_cells));

    for(Cell_iterator it = incident_cells.begin();
        it != incident_cells.end();
        ++it)
      {
        const int index = (*it)->index(v);
        // for (i,j) from 0 to 3, i<j, i!=index, j!=index
        for(int i = 0; i<3; ++i)
          {
            if(i == index) continue;
            for(int j = i; j<4; ++j)
              {
                if(j == index) continue;
                const int k = 6-index-i-j; // the last vertex
                // Consider the edge (i, j).
                // This edge will be twice in the zone. So one must test if
                // the orientation (index, i, j, k) is positive. If not,
                // the other occurence will be positively oriented.
                typename Geom_traits::Orientation_3 orientation = 
                  tr.geom_traits().orientation_3_object();
                const Vertex_handle& vindex = (*it)->vertex(index);
                const Vertex_handle& vi = (*it)->vertex(i);
                const Vertex_handle& vj = (*it)->vertex(j);
                const Vertex_handle& vk = (*it)->vertex(k);
                if( orientation(vindex->point().point(),
                                vi->point().point(),
                                vj->point().point(),
                                vk->point().point())
                    == CGAL::POSITIVE )
                  if(complex_2.is_in_complex(*it, i, j))
                    test_if_encroached(vi, vj);
              }
          }
      }
  }

}; // end Refine_edges_base  

template <typename Tr,
          typename Base = Refine_edges_base<Tr> >
struct Refine_edges : 
  public Base, 
  public Mesher_level <
    Triangulation_mesher_level_traits_3<Tr>,
    Refine_edges<Tr, Base>,
    std::pair<typename Tr::Vertex_handle,
              typename Tr::Vertex_handle>,
    Null_mesher_level
  >
{
  typedef Refine_edges<Tr, Base> Self;
  typedef Mesher_level <
    Triangulation_mesher_level_traits_3<Tr>,
    Refine_edges<Tr, Base>,
    std::pair<typename Tr::Vertex_handle,
              typename Tr::Vertex_handle>,
    Null_mesher_level
  > Mesher;

  Refine_edges(Tr& t): Base(t), Mesher(&null_mesher_level) {}
}; // end Refine_edges

}; // end namespace Mesh_3
}; // end namespace CGAL

#endif // CGAL_MESH_3_REFINE_EDGES_H
