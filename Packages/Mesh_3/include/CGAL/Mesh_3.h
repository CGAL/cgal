#ifndef CGAL_MESH_3_H
#define CGAL_MESH_3_H

#include <CGAL/basic.h>
#include <list>
#include <utility> //std::pair

#include <CGAL/Triangulation_2_traits_3.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Conforming_Delaunay_triangulation_2.h>

namespace CGAL {

// Tr: a conform Delaunay triangulation
template <class Tr>
class Mesh_3 : public Constrained_Delaunay_triangulation_3<Tr> {
public:
  typedef Tr Triangulation;
  typedef Constrained_Delaunay_triangulation_3<Tr> CDT;
  typedef Mesh_3<Triangulation> Self;
  
  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Tr::Triangulation_data_structure Tds;

  typedef typename Tr::Cell_handle Cell_handle;

  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;

  typedef typename CDT::Constrained_edge Constrained_edge;

  typedef typename CDT::Geom_traits Geom_traits;
  typedef typename Geom_traits::Point_3 Point;

  typedef std::pair<Point, Point> Constraint_2;
  typedef std::list<Constraint_2> List_of_constraints_2;
  typedef std::list<Point> List_of_seeds;

public:
  // constructor from a Triangulation& t
  explicit 
  Mesh_3(const Geom_traits& gt = Geom_traits()) 
    : CDT(gt) {};

  void
  create_clusters();

  void
  create_clusters_of_vertex(Vertex_handle v);


  typedef Triangulation_2_traits_3<Geom_traits> Geom_traits_2_3;

  struct Geom_traits_2 : public Geom_traits_2_3, public Geom_traits
  {
    typedef typename Geom_traits_2_3::Compare_x_2 Compare_x_2;
    typedef typename Geom_traits_2_3::Compare_y_2 Compare_y_2;
    typedef typename Geom_traits_2_3::Orientation_2 Orientation_2;
    typedef typename Geom_traits_2_3::Point_2 Point_2;
    typedef typename Geom_traits_2_3::Segment_2 Segment_2;
    typedef typename Geom_traits_2_3::Triangle_2 Triangle_2;
  };
  
  void insert_facet(List_of_constraints_2 contraints,
		    List_of_seeds seeds,
		    bool seeds_mark = false)
  {
    typedef Triangulation_vertex_base_2<Geom_traits> Vb;
    typedef Delaunay_mesh_face_base_2<Geom_traits> Fb;
    typedef Triangulation_data_structure_2<Vb, Fb> Tds;
    typedef Constrained_Delaunay_triangulation_2<Geom_traits_2, Tds,
      CGAL::Exact_predicates_tag> CDT_2;
    typedef Constrained_triangulation_plus_2<CDT_2> CDTP_2;
    typedef Conforming_Delaunay_triangulation_2<CDTP_2> Conf_DT_2;

    typedef typename List_of_constraints_2::iterator Constraints_2_iterator;
    typedef typename List_of_seeds::iterator Seeds_iterator;

    typedef typename Conf_DT_2::Edge_iterator Edge_iterator;

    Conf_DT_2 cdt2;

    for(Constraints_2_iterator it = contraints.begin();
	it != contraints.end();
	++it)
      cdt2.insert_constraint(it->first, it->second);

    cdt2.mark_facets(seeds.begin(), seeds.end(), seeds_mark);

    cdt2.make_conforming_Delaunay();
  };

 private:
   // PRIVATE TYPES
   class Is_really_a_constrained_edge {
     const Self& _m;
   public:
     explicit Is_really_a_constrained_edge(const Self& m) : _m(m) {};
     bool operator()(const Constrained_edge& ce) const
       {
	 return ce.first->is_adjacent_by_constraint(ce.second);
       }
   };

   template <class Cont, class Pred>
   class Filter_insert_iterator {
   protected:
     Cont* container;
     Pred predicate;
   public:
     typedef Cont container_type;
     typedef std::output_iterator_tag iterator_category;
     typedef void value_type;
     typedef void difference_type;
    typedef void pointer;
    typedef void reference;
      
    explicit Filter_insert_iterator(Cont& x, Pred p = Pred()) :
      container(&x), predicate(p) {}

    Filter_insert_iterator<Cont, Pred>&
    operator=(const typename Cont::value_type& value) { 
      if( predicate(value) )
	container->insert(value);
      return *this;
    }

    Filter_insert_iterator<Cont, Pred>& 
    operator*() { return *this; }

    Filter_insert_iterator<Cont, Pred>& 
    operator++() { return *this; }

    Filter_insert_iterator<Cont, Pred>& 
    operator++(int) { return *this; }
  };

  template <class Cont, class Pred>
  inline 
  Filter_insert_iterator<Cont, Pred> 
  filter_inserter(Cont& x, Pred p = Pred()) {
    return Filter_insert_iterator<Cont, Pred>(x, p);
  };

  struct Cluster {
    bool reduced; // Is the cluster reduced

    typedef std::map<Vertex_handle, bool> Vertices_map;
    Vertices_map vertices;

     bool is_reduced() const {
      return reduced;
    }

    bool is_reduced(const Vertex_handle v) {
      return vertices[v];
    }
  };

  typedef std::multimap<Vertex_handle, Cluster> Cluster_map;
private:
  Cluster_map cluster_map;

}; // end Mesh_3

template <class Tr>
void 
Mesh_3<Tr>::
create_clusters()
{
  for(Finite_vertices_iterator vit = finite_vertices_begin();
      vit != finite_vertices_end();
      vit++)
    create_clusters_of_vertex(vit);
}

template <class Tr>
void 
Mesh_3<Tr>::
create_clusters_of_vertex(Vertex_handle v)
{
  typedef std::set<Vertex_handle> Set_of_vertices;
  typedef std::list<Vertex_handle> Clustered_vertices;
  typedef std::list<Clustered_vertices> List_of_baby_clusters;
  typedef typename List_of_baby_clusters::iterator Cl_it;

  List_of_baby_clusters futur_clusters;
  Set_of_vertices vh_set;

  // store incident vertices into vh_set;
  incident_vertices(v, std::inserter(vh_set, vh_set.end()));
}

} //end namespace CGAL

#endif
