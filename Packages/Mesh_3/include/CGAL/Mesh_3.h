#ifndef CGAL_MESH_3_H
#define CGAL_MESH_3_H

#include <CGAL/basic.h>
#include <list>

namespace CGAL {

// Tr: a conform Delaunay triangulation
template <class Tr>
class Mesh_3 : public Tr {
public:
  typedef Tr Triangulation;
  typedef Mesh_3<Triangulation> Self;
  
  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Tr::Triangulation_data_structure Tds;

  typedef typename Tr::Face_handle Face_handle;

  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
public:
    // constructor from a Triangulation& t
  explicit 
  Mesh_3(Triangulation& t, const Geom_traits& gt = Geom_traits(),
	 bool dont_refine = false);

  void create_clusters();
private:
  // PRIVATE TYPES
  class Is_really_a_constrained_edge {
    const Self& _m;
  public:
    explicit Is_really_a_constrained_edge(const Self& m) : _m(m) {};
    bool operator()(const Constrained_edge& ce) const
      {
	Face_handle fh;
	int i;
	return _m.is_edge(ce.first, ce.second, fh,i) &&
	  fh->is_constrained(i);
      }
  };

  
  template <class Cont, class Pred>
  class Filter_insert_iterator {
  protected:
    Cont* container;
    Predicate p;
  public:
    typedef Cont container_type;
    typedef std::output_iterator_tag iterator_category;
    typedef void value_type;
    typedef void difference_type;
    typedef void pointer;
    typedef void reference;
      
    explicit Filter_insert_iterator(Cont& x, Pred p = Pred()) :
      container(&x), predicate(p) {}

    Filter_insert_iterator<Cont>&
    operator=(const typename Cont::value_type& value) { 
      if( predicate(value) )
	container->insert(value);
      return *this;
    }

    Filter_insert_iterator<Cont>& operator*() { return *this; }
    Filter_insert_iterator<Cont>& operator++() { return *this; }
    Filter_insert_iterator<Cont>& operator++(int) { return *this; }
  };

  inline 
  template <class Cont, class Pred>
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
void Mesh_3<Tr>::
create_clusters()
{
  for(Finite_vertices_iterator vit = finite_vertices_begin();
      vit != finite_vertices_end();
      vit++)
    create_clusters_of_vertex(vit);
}

template <class Tr>
void Mesh_2<Tr>::
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
