#ifndef CGAL_MESH_3_H
#define CGAL_MESH_3_H

#include <CGAL/basic.h>

namespace CGAL {

// Tr: a conform Delaunay triangulation
template <class Tr>
class Mesh_3 : public Tr {
public:
  typedef Tr Triangulation;
  typedef Mesh_3<Triangulation> Self;
  
  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Tr::Triangulation_data_structure Tds;

  typedef typename Tr::Face_handle            Face_handle;
public:
    // constructor from a Triangulation& t
  explicit 
  Mesh_3(Triangulation& t, const Geom_traits& gt = Geom_traits(),
	 bool dont_refine = false);

  void fill_facette_map();

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
  
private:
  
  edges_to_be_conform 

}; // end Mesh_3

} //end namespace CGAL

#endif
