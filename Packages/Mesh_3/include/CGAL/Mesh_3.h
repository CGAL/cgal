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
class Mesh_3 : public Constrained_regular_triangulation_3<Tr> {
public:
  typedef Tr Triangulation;
  typedef Constrained_regular_triangulation_3<Tr> CDT;
  typedef Mesh_3<Triangulation> Self;
  
  typedef typename Tr::Geom_traits Geom_traits;
  typedef typename Tr::Triangulation_data_structure Tds;

  typedef typename Tr::Cell_handle Cell_handle;

  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;

  typedef typename CDT::Constrained_edge Constrained_edge;

  typedef typename CDT::Geom_traits Geom_traits;
  typedef typename Geom_traits::Bare_point Point_3;

  typedef std::pair<Point_3, Point_3> Constraint_2;
  typedef std::list<Constraint_3> List_of_constraints_3;
  typedef std::list<Point> List_of_seeds;

public:
  // constructor from a Triangulation& t
  explicit 
  Mesh_3(const Geom_traits& gt = Geom_traits()) 
    : CDT(gt) {};

  struct Geom_traits_2 : public Geom_traits_2_3, public Geom_traits
  {
    typedef typename Geom_traits_2_3::Compare_x_2 Compare_x_2;
    typedef typename Geom_traits_2_3::Compare_y_2 Compare_y_2;
    typedef typename Geom_traits_2_3::Orientation_2 Orientation_2;
    typedef typename Geom_traits_2_3::Point_2 Point_2;
    typedef typename Geom_traits_2_3::Segment_2 Segment_2;
    typedef typename Geom_traits_2_3::Triangle_2 Triangle_2;
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

}; // end Mesh_3

} //end namespace CGAL

#endif
