#ifndef POLYHEDRON_TYPE_H
#define POLYHEDRON_TYPE_H

// CGAL
// kernel
#include "Kernel_type.h"

// surface mesh
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_3.h>

#include <set>

template <typename Refs, typename Tag, typename Point, typename Patch_id>
class Polyhedron_demo_vertex : 
  public CGAL::HalfedgeDS_vertex_base<Refs, Tag, Point>
{
public:
  typedef std::set<Patch_id> Set_of_indices;

private:
  typedef CGAL::HalfedgeDS_vertex_base<Refs, Tag, Point> Pdv_base;

  Set_of_indices indices;
  std::size_t mID;
public:
  int nb_of_feature_edges;

  bool is_corner() const {
    return nb_of_feature_edges > 2;
  }

  bool is_feature_vertex() const {
    return nb_of_feature_edges != 0;
  }

  void add_incident_patch(const Patch_id i) {
    indices.insert(i);
  }

  const Set_of_indices&
  incident_patches_ids_set() const {
    return indices;
  }

  std::size_t& id()       { return mID; }
  std::size_t  id() const { return mID; }
  
  Polyhedron_demo_vertex() : Pdv_base(), mID(-1), nb_of_feature_edges(0) {}
  Polyhedron_demo_vertex(const Point& p) : Pdv_base(p), mID(-1), nb_of_feature_edges(0) {}
};

template <class Refs, class Tprev, class Tvertex, class Tface>
class Polyhedron_demo_halfedge : 
  public CGAL::HalfedgeDS_halfedge_base<Refs,Tprev,Tvertex,Tface>
{
private:
  bool feature_edge;
  std::size_t mID;
public:

  Polyhedron_demo_halfedge() 
    : feature_edge(false), mID(-1) {};

  bool is_feature_edge() const {
    return feature_edge;
  }

  void set_feature_edge(const bool b) {
    feature_edge = b;
    this->opposite()->feature_edge = b;
  }
  
  std::size_t& id()       { return mID; }
  std::size_t  id() const { return mID; }

};

template <class Refs, class T_, class Pln_, class Patch_id_>
class Polyhedron_demo_face : 
  public CGAL::HalfedgeDS_face_base<Refs,T_,Pln_>
{
private:
  Patch_id_ patch_id_;
  std::size_t mID;
public:
  typedef Patch_id_ Patch_id;
  
  Polyhedron_demo_face() 
    : patch_id_(1), mID(-1) {}
  
  int patch_id() const {
    return patch_id_;
  }
  
  void set_patch_id(const int i) {
    patch_id_ = i;
  }
  
  std::size_t& id()       { return mID; }
  std::size_t  id() const { return mID; }

};

template <typename Patch_id>
class Polyhedron_demo_items : public CGAL::Polyhedron_items_3 {
public:
  // wrap vertex
  template<class Refs, class Traits> struct Vertex_wrapper
  {
    typedef typename Traits::Point_3 Point;
    typedef Polyhedron_demo_vertex<Refs,
      CGAL::Tag_true,
      Point,
      Patch_id> Vertex;
  };

  // wrap face
  template<class Refs, class Traits> struct Face_wrapper
  {
    typedef Polyhedron_demo_face<Refs,
      CGAL::Tag_true,
      typename Traits::Plane_3,
      Patch_id> Face;
  };

  // wrap halfedge
  template<class Refs, class Traits> struct Halfedge_wrapper
  {
    typedef Polyhedron_demo_halfedge<Refs,
      CGAL::Tag_true,
      CGAL::Tag_true,
      CGAL::Tag_true> Halfedge;
  };
};

#include "Polyhedron_type_fwd.h"

// surface mesh
typedef CGAL::Polyhedron_3<Kernel, Polyhedron_demo_items<int> > Polyhedron;

#endif // POLYHEDRON_TYPE_H
