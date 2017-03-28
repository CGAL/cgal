#ifndef POLYHEDRON_TYPE_H
#define POLYHEDRON_TYPE_H

// CGAL
// kernel
#include "Kernel_type.h"

// surface mesh
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Has_timestamp.h>
#include <CGAL/tags.h>

#include <CGAL/Polygon_mesh_processing/properties.h>

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
  std::size_t time_stamp_;

public:
  int nb_of_feature_edges;

  bool is_corner() const {
    return nb_of_feature_edges > 2;
  }

  bool is_feature_vertex() const {
    return nb_of_feature_edges != 0;
  }

  void add_incident_patch(const Patch_id& i) {
    indices.insert(i);
  }

  /// For the determinism of Compact_container iterators
  ///@{
  typedef CGAL::Tag_true Has_timestamp;

  std::size_t time_stamp() const {
    return time_stamp_;
  }
  void set_time_stamp(const std::size_t& ts) {
    time_stamp_ = ts;
  }
  ///}@

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
  std::size_t time_stamp_;

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

  /// For the determinism of Compact_container iterators
  ///@{
  typedef CGAL::Tag_true Has_timestamp;

  std::size_t time_stamp() const {
    return time_stamp_;
  }
  void set_time_stamp(const std::size_t& ts) {
    time_stamp_ = ts;
  }
  ///@}
};



template <typename Integral>
inline std::pair<Integral, Integral>
patch_id_default_value(std::pair<Integral, Integral>)
{
  return std::pair<Integral, Integral>(1, 0);
}

template <typename Integral>
inline Integral patch_id_default_value(Integral)
{
  return Integral(1);
}

template <class Refs, class T_, class Pln_, class Patch_id_>
class Polyhedron_demo_face : 
  public CGAL::HalfedgeDS_face_base<Refs,T_,Pln_>
{
private:
  Patch_id_ patch_id_;
  std::size_t mID;
  std::size_t time_stamp_;

public:
  typedef Patch_id_ Patch_id;
  
  Polyhedron_demo_face()
    : patch_id_(patch_id_default_value(Patch_id())), mID(-1) {}
  
  const Patch_id& patch_id() const {
    return patch_id_;
  }
  
  void set_patch_id(const Patch_id& i) {
    patch_id_ = i;
  }
  
  std::size_t& id()       { return mID; }
  std::size_t  id() const { return mID; }

  /// For the determinism of Compact_container iterators
  ///@{
  typedef CGAL::Tag_true Has_timestamp;

  std::size_t time_stamp() const {
    return time_stamp_;
  }
  void set_time_stamp(const std::size_t& ts) {
    time_stamp_ = ts;
  }
  ///@}
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
typedef Polyhedron_demo_items<Patch_id>              Polyhedron_items;
typedef CGAL::Polyhedron_3<Kernel, Polyhedron_items> Polyhedron;



struct Is_feature_pmap {

  typedef Polyhedron::Halfedge_handle key_type;
  typedef bool value_type;
  friend bool get(const Is_feature_pmap&, Polyhedron::Halfedge_handle h)
  {
    return h->is_feature_edge();
  }
  
  friend void put(const Is_feature_pmap&, Polyhedron::Halfedge_handle h, bool b)
  {
     h->set_feature_edge(b);
  }

};

inline Is_feature_pmap get(CGAL::halfedge_is_feature_t, const Polyhedron&)
{
  return Is_feature_pmap();
}

struct vertex_num_feature_edges_pmap {

  typedef Polyhedron::Vertex_handle value_type;
  typedef int key_type;
  friend int get(const vertex_num_feature_edges_pmap&, Polyhedron::Vertex_handle h)
  {
    return h->nb_of_feature_edges;
  }
  
  friend void put(const vertex_num_feature_edges_pmap&, Polyhedron::Vertex_handle h, int n)
  {
     h->nb_of_feature_edges = n;
  }

};

inline vertex_num_feature_edges_pmap get(CGAL::vertex_num_feature_edges_t, const Polyhedron&)
{
  return vertex_num_feature_edges_pmap();
}


struct Patch_id_pmap {

  typedef Polyhedron::Face_handle key_type;
  typedef Polyhedron::Facet::Patch_id value_type;
  typedef value_type reference;
  typedef boost::writable_property_map_tag                          category;

  friend Polyhedron::Facet::Patch_id get(const Patch_id_pmap&, Polyhedron::Face_handle h)
  {
    return h->patch_id();
  }
  
  friend void put(const Patch_id_pmap&, Polyhedron::Face_handle h,
                  Polyhedron::Facet::Patch_id pid)
  {
     h->set_patch_id(pid);
  }

};



inline Patch_id_pmap get(CGAL::face_patch_id_t, const Polyhedron&)
{
  return Patch_id_pmap();
}

inline
boost::property_map<Polyhedron, boost::vertex_index_t>::type
get(CGAL::vertex_selection_t, const Polyhedron& p)
{
  return get(boost::vertex_index,p);
}

inline
boost::property_map<Polyhedron, boost::face_index_t>::type
get(CGAL::face_selection_t, const Polyhedron& p)
{
  return get(boost::face_index,p);
}
namespace boost {
  
  template <>
  struct property_map<Polyhedron, CGAL::halfedge_is_feature_t>
  {
    typedef Is_feature_pmap type;
  };
  
  template <>
  struct property_map<Polyhedron, CGAL::face_patch_id_t>
  {
    typedef Patch_id_pmap type;
  };
  
  template <>
  struct property_map<Polyhedron, CGAL::face_selection_t>
  {
    typedef boost::property_map<Polyhedron,boost::face_index_t>::type type;
    typedef boost::property_map<Polyhedron,boost::face_index_t>::const_type const_type;
  };

  template <>
  struct property_map<Polyhedron, CGAL::vertex_selection_t>
  {
    typedef boost::property_map<Polyhedron,boost::vertex_index_t>::type type;
    typedef boost::property_map<Polyhedron,boost::vertex_index_t>::const_type const_type;
  };
  template <>
  struct property_map<Polyhedron, CGAL::vertex_num_feature_edges_t>
  {
    typedef vertex_num_feature_edges_pmap type;
  };
}
#endif // POLYHEDRON_TYPE_H
