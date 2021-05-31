// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Miguel Granados <granados@mpi-sb.mpg.de>
//                 Peter Hachenberger <hachenb@mpi-sb.mpg.de>

#ifndef CGAL_NEF_K3_TREE_H
#define CGAL_NEF_K3_TREE_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_3/quotient_coordinates_to_homogeneous_point.h>
#include <CGAL/Lazy_kernel.h>
#include <CGAL/Cartesian.h>
#include <boost/container/deque.hpp>


#include <deque>
#include <sstream>
#include <string>
#include <map>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 503
#include <CGAL/Nef_2/debug.h>

namespace CGAL {


template <class Traits>
class K3_tree
{

template <typename Kernel, typename Object,
          typename Vertex, typename Coordinate>
class Smaller_than
{
public:
  Smaller_than(Coordinate c) : coord(c) {
    CGAL_assertion( c >= 0 && c <=2);
  }
  bool operator()( const Vertex& v1, const Vertex& v2) {
    switch(coord) {
    case 0: return CGAL::compare_x(v1->point(), v2->point()) == SMALLER;
    case 1: return CGAL::compare_y(v1->point(), v2->point()) == SMALLER;
    case 2: return CGAL::compare_z(v1->point(), v2->point()) == SMALLER;
    default: CGAL_error();
    }
    return false;
  }

  bool operator()( const Object& o1, const Object& o2) {
    Vertex v1,v2;
    CGAL::assign(v1,o1);
    CGAL::assign(v2,o2);
    switch(coord) {
    case 0: return CGAL::compare_x(v1->point(), v2->point()) == SMALLER;
    case 1: return CGAL::compare_y(v1->point(), v2->point()) == SMALLER;
    case 2: return CGAL::compare_z(v1->point(), v2->point()) == SMALLER;
    default: CGAL_error();
    }
    return false;
  }
private:
  Coordinate coord;
};


template <typename Object, typename Vertex,
          typename Coordinate, typename EK>
  class Smaller_than<CGAL::Lazy_kernel<EK>, Object, Vertex, Coordinate>
{
public:
  Smaller_than(Coordinate c) : coord(c) {
    CGAL_assertion( c >= 0 && c <=2);
  }
  bool operator()( const Vertex& v1, const Vertex& v2) {
    switch(coord) {
    case 0: return CGAL::to_interval(v1->point().x()).second <
                           CGAL::to_interval(v2->point().x()).first;
    case 1: return CGAL::to_interval(v1->point().y()).second <
                           CGAL::to_interval(v2->point().y()).first;
    case 2: return CGAL::to_interval(v1->point().z()).second <
                           CGAL::to_interval(v2->point().z()).first;
    default: CGAL_error();
    }
    return false;
  }

  bool operator()( const Object& o1, const Object& o2) {
    Vertex v1,v2;
    CGAL::assign(v1,o1);
    CGAL::assign(v2,o2);
    switch(coord) {
    case 0: return CGAL::to_interval(v1->point().x()).second <
                           CGAL::to_interval(v2->point().x()).first;
    case 1: return CGAL::to_interval(v1->point().y()).second <
                           CGAL::to_interval(v2->point().y()).first;
    case 2: return CGAL::to_interval(v1->point().z()).second <
                           CGAL::to_interval(v2->point().z()).first;
    default: CGAL_error();
    }
    return false;
  }
private:
  Coordinate coord;
};


public:
  friend class Objects_along_ray;
  friend class Objects_around_segment;
  friend class Objects_around_box;

public:

typedef typename Traits::SNC_decorator SNC_decorator;
typedef typename Traits::Infimaximal_box Infimaximal_box;
typedef typename Traits::Vertex_handle Vertex_handle;
typedef typename Traits::Halfedge_handle Halfedge_handle;
typedef typename Traits::Halffacet_handle Halffacet_handle;
typedef typename Traits::Object_handle Object_handle;
typedef std::vector<Object_handle> Object_list;
typedef typename Object_list::const_iterator Object_const_iterator;
typedef typename Object_list::iterator Object_iterator;
typedef typename Object_list::size_type size_type;

typedef typename Traits::Point_3 Point_3;
typedef typename Traits::Segment_3 Segment_3;
typedef typename Traits::Ray_3 Ray_3;
typedef typename Traits::Vector_3 Vector_3;
typedef typename Traits::Plane_3 Plane_3;
typedef typename Traits::Triangle_3 Triangle_3;
typedef typename Traits::Aff_transformation_3 Aff_transformation_3;

typedef typename Traits::Bounding_box_3 Bounding_box_3;
typedef typename Traits::Side_of_plane Side_of_plane;
typedef typename Traits::Objects_bbox Objects_bbox;

typedef typename Traits::Kernel Kernel;
typedef typename Kernel::RT RT;
typedef typename Kernel::FT FT;

typedef Smaller_than<
  Kernel,
  Object_handle,
  Vertex_handle,
  int> Smaller_;

  class Node  {
  friend class K3_tree<Traits>;
public:
    typedef Node* Node_handle;

  Node( Node_handle p, Node_handle l, Node_handle r, Plane_3 pl, const Object_list& O) :
    parent_node(p), left_node(l), right_node(r), splitting_plane(pl),
        object_list(O) {
    if(l == nullptr)
      point_on_plane = Point_3();
    else
      point_on_plane = pl.point();
  }
  bool is_leaf() const {
    CGAL_assertion( (left_node != nullptr && right_node != nullptr) ||
                    (left_node == nullptr && right_node == nullptr));
    return (left_node == nullptr && right_node == nullptr);
  }
  Node_handle parent() const { return parent_node; }
  Node_handle left() const { return left_node; }
  Node_handle right() const { return right_node; }
  const Plane_3& plane() const { return splitting_plane; }
  const Object_list& objects() const { return object_list; }

  void transform(const Aff_transformation_3& t) {
    if(left_node != nullptr) {
        CGAL_assertion(right_node != nullptr);
        left_node->transform(t);
         right_node->transform(t);
          splitting_plane = splitting_plane.transform(t);
    }
  }

  std::size_t bytes() {
    // bytes used for the Kd-tree
    std::size_t s = sizeof(Node);
    if(left_node != nullptr)
      s += left_node->bytes();
    if(right_node != nullptr)
      s += right_node->bytes();
    typename Object_list::iterator o;
    for(o = object_list.begin(); o != object_list.end(); ++o)
      s += sizeof(*o);
    return s;
  }

  std::size_t leafs(int mask = 255, int lower_limit = 0) {
    std::size_t s = 0;
    Halffacet_handle f;
    Halfedge_handle e;
    Vertex_handle v;
    typename Object_list::iterator o;
    if(mask == 0)
      s = 1;
    else {
      for(o = object_list.begin(); o != object_list.end(); ++o) {
        if((mask & 1) && assign(v,*o))
          ++s;
        else if((mask&2) && assign(e,*o))
          ++s;
        else if(((mask&4) || (mask&8)) && assign(f,*o)) {
          if(mask&4)
            ++s;
          else {
            int length = 0;
            typename Traits::SHalfedge_around_facet_circulator safc(f->facet_cycles_begin()),
              send(safc);
            while(++length < lower_limit && ++safc != send) ;
            if(length >= lower_limit)
              ++s;
          }
        }
      }
    }

    if(left_node != nullptr)
      s += left_node->leafs(mask, lower_limit);
    if(right_node != nullptr)
      s += right_node->leafs(mask, lower_limit);
    return s;
  }

  template<typename Depth>
  void add_facet(Halffacet_handle f, Depth depth) {
    if(left_node == nullptr) {
      object_list.push_back(make_object(f));
      return;
    }

    Side_of_plane sop;
    Oriented_side side = sop(splitting_plane.point(), f, depth);
    if( side == ON_NEGATIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
      left_node->add_facet(f, depth+1);
    if( side == ON_POSITIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
      right_node->add_facet(f, depth+1);
  }

  template<typename Depth>
  void add_edge(Halfedge_handle e, Depth depth) {
    if(left_node == nullptr) {
      object_list.push_back(make_object(e));
      return;
    }

    Side_of_plane sop;
    Oriented_side side = sop(splitting_plane.point(), e, depth);
    if( side == ON_NEGATIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
      left_node->add_edge(e, depth+1);
    if( side == ON_POSITIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
      right_node->add_edge(e, depth+1);
  }

  template<typename Depth>
  void add_vertex(Vertex_handle v, Depth depth) {
    if(left_node == nullptr) {
      object_list.push_back(make_object(v));
      return;
    }

    Side_of_plane sop;
    Oriented_side side = sop(splitting_plane.point(), v, depth);
    if( side == ON_NEGATIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
      left_node->add_vertex(v, depth+1);
    if( side == ON_POSITIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
      right_node->add_vertex(v, depth+1);
  }


friend std::ostream& operator<<
  (std::ostream& os, const Node_handle node) {
  CGAL_assertion( node != nullptr);
  if( node->is_leaf())
    os <<  node->objects().size();
  else {
    os << " ( ";
    if( !node->left()) os << '-';
    else os << node->left();
    os << " , ";
    if( !node->right()) os << '-';
    else os << node->right();
    os << " ) ";
  }
  return os;
}

    /*
~Node() noexcept(!CGAL_ASSERTIONS_ENABLED)
{
  CGAL_NEF_TRACEN("~Node: deleting node...");
  CGAL_destructor_assertion_catch(
    if( !is_leaf()) {
      delete left_node;
      delete right_node;
    }
  );
}
    */
private:



  Node_handle parent_node;
  Node_handle left_node;
  Node_handle right_node;
  Plane_3 splitting_plane;
  Point_3 point_on_plane;
  Object_list object_list;
};

  typedef boost::container::deque<Node> Node_range;
  typedef Node* Node_handle;


public:
  class Objects_around_segment
  {
   public:
    class Iterator;
  protected:
    Traits traits;
    Node_handle root_node;
    Segment_3 segment;
    bool initialized;
  public:
    Objects_around_segment() : initialized(false) {}
    Objects_around_segment( const K3_tree& k, const Segment_3& s) :
      root_node(k.root), segment(s), initialized(true) {
      CGAL_NEF_TRACEN("Objects_around_segment: input segment: "<<segment);
    }
    void initialize( const K3_tree& k, const Segment_3& s) {
      root_node = k.root;
      segment = s;
      initialized = true;
      CGAL_NEF_TRACEN("Objects_around_segment: input segment: "<<s<<" (initialize)");
    }
  public:
    Iterator begin() const {
      CGAL_assertion( initialized == true);
      return Iterator( root_node, segment);
    }
    Iterator end() const {
      return Iterator();
    }
    class Iterator
    {
      friend class K3_tree;
      typedef Iterator Self;
      typedef std::pair< const Node_handle, Segment_3> Candidate;
    protected:
      std::list<Candidate> S;
      Node_handle node;
      Traits traits;
      CGAL_assertion_code( Segment_3 prev_segment;)
      CGAL_assertion_code( bool first_segment;)
    public:
      Iterator() : node() {}
      Iterator( const Node_handle root, const Segment_3& s) {
        CGAL_assertion_code( first_segment = true);
        S.push_front( Candidate( root, s));
        ++(*this); // place the interator in the first intersected cell
      }
      Iterator( const Self& i) : S(i.S), node(i.node) {}
      const Object_list& operator*() const {
        CGAL_assertion( node != nullptr);
        return node->objects();
      }
      Self& operator++() {

if( S.empty())
  node = nullptr; // end of the iterator
else {
  while( !S.empty()) {
    Node_handle n = S.front().first;
    Segment_3 s = S.front().second;
    S.pop_front();
    if( n->is_leaf()) {

      CGAL_assertion_code(
        if( first_segment) {
          first_segment = false;
          CGAL_NEF_TRACEN("operator++: prev_segment=(none), segment="<<s);
        } else {
          CGAL_assertion( prev_segment.target() == s.source());
          CGAL_assertion( prev_segment.direction() == s.direction());
          CGAL_NEF_TRACEN("operator++: prev_segment="<<prev_segment<<", segment="<<s);
        }
        prev_segment = s);

      node = n;
      break;
    }
    else {
      CGAL_NEF_TRACEN("find next intersected cell: segment: "<<s);
      CGAL_NEF_TRACEN("find next intersected cell: node plane: "<<n->plane() <<
             ", point: "<<n->plane().point());
      Oriented_side src_side = n->plane().oriented_side(s.source());
      Oriented_side tgt_side = n->plane().oriented_side(s.target());
      if( src_side == ON_ORIENTED_BOUNDARY && tgt_side == ON_ORIENTED_BOUNDARY)
        src_side = tgt_side = ON_NEGATIVE_SIDE;
      else if( src_side == ON_ORIENTED_BOUNDARY)
        src_side = tgt_side;
      else if( tgt_side == ON_ORIENTED_BOUNDARY)
        tgt_side = src_side;
      if( src_side == tgt_side)
        S.push_front( Candidate( get_child_by_side( n, src_side), s));
      else {
        Segment_3 s1, s2;
        divide_segment_by_plane( s, n->plane(), s1, s2);
        S.push_front( Candidate( get_child_by_side( n, tgt_side), s2)); // cell on target pushed first
        S.push_front( Candidate( get_child_by_side( n, src_side), s1));
      }
    }
  }
}

        return *this;
      }
      bool operator==(const Self& i) const {
        return (node == i.node);
      }
      bool operator!=(const Self& i) const {
        return !(*this == i);
      }
    private:
      Node_handle get_node() const {
        CGAL_assertion( node != nullptr);
        return node;
      }

inline
Node_handle get_child_by_side( const Node_handle node, Oriented_side side) {
  CGAL_assertion( node != nullptr);
  CGAL_assertion( side != ON_ORIENTED_BOUNDARY);
  if( side == ON_NEGATIVE_SIDE) {
    return node->left();
  }
  CGAL_assertion( side == ON_POSITIVE_SIDE);
  return node->right();
}

void divide_segment_by_plane( Segment_3 s, Plane_3 pl,
                              Segment_3& s1, Segment_3& s2) {
  Object o = traits.intersect_object()( pl, s);
  Point_3 ip;
  CGAL_assertion( CGAL::assign( ip, o));
  CGAL::assign( ip, o);
  ip = normalized(ip);
  s1 = Segment_3( s.source(), ip);
  s2 = Segment_3( ip, s.target());
  CGAL_assertion( s1.target() == s2.source());
  CGAL_assertion( s1.direction() == s.direction());
  CGAL_assertion( s2.direction() == s.direction());
}

    };
  };

  class Objects_along_ray : public Objects_around_segment
  {
    typedef Objects_around_segment Base;
  protected:
    Traits traits;
  public:
    Objects_along_ray( const K3_tree& k, const Ray_3& r) {
      CGAL_NEF_TRACEN("Objects_along_ray: input ray: "<<r);
      Vector_3 vec(r.to_vector());
      // First of all, we need to find out wheather we are working over an extended kernel or on a standard kernel. As precondition we have that ray is oriented in the minus x axis direction.  When having an extended kernel, the ray can be subtituted by a segment with the endpoint on the 'intersection' between the ray and the bounding infimaximal box.  In the presence of a standard kernel, the intersection is computed with the bounding box with the vertices of the Nef polyhedron.
      Point_3 p(r.source()), q;
      Bounding_box_3 b = k.bounding_box;
      int c = (CGAL::abs(vec[0]) > CGAL::abs(vec[1]) ? 0 : 1);
      c = (CGAL::abs(vec[2]) > CGAL::abs(vec[c]) ? 2 : c);

      Point_3 pt_on_minus_x_plane = vec[c] < 0 ?
        Point_3(FT(b.min_coord(0)), FT(b.min_coord(1)),FT(b.min_coord(2))) :
        Point_3(FT(b.max_coord(0)), FT(b.max_coord(1)),FT(b.max_coord(2)));
      // We compute the intersection between a plane with normal vector in
      // the minus x direction and located at the minimum point of the bounding box, and the input ray.  When the ray does not intersect the bounding volume, there won't be any object hit, so it is safe to construct a segment that simply lay in the unbounded side of the bounding box.  This approach is taken instead of somehow (efficiently) report that there was no hit object, in order to mantain a clear interface with the Iterator class.
      Plane_3 pl_on_minus_x;
      if(c==0)
        pl_on_minus_x = Plane_3(pt_on_minus_x_plane, Vector_3( 1, 0, 0));
      else if(c==1)
        pl_on_minus_x = Plane_3(pt_on_minus_x_plane, Vector_3( 0, 1, 0));
      else {
        CGAL_assertion_msg(c==2, "wrong value");
         pl_on_minus_x = Plane_3(pt_on_minus_x_plane, Vector_3( 0, 0, 1));
      }
      Object o = traits.intersect_object()( pl_on_minus_x, r);
      if( !CGAL::assign( q, o) || pl_on_minus_x.has_on(p))
        q = r.source() + vec;
      else
        q = normalized(q);
      Base::initialize( k, Segment_3( p, q));
    }
  };

class Objects_around_box {

 public:
  class Iterator;
 protected:
  Node_handle root_node;
  Bounding_box_3 box;
  bool initialized;

 public:
  Objects_around_box() : initialized(false) {}
  Objects_around_box(const K3_tree& k, const Bounding_box_3& b) :
    root_node(k.root), box(b), initialized(true) {}

  void initialize( const K3_tree& k, const Bounding_box_3& b) {
    root_node = k.root;
    box = b;
    initialized = true;
  }

 public:
  Iterator begin() const {
    CGAL_assertion( initialized == true);
    return Iterator( root_node, box);
  }

  Iterator end() const {
    return Iterator();
  }

  class Iterator {

    friend class K3_tree;
    typedef Iterator Self;
    typedef std::pair< const Node_handle, Bounding_box_3> Candidate;

  protected:
    std::list<Candidate> S;
    const Node_handle node;

  public:
    Iterator() : node(nullptr) {}

    Iterator( const Node_handle root, const Bounding_box_3& s) {
      S.push_front( Candidate( root, s));
      ++(*this); // place the interator in the first intersected cell
    }

    Iterator( const Self& i) : S(i.S), node(i.node) {}

    const Object_list& operator*() const {
      CGAL_assertion( node != nullptr);
      return node->objects();
    }

    Self& operator++() {

      if(S.empty())
        node = nullptr; // end of the iterator
      else {
        while( !S.empty()) {
          const Node_handle n = S.front().first;
          Bounding_box_3 b = S.front().second;
          S.pop_front();
          if( n->is_leaf()) {
            node = n;
            break;
          } else {
            Point_3 pmin(b.min_coord(0), b.min_coord(1), b.min_coord(2));
            Point_3 pmax(b.max_coord(0), b.max_coord(1), b.max_coord(2));
            Oriented_side src_side =
              n->plane().oriented_side(pmax);
            Oriented_side tgt_side =
              n->plane().oriented_side(pmin);
            if( src_side == tgt_side &&
                src_side != ON_ORIENTED_BOUNDARY)
              S.push_front( Candidate( get_child_by_side( n, src_side), b));
            else {
              S.push_front( Candidate( get_child_by_side( n, tgt_side), b));
              S.push_front( Candidate( get_child_by_side( n, src_side), b));
            }
          }
        }
      }
      return *this;
    }

    bool operator==(const Self& i) const {
      return (node == i.node);
    }

    bool operator!=(const Self& i) const {
      return !(*this == i);
    }

  private:
    Node_handle get_node() const {
      CGAL_assertion( node != nullptr);
      return node;
    }

    inline
    Node_handle get_child_by_side( const Node_handle node, Oriented_side side) {
      CGAL_assertion( node != nullptr);
      CGAL_assertion( side != ON_ORIENTED_BOUNDARY);
      if( side == ON_NEGATIVE_SIDE) {
        return node->left();
      }
      CGAL_assertion( side == ON_POSITIVE_SIDE);
      return node->right();
    }
  };
};

private:
#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
  bool reference_counted;
#endif
  Traits traits;


  Node_handle root;
  boost::container::deque<Node> nodes;

  int max_depth;
  Bounding_box_3 bounding_box;
public:
  template<typename SNC_structure>
  K3_tree(SNC_structure* W)
#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
    : reference_counted(false)
#endif
    {

    typedef typename SNC_structure::Vertex_iterator Vertex_iterator;
    typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
    typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator;

    CGAL_assertion( W != nullptr);
    Object_list objects;
    Vertex_iterator v;
    Halfedge_iterator e;
    Halffacet_iterator f;
    CGAL_forall_vertices( v, *W)
      objects.push_back(make_object(Vertex_handle(v)));
    typename Object_list::difference_type v_end = objects.size();
    CGAL_forall_edges( e, *W)
      objects.push_back(make_object(Halfedge_handle(e)));
    CGAL_forall_facets( f, *W) {
      objects.push_back(make_object(Halffacet_handle(f)));
    }
    Object_iterator oli=objects.begin()+v_end;
    root = build_kdtree( objects, oli, 0);
  }

  K3_tree(Object_list& objects, Object_iterator& v_end) {

typename Object_list::difference_type n_vertices = std::distance(objects.begin(),v_end);
 CGAL_NEF_TRACEN("K3_tree(): n_vertices = " << std::distance(objects.begin(),v_end));
 std::frexp( (double) n_vertices, &max_depth);

 // TODO: in the presence of a infimaximal bounding box, the bounding box does not have to be computed
 Objects_bbox objects_bbox = traits.objects_bbox_object();
 bounding_box = objects_bbox(objects);
 //CGAL_NEF_TRACEN("bounding box:"<<objects_bbox);

#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
    Point_3 p1(1,2,7), p2(p1);
    reference_counted = (&(p1.hx()) == &(p2.hx()));
    CGAL_NEF_TRACEN("reference counted " << reference_counted);
#endif
    root = build_kdtree( objects, v_end, 0);
  }
  const Object_list& objects_around_point( const Point_3& p) const {
    return locate( p, root);
  }
  Objects_along_ray objects_along_ray( const Ray_3& r) const {
    return Objects_along_ray( *this, r);
  }
  Object_list objects_around_segment( const Segment_3& s) const {
    Object_list O;

    Objects_around_segment objects( *this, s);
    Unique_hash_map< Vertex_handle, bool> v_mark(false);
    Unique_hash_map< Halfedge_handle, bool> e_mark(false);
    Unique_hash_map< Halffacet_handle, bool> f_mark(false);
    for( typename Objects_around_segment::Iterator oar = objects.begin();
         oar != objects.end(); ++oar) {
      for( typename Object_list::const_iterator o = (*oar).begin();
           o != (*oar).end(); ++o) { // TODO: implement operator->(...)
        Vertex_handle v;
        Halfedge_handle e;
        Halffacet_handle f;
        if( CGAL::assign( v, *o)) {
          if( !v_mark[v]) {
            O.push_back(*o);
            v_mark[v] = true;
          }
        }
        else if( CGAL::assign( e, *o)) {
          if( !e_mark [e]) {
            O.push_back(*o);
            e_mark[e] = true;
          }
        }
        else if( CGAL::assign( f, *o)) {
          if( !f_mark[f]) {
            O.push_back(*o);
            f_mark[f] = true;
          }
        }
        else
          CGAL_error_msg( "wrong handle");
      }
    }
    return O;
  }

  bool is_point_on_cell( const Point_3& p, const typename Objects_around_segment::Iterator& target) const {
    return is_point_on_cell( p, target.get_node(), root);
  }

  void add_facet(Halffacet_handle f) {
    root->add_facet(f,0);
  }

  void add_edge(Halfedge_handle e) {
    root->add_edge(e,0);
  }

  void add_vertex(Vertex_handle v) {
    root->add_vertex(v,0);
  }

  class BBox_updater {
    Bounding_box_3 b;

  public:
    BBox_updater() {}

    void pre_visit(const Node_handle) {}
    void post_visit(const Node_handle n) {
      typename Object_list::const_iterator o;
      for( o = n->objects().begin();
           o != n->objects().end(); ++o) {
        Vertex_handle v;
        if( CGAL::assign( v, *o))
          b.extend(v->point());
      }
    }

    Bounding_box_3 box() const{
      return b;
    }

  };

  template <typename Visitor>
  void visit_k3tree(const Node_handle current, Visitor& V) const {
    V.pre_visit(current);
    if(current->left() != nullptr) {
      visit_k3tree(current->left(), V);
      visit_k3tree(current->right(), V);
    }
    V.post_visit(current);
  }

  size_t bytes() { return root->bytes();}
  size_t leafs(int mask = 255, int lower_limit=0) { return root->leafs(mask, lower_limit);}

  void transform(const Aff_transformation_3& t) {
    if(root == nullptr){
      return;
    }
    root->transform(t);

    BBox_updater bbup;
    visit_k3tree(root, bbup);
    bounding_box = bbup.box();
  }


#ifdef CODE_DOES_NOT_WORK_WITH_BOTH_KERNELS_AT_THE_SAME_TIME
template <typename T>
friend std::ostream& operator<<
(std::ostream& os, const K3_tree<T>& k3_tree) {
  os << (const Node_handle)k3_tree.root; // no default conversion to const Node_handle?
  return os;
}
#endif
std::string dump_object_list( const Object_list& O, int level = 0) {
  std::stringstream os;
  typename Object_list::size_type v_count = 0, e_count = 0, f_count = 0;
  typename Object_list::const_iterator o;
  Vertex_handle v;
  Halfedge_handle e;
  Halffacet_handle f;
  for( o = O.begin(); o != O.end(); ++o) {
    if( CGAL::assign( v, *o)) {
      if( level) os << v->point() << std::endl;
      ++v_count;
    }
    else if( CGAL::assign( e, *o)) {
      if( level) os << e->source()->point() << "->"
                    << e->twin()->source()->point() << std::endl;
      ++e_count;
    }
    else if( CGAL::assign( f, *o)) {
      if( level) os << "facet" << std::endl;
      ++f_count;
    }
    else
      CGAL_error_msg( "wrong handle");
  }
  os << v_count << "v " << e_count << "e " << f_count << "f ";
  return os.str();
 }

bool update( Unique_hash_map<Vertex_handle, bool>& V,
             Unique_hash_map<Halfedge_handle, bool>& E,
             Unique_hash_map<Halffacet_handle, bool>& F) {
  return update( root, V, E, F);
}

bool update( Node_handle node,
             Unique_hash_map<Vertex_handle, bool>& V,
             Unique_hash_map<Halfedge_handle, bool>& E,
             Unique_hash_map<Halffacet_handle, bool>& F) {
  CGAL_assertion( node != nullptr);
  if( node->is_leaf()) {
    bool updated = false;
    Object_list* O = &node->object_list;
    typename Object_list::iterator onext, o = O->begin();
    while( o != O->end()) {
      onext = o;
      onext++;
      Vertex_handle v;
      Halfedge_handle e;
      Halffacet_handle f;
      if( CGAL::assign( v, *o)) {
        if( !V[v]) {
          O->erase(o);
          updated = true;
        }
      }
      else if( CGAL::assign( e, *o)) {
        if( !E[e]) {
          O->erase(o);
          updated = true;
        }
      }
      else if( CGAL::assign( f, *o)) {
        if( !F[f]) {
          O->erase(o);
          updated = true;
        }
      }
      else CGAL_error_msg( "wrong handle");
      o = onext;
    }
    return updated;
  }
  // TODO: protect the code below from optimizations!
  bool left_updated = update( node->left_node, V, E, F);
  CGAL_NEF_TRACEN("k3_tree::update(): left node updated? "<<left_updated);
  bool right_updated = update( node->right_node, V, E, F);
  CGAL_NEF_TRACEN("k3_tree::update(): right node updated? "<<right_updated);
  return (left_updated || right_updated);
}
  /*
~K3_tree() noexcept(!CGAL_ASSERTIONS_ENABLED)
{
  CGAL_NEF_TRACEN("~K3_tree: deleting root...");
  CGAL_destructor_assertion_catch(
    delete root;
  );
}
  */

private:

template <typename Depth>
Node_handle build_kdtree(Object_list& O, Object_iterator v_end,
                   Depth depth, Node_handle parent=nullptr, int non_efective_splits=0) {
  CGAL_precondition( depth >= 0);
  CGAL_NEF_TRACEN( "build_kdtree: "<<O.size()<<" objects, "<<"depth "<<depth);
  CGAL_NEF_TRACEN( "build_kdtree: "<<dump_object_list(O,1));
  if( !can_set_be_divided(O.begin(), v_end, depth)) {
    CGAL_NEF_TRACEN("build_kdtree: set cannot be divided");
    nodes.push_back(Node( parent, nullptr, nullptr, Plane_3(), O));
    return &(nodes.back());
  }
  Object_iterator median;
  Plane_3 partition_plane = construct_splitting_plane(O.begin(), v_end, median, depth);
  CGAL_NEF_TRACEN("build_kdtree: plane: "<<partition_plane<< " " << partition_plane.point());

  Object_list O1, O2;
  Vertex_handle vm,vx;
  CGAL::assign(vm,*median);
  Smaller_ smaller(depth%3);
  for(Object_iterator oi=O.begin();oi!=median;++oi) {
    O1.push_back(*oi);
    CGAL::assign(vx,*oi);
    if(!smaller(vx, vm))
      O2.push_back(*oi);
  }

  O1.push_back(*median);
  O2.push_back(*median);

  for(Object_iterator oi=median+1;oi!=v_end;++oi) {
    O2.push_back(*oi);
    CGAL::assign(vx,*oi);
    if(!smaller(vm, vx))
      O1.push_back(*oi);
  }

#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
  Side_of_plane sop(reference_counted);
#else
  Side_of_plane sop;
#endif
  typename Object_list::size_type v_end1 = O1.size();
  typename Object_list::size_type v_end2 = O2.size();
  bool splitted = classify_objects( v_end, O.end(), partition_plane, sop,
                                    std::back_inserter(O1),
                                    std::back_inserter(O2), depth);

  bool non_efective_split = false;
  if( !splitted) {
    CGAL_NEF_TRACEN("build_kdtree: splitting plane not found");
    //    if(depth > max_depth)
    nodes.push_back(Node( parent, nullptr, nullptr, Plane_3(), O));
    return &(nodes.back());
  } else {
    CGAL_NEF_TRACEN("Sizes " << O1.size() << ", " << O2.size() << ", " << O.size());
    CGAL_assertion( O1.size() <= O.size() && O2.size() <= O.size());
    CGAL_assertion( O1.size() + O2.size() >= O.size());
    non_efective_split = ((O1.size() == O.size()) || (O2.size() == O.size()));
  }
  if( non_efective_split)
    non_efective_splits++;
  else
    non_efective_splits = 0;
  if( non_efective_splits > 2) {
    CGAL_NEF_TRACEN("build_kdtree: non efective splits reached maximum");
    nodes.push_back(Node( parent, nullptr, nullptr, Plane_3(), O));
    return &(nodes.back());
  }
  nodes.push_back(Node( parent, nullptr, nullptr, partition_plane, Object_list()));
  Node_handle node = &(nodes.back());
  node->left_node = build_kdtree( O1, O1.begin()+v_end1, depth + 1, node, non_efective_splits);
  node->right_node = build_kdtree( O2, O2.begin()+v_end2, depth + 1, node, non_efective_splits);
  return node;
}

template <typename Depth>
bool can_set_be_divided(Object_iterator start, Object_iterator end, Depth depth) {
  if(depth >= max_depth)
    return false;
  if(std::distance(start,end)<2)
    return false;
  return true;
}

template <typename OutputIterator, typename Depth>
bool classify_objects(Object_iterator start, Object_iterator end,
                      Plane_3 partition_plane, Side_of_plane& sop,
                      OutputIterator o1, OutputIterator o2, Depth depth) {
  typename Object_list::difference_type on_oriented_boundary = 0;
  typename Object_list::const_iterator o;

  Point_3 point_on_plane(partition_plane.point());

  for( o = start; o != end; ++o) {
    Oriented_side side = sop( point_on_plane, *o, depth);
    if( side == ON_NEGATIVE_SIDE || side == ON_ORIENTED_BOUNDARY) {
      *o1 = *o;
      ++o1;
    }
    if( side == ON_POSITIVE_SIDE || side == ON_ORIENTED_BOUNDARY) {
      *o2 = *o;
      ++o2;
    }
    if( side == ON_ORIENTED_BOUNDARY)
      ++on_oriented_boundary;
  }
  return (on_oriented_boundary != std::distance(start,end));
}


template <typename Depth>
Plane_3 construct_splitting_plane(Object_iterator start, Object_iterator end,
                                  Object_iterator& median, Depth depth) {
  CGAL_precondition( depth >= 0);
  typename Object_list::difference_type n=std::distance(start,end);
  CGAL_assertion(n>1);

  std::nth_element(start, start+n/2, end,
                     Smaller_(depth%3));

  Vertex_handle v;
  median = start+n/2;
  CGAL::assign(v,*median);
  switch( depth % 3) {
  case 0: return Plane_3( v->point(), Vector_3( 1, 0, 0)); break;
  case 1: return Plane_3( v->point(), Vector_3( 0, 1, 0)); break;
  case 2: return Plane_3( v->point(), Vector_3( 0, 0, 1)); break;
  }

  CGAL_error_msg( "never reached");
  return Plane_3();
}

Node_handle locate_cell_containing( const Point_3& p, const Node_handle node) const {
  CGAL_precondition( node != nullptr);
  if( node->is_leaf())
    return node;
  else {
    Oriented_side side = node->plane().oriented_side(p);
    if( side == ON_NEGATIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
      return locate_cell_containing( p, node->left());
    else { // side == ON_POSITIVE_SIDE
      CGAL_assertion( side == ON_POSITIVE_SIDE);
      return locate_cell_containing( p, node->right());
    }
  }
}

const Object_list& locate( const Point_3& p, const Node_handle node) const {
  CGAL_precondition( node != nullptr);
  return locate_cell_containing( p, node)->objects();
}

bool is_point_on_cell( const Point_3& p, const Node_handle target, const Node_handle current) const {
  CGAL_precondition( target != nullptr && current != nullptr);
  if( current->is_leaf())
    return (current == target);
  Oriented_side side = current->plane().oriented_side(p);
  if( side == ON_NEGATIVE_SIDE)
    return is_point_on_cell( p, target, current->left());
  else if( side == ON_POSITIVE_SIDE)
    return is_point_on_cell( p, target, current->right());
  CGAL_assertion( side == ON_ORIENTED_BOUNDARY);
  return (is_point_on_cell( p, target, current->left()) ||
          is_point_on_cell( p, target, current->right()));
}

};

} //namespace CGAL

#endif // CGAL_NEF_K3_TREE_H
