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
#include <CGAL/Nef_3/SNC_iteration.h>
#include <CGAL/Lazy_kernel.h>
#include <CGAL/Cartesian.h>

#include <boost/container/deque.hpp>

#include <sstream>
#include <string>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 503
#include <CGAL/Nef_2/debug.h>

namespace CGAL {


template <class Traits>
class K3_tree
{

template <typename Kernel, typename Vertex, typename Coordinate>
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

private:
  Coordinate coord;
};

public:
  friend class Objects_along_ray;
  friend class Objects_around_segment;

public:

typedef typename Traits::SNC_decorator SNC_decorator;
typedef typename Traits::Infimaximal_box Infimaximal_box;
typedef typename Traits::Vertex_handle Vertex_handle;
typedef typename Traits::Vertex_list Vertex_list;
typedef typename Vertex_list::iterator Vertex_iterator;
typedef typename Vertex_list::const_iterator Vertex_const_iterator;
typedef typename Traits::Halfedge_handle Halfedge_handle;
typedef typename Traits::Halffacet_handle Halffacet_handle;
typedef typename Traits::Object_handle Object_handle;
typedef typename Traits::Object_list Object_list;
typedef typename Traits::Halfedge_list Halfedge_list;
typedef typename Halfedge_list::const_iterator Halfedge_const_iterator;
typedef typename Traits::Halffacet_list Halffacet_list;
typedef typename Halffacet_list::const_iterator Halffacet_const_iterator;

typedef typename Traits::Point_3 Point_3;
typedef typename Traits::Segment_3 Segment_3;
typedef typename Traits::Ray_3 Ray_3;
typedef typename Traits::Vector_3 Vector_3;
typedef typename Traits::Plane_3 Plane_3;
typedef typename Traits::Triangle_3 Triangle_3;
typedef typename Traits::Aff_transformation_3 Aff_transformation_3;

typedef typename Traits::Bounding_box_3 Bounding_box_3;
typedef typename Traits::Side_of_plane Side_of_plane;

typedef typename Traits::Kernel Kernel;
typedef typename Kernel::RT RT;
typedef typename Kernel::FT FT;

typedef Smaller_than<
  Kernel,
  Vertex_handle,
  int> Smaller;

  class Node  {
  friend class K3_tree<Traits>;
public:
    typedef Node* Node_handle;
  Node(const Vertex_list& V, const Halfedge_list& E, const Halffacet_list& F) :
    left_node(nullptr), right_node(nullptr), vertex_list(V), edge_list(E), facet_list(F)
  {
  }

  Node(Node_handle l, Node_handle r, const Plane_3& pl) :
    left_node(l), right_node(r), splitting_plane(pl)
  {
  }

  bool is_leaf() const {
    CGAL_assertion( (left_node != nullptr && right_node != nullptr) ||
                    (left_node == nullptr && right_node == nullptr));
    return (left_node == nullptr && right_node == nullptr);
  }

  Node_handle left() const { return left_node; }
  Node_handle right() const { return right_node; }
  const Plane_3& plane() const { return splitting_plane; }

  bool empty() { return vertex_list.empty() && edge_list.empty() && facet_list.empty(); }
  Vertex_const_iterator vertices_begin() { return vertex_list.begin(); }
  Vertex_const_iterator vertices_end() { return vertex_list.end(); }
  Halfedge_const_iterator edges_begin() { return edge_list.begin(); }
  Halfedge_const_iterator edges_end() { return edge_list.end(); }
  Halffacet_const_iterator facets_begin() { return facet_list.begin(); }
  Halffacet_const_iterator facets_end() { return facet_list.end(); }

  void transform(const Aff_transformation_3& t) {
    if(left_node != nullptr) {
        CGAL_assertion(right_node != nullptr);
        left_node->transform(t);
        right_node->transform(t);
        splitting_plane = splitting_plane.transform(t);
    }
  }

  void add_facet(Halffacet_handle f, int depth) {
    if(left_node == nullptr) {
      facet_list.push_back(f);
      return;
    }

    Side_of_plane sop(splitting_plane.point(), depth%3);
    Oriented_side side = sop(f);
    if( side == ON_NEGATIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
      left_node->add_facet(f, depth+1);
    if( side == ON_POSITIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
      right_node->add_facet(f, depth+1);
  }

  void add_edge(Halfedge_handle e, int depth) {
    if(left_node == nullptr) {
      edge_list.push_back(e);
      return;
    }

    Side_of_plane sop(splitting_plane.point(), depth%3);
    Oriented_side side = sop(e);
    if( side == ON_NEGATIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
      left_node->add_edge(e, depth+1);
    if( side == ON_POSITIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
      right_node->add_edge(e, depth+1);
  }

  void add_vertex(Vertex_handle v, int depth) {
    if(left_node == nullptr) {
      vertex_list.push_back(v);
      return;
    }

    Side_of_plane sop(splitting_plane.point(), depth%3);
    Oriented_side side = sop(v);
    if( side == ON_NEGATIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
      left_node->add_vertex(v, depth+1);
    if( side == ON_POSITIVE_SIDE || side == ON_ORIENTED_BOUNDARY)
      right_node->add_vertex(v, depth+1);
  }


private:

  Node_handle left_node;
  Node_handle right_node;
  Plane_3 splitting_plane;
  Vertex_list vertex_list;
  Halfedge_list edge_list;
  Halffacet_list facet_list;
};

  typedef boost::container::deque<Node> Node_range;
  typedef Node* Node_handle;
  typedef std::vector<Node_handle> Node_list;


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
    Objects_around_segment() : root_node(nullptr), initialized(false) {}
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
      Iterator() : node()
      { CGAL_assertion_code( first_segment = false); }
      Iterator( const Node_handle root, const Segment_3& s) {
        CGAL_assertion_code( first_segment = true);
        S.push_front( Candidate( root, s));
        ++(*this); // place the iterator in the first intersected cell
      }
      Iterator( const Self& i) : S(i.S), node(i.node) {}
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


private:
#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
  bool reference_counted;
#endif
  Traits traits;


  Node_handle root;
  Node_range nodes;

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

    Vertex_list vertices;
    vertices.reserve(W->number_of_vertices());
    Vertex_iterator v;
    CGAL_forall_vertices(v, *W)
      vertices.push_back(v);

    Halfedge_list edges;
    edges.reserve(W->number_of_halfedges());
    Halfedge_iterator e;
    CGAL_forall_edges(e, *W)
      edges.push_back(e);

    Halffacet_list facets;
    facets.reserve(W->number_of_halffacets());
    Halffacet_iterator f;
    CGAL_forall_facets(f, *W)
      facets.push_back(f);

    CGAL_NEF_TRACEN("K3_tree(): n_vertices = " << vertices.size());
    std::frexp( (double) vertices.size(), &max_depth);

    // TODO: in the presence of a infimaximal bounding box, the bounding box does not have to be computed
    bounding_box = Bounding_box_3();
    for(typename Vertex_list::iterator vi=vertices.begin(); vi!=vertices.end(); ++vi)
        bounding_box.extend((*vi)->point());
    //CGAL_NEF_TRACEN("bounding box:"<<objects_bbox);

#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
    Point_3 p1(1,2,7), p2(p1);
    reference_counted = (&(p1.hx()) == &(p2.hx()));
    CGAL_NEF_TRACEN("reference counted " << reference_counted);
#endif
    non_efective_splits=0;
    root = build_kdtree(vertices, edges, facets, 0);
  }
  Node_handle locate_node_containing( const Point_3& p) const {
    return locate_node_containing( p, root);
  }
  Node_list nodes_along_ray( const Ray_3& r) const {
    Segment_3 s = ray_to_segment(r);
    return nodes_around_segment(s);
  }
  Node_list nodes_around_segment( const Segment_3& s) const {
    Node_list result;
    Objects_around_segment objects( *this, s);
    for(typename Objects_around_segment::Iterator oas = objects.begin(); oas != objects.end(); ++oas) {
      result.push_back(oas.get_node());
    }
    return result;
  }
  bool is_point_in_node( const Point_3& p, const Node_handle target) const {
    return is_point_in_node( p, target, root);
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
      for(Vertex_const_iterator vi = n->vertex_list.begin(); vi!=n->vertex_list.end(); ++vi) {
          b.extend((*vi)->point());
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

private:

int non_efective_splits;

Node_handle build_kdtree(Vertex_list& V, Halfedge_list& E, Halffacet_list& F,
                         int depth) {
  CGAL_precondition( depth >= 0);

  if( !can_set_be_divided(depth, V.size())) {
    CGAL_NEF_TRACEN("build_kdtree: set cannot be divided");
    nodes.push_back(Node(V, E, F));
    return &(nodes.back());
  }

  int coord = depth%3;
  Point_3 point_on_plane = find_median_point(V, coord);
  CGAL_NEF_TRACEN("build_kdtree: plane: "<<partition_plane<< " " << point_on_plane);

#ifdef CGAL_NEF_EXPLOIT_REFERENCE_COUNTING
  Side_of_plane sop(point_on_plane, coord, reference_counted);
#else
  Side_of_plane sop(point_on_plane, coord);
#endif
  sop.reserve(V.size());

  Vertex_list V1,V2;
  classify_objects(V, sop, V1, V2);

  Halfedge_list E1,E2;
  bool splitted = classify_objects(E, sop, E1, E2);

  Halffacet_list F1,F2;
  splitted = splitted && classify_objects(F, sop, F1, F2);

  if(!splitted) {
    CGAL_NEF_TRACEN("build_kdtree: splitting plane not found");
    nodes.push_back(Node(V, E, F));
    return &(nodes.back());
  }
  auto O_size = V.size() + E.size() + F.size();
  auto O1_size = V1.size() + E1.size() + F1.size();
  auto O2_size = V2.size() + E2.size() + F2.size();
  CGAL_NEF_TRACEN("Sizes " << O1_size << ", " << O2_size << ", " << O_size);
  CGAL_assertion( O1_size <= O_size && O2_size <= O_size);
  CGAL_assertion( O1_size + O2_size >= O_size);
  if((O1_size == O_size) || (O2_size == O_size))
    non_efective_splits++;
  else
    non_efective_splits = 0;

  if(non_efective_splits > 2) {
    CGAL_NEF_TRACEN("build_kdtree: non effective splits reached maximum");
    nodes.push_back(Node(V, E, F));
    return &(nodes.back());
  }

  Node_handle left_node = build_kdtree(V1, E1, F1, depth + 1);
  Node_handle right_node = build_kdtree(V2, E2, F2, depth + 1);
  nodes.push_back(Node(left_node, right_node, construct_splitting_plane(point_on_plane, coord, typename Traits::Kernel::Kernel_tag())));
  return &(nodes.back());
}

bool can_set_be_divided(int depth, typename Vertex_list::size_type size) {
  if(depth >= max_depth)
    return false;
  if(size <= 2)
    return false;
  return true;
}

template <typename List>
static bool classify_objects(const List& l,Side_of_plane& sop,
                      List& l1, List& l2) {
  typename List::size_type on_oriented_boundary = 0;
  for(typename List::const_iterator i=l.begin(); i!=l.end(); ++i) {
    Oriented_side side = sop(*i);
    switch(side) {
      case ON_NEGATIVE_SIDE:
        l1.push_back(*i);
        break;
      case ON_POSITIVE_SIDE:
        l2.push_back(*i);
        break;
      case ON_ORIENTED_BOUNDARY:
        ++on_oriented_boundary;
        l1.push_back(*i);
        l2.push_back(*i);
    }
  }
  return (on_oriented_boundary != l.size());
}

static Point_3 find_median_point(Vertex_list& V, int coord) {
  CGAL_assertion(V.size() > 1);

  Smaller smaller(coord);
  Vertex_iterator begin = V.begin();
  Vertex_iterator median = begin + V.size()/2;
  std::nth_element(begin, median, V.end(), smaller);
  Vertex_iterator prev = std::prev(median);
  std::nth_element(begin, prev, median, smaller);

  return CGAL::midpoint((*median)->point(), (*prev)->point());
}

static Plane_3 construct_splitting_plane(const Point_3& pt, int coord, const Homogeneous_tag&)
{
  switch(coord) {
  case 0: return Plane_3(pt, Vector_3(1, 0, 0));
  case 1: return Plane_3(pt, Vector_3(0, 1, 0));
  case 2: return Plane_3(pt, Vector_3(0, 0, 1));
  }

  CGAL_error_msg( "never reached");
  return Plane_3();
}

static Plane_3 construct_splitting_plane(const Point_3& pt, int coord, const Cartesian_tag&)
{
  switch(coord) {
  case 0: return Plane_3(1, 0, 0, -pt.x());
  case 1: return Plane_3(0, 1, 0, -pt.y());
  case 2: return Plane_3(0, 0, 1, -pt.z());
  }

  CGAL_error_msg( "never reached");
  return Plane_3();
}

static Node_handle get_child_by_side( const Node_handle node, Oriented_side side) {
  CGAL_assertion( node != nullptr);
  CGAL_assertion( side != ON_ORIENTED_BOUNDARY);
  if( side == ON_NEGATIVE_SIDE) {
    return node->left();
  }
  CGAL_assertion( side == ON_POSITIVE_SIDE);
  return node->right();
}

Node_handle locate_node_containing( const Point_3& p, const Node_handle node) const {
  CGAL_precondition( node != nullptr);
  if( node->is_leaf())
    return node;

  Oriented_side side = node->plane().oriented_side(p);
  if(side == ON_ORIENTED_BOUNDARY)
    side = ON_NEGATIVE_SIDE;
  return locate_node_containing(p, get_child_by_side(node, side));
}


bool is_point_in_node( const Point_3& p, const Node_handle target, const Node_handle current) const {
  CGAL_precondition( target != nullptr && current != nullptr);
  if( current->is_leaf())
    return (current == target);
  Oriented_side side = current->plane().oriented_side(p);
  if( side == ON_NEGATIVE_SIDE)
    return is_point_in_node( p, target, current->left());
  else if( side == ON_POSITIVE_SIDE)
    return is_point_in_node( p, target, current->right());
  CGAL_assertion( side == ON_ORIENTED_BOUNDARY);
  return (is_point_in_node( p, target, current->left()) ||
          is_point_in_node( p, target, current->right()));
}

Segment_3 ray_to_segment(const Ray_3& r) const
{
  CGAL_NEF_TRACEN("Objects_along_ray: input ray: "<<r);
  Vector_3 vec(r.to_vector());
  /* First of all, we need to find out whether we are working over an extended
   * kernel or on a standard kernel. As precondition we have that ray is oriented
   * in the minus x axis direction.  When having an extended kernel, the ray can
   * be substituted by a segment with the endpoint on the 'intersection' between
   * the ray and the bounding infimaximal box.  In the presence of a standard
   * kernel, the intersection is computed with the bounding box with the vertices
   * of the Nef polyhedron.*/
  Point_3 p(r.source()), q;
  Bounding_box_3 b = bounding_box;
  typename Kernel::Non_zero_coordinate_index_3 non_zero_coordinate_index_3;
  int c = non_zero_coordinate_index_3(vec);

  Point_3 pt_on_minus_x_plane = vec[c] < 0 ?
    Point_3(FT(b.min_coord(0)), FT(b.min_coord(1)),FT(b.min_coord(2))) :
    Point_3(FT(b.max_coord(0)), FT(b.max_coord(1)),FT(b.max_coord(2)));
  /* We compute the intersection between a plane with normal vector in the minus x
   * direction and located at the minimum point of the bounding box, and the input
   * ray.  When the ray does not intersect the bounding volume, there won't be any
   * object hit, so it is safe to construct a segment that simply lay in the
   * unbounded side of the bounding box.  This approach is taken instead of somehow
   * (efficiently) report that there was no hit object, in order to maintain a clear
   * interface with the Iterator class.*/
  Plane_3 pl_on_minus_x = K3_tree::construct_splitting_plane(pt_on_minus_x_plane, c, typename Traits::Kernel::Kernel_tag());
  Object o = traits.intersect_object()( pl_on_minus_x, r);
  if( !CGAL::assign( q, o) || pl_on_minus_x.has_on(p))
    q = r.source() + vec;
  else
    q = normalized(q);
  return Segment_3( p, q);
}

};

} //namespace CGAL

#endif // CGAL_NEF_K3_TREE_H
