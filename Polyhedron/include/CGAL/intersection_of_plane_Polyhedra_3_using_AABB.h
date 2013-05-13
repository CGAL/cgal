#ifndef CGAL_INTERSECTION_OF_PLANE_POLYHEDRA_3_USING_AABB_H
#define CGAL_INTERSECTION_OF_PLANE_POLYHEDRA_3_USING_AABB_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_segment_primitive.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>

#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Plane_3.h>

#include <CGAL/intersection_of_Polyhedra_3.h>

#include <vector>
#include <map>
#include <algorithm>

#include <boost/tuple/tuple.hpp>

namespace CGAL {

template<class Polyhedron, class PK>
class AABB_to_polyline_imp
{
private:
  typedef typename Polyhedron::Traits::Kernel K;
  typedef AABB_polyhedron_triangle_primitive<K, Polyhedron>       AABB_primitive;
  typedef AABB_traits<K, AABB_primitive>                          AABB_traits;
  typedef AABB_tree<AABB_traits>                                  AABB_tree;

  typedef typename AABB_tree::Object_and_primitive_id             Object_and_primitive_id;
  typedef typename AABB_tree::Primitive_id                        Primitive_id;

  typedef typename PK::Plane_3    Plane;
  typedef typename PK::Segment_3  Segment;
  typedef typename PK::Point_3    Point;

  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Facet_handle Facet_handle;

  // to store two intersecting halfedges of a facet
  struct Halfedge_pair {
    Halfedge_handle h1, h2;
    int halfedge_count;
    int halfedge_intersections;

    Halfedge_pair() : halfedge_count(0), halfedge_intersections(0) { }
    Halfedge_pair(Halfedge_handle h1, Halfedge_handle h2) : h1(h1), h2(h2), halfedge_count(2), halfedge_intersections(0)
    { }

    void put(Halfedge_handle h) {
      if(halfedge_count > 2) { CGAL_assertion(false && "degenerate case occured"); }
      (halfedge_count == 0 ? h1 : h2) = h;
      ++halfedge_count;
    }

    // halfedge `h` points to a neighbor facet
    Halfedge_handle get_other(Halfedge_handle h) {
      if(halfedge_count != 2) { CGAL_assertion(false && "degenerate case occured"); }
      return h->opposite() == h1 ? h2 : h1;
    }
    
    void halfedge_intersection_found() { 
      ++halfedge_intersections; 
      if(halfedge_intersections > 2) {
        std::cout << "halfedge_intersections: " << halfedge_intersections << std::endl;
        CGAL_assertion(false);        
      }
    }
  };

  bool halfedge_do_intersect(Halfedge_handle hf, const Plane& plane) {
    const Point& s = hf->vertex()->point();
    const Point& t = hf->opposite()->vertex()->point();
    return CGAL::do_intersect(plane, Segment(s,t));
  }

  Point halfedge_intersection(Halfedge_handle hf, const Plane& plane) {
    const Point& s = hf->vertex()->point();
    const Point& t = hf->opposite()->vertex()->point();
    Object intersection = CGAL::intersection(plane, Segment(s,t));

    const Point* i_point = 0; 
    if(!(i_point = object_cast<Point>(&intersection))) { 
      CGAL_assertion(false && "degenarate case occured");
    }

    return Point(*i_point);
  }

  typedef std::map<Facet_handle, Halfedge_pair> Facet_intersection_map;
  typedef typename Facet_intersection_map::iterator Facet_intersection_map_iterator;
  // members
  Facet_intersection_map facet_intersection_map;
  AABB_tree tree;

  void process_intersection(Halfedge_handle hf) {
    if(!hf->is_border()) {
      facet_intersection_map.find(hf->facet())
        ->second.halfedge_intersection_found();
    }
    if(!hf->opposite()->is_border()) {
      facet_intersection_map.find(hf->opposite()->facet())
        ->second.halfedge_intersection_found();
    }
  }

  template<class OutputIterator>
  void intersect_plane(const Plane& plane, OutputIterator out) 
  {
    facet_intersection_map.clear();

    std::vector<Primitive_id> intersected_triangles;
    tree.all_intersected_primitives(plane, std::back_inserter(intersected_triangles));

    // put two intersected edge for each facet
    for(typename std::vector<Primitive_id>::iterator it = intersected_triangles.begin();
      it != intersected_triangles.end(); ++it)
    {
      Halfedge_pair& halfedge_intersections = facet_intersection_map[*it];

      typename Polyhedron::Facet::Halfedge_around_facet_circulator edge_it = (*it)->facet_begin();
      do {
        if(halfedge_do_intersect(edge_it, plane))
        { halfedge_intersections.put(edge_it); }

      } while(++edge_it != (*it)->facet_begin());
    }

    // for each intersected facet, make it head and try to create a loop (works fine with border determination cases)
    for(typename std::vector<Primitive_id>::iterator it = intersected_triangles.begin();
      it != intersected_triangles.end(); ++it)
    {
      Facet_handle active_facet = *it;
      Facet_intersection_map_iterator active_facet_it = facet_intersection_map.find(active_facet);
      if(active_facet_it->second.halfedge_intersections == 2) { continue; } // we previously process this facet, continue
      CGAL_assertion(active_facet_it->second.halfedge_intersections == 0);

      std::vector<Point> polyline;

      // now create one polyline
      Facet_handle head_facet = active_facet;
      Halfedge_handle head_edge = active_facet_it->second.h1;
      Halfedge_handle active_edge = head_edge;

      bool loop_ok = true;
      while( true ) {
        // process active edge and delete active facet
        polyline.push_back(halfedge_intersection(active_edge, plane));
        process_intersection(active_edge); // marks neighbor facets (i.e. ++halfedge_intersections)

        // proceed to next facet, if edge is border then we know that we can not complete the loop
        if(active_edge->opposite()->is_border()) { 
          // this is the first border, so go to the other direction from head_facet
          if(loop_ok) {
            active_edge = facet_intersection_map.find(head_facet)
              ->second.get_other(head_edge->opposite());
            loop_ok = false;
            std::reverse(polyline.begin(), polyline.end());
            continue;
          }
          else { break; }
        }
        active_facet = active_edge->opposite()->facet();
        // check whether the loop is completed
        if(active_facet == head_facet) { 
          polyline.push_back(halfedge_intersection(head_edge, plane)); // put first edge again
          break;
        }

        Facet_intersection_map_iterator tmp_it = facet_intersection_map.find(active_facet);
        // CGAL_assertion(tmp_it != facet_intersection_map.end());
        active_edge = tmp_it->second.get_other(active_edge); 
      }

      *out++ = polyline;
    }
  }

public:
  template<class InputIterator, class OutputIterator>
  void operator()(Polyhedron& mesh, 
    InputIterator plane_begin,
    InputIterator plane_end,
    OutputIterator out) 
  {
    tree.rebuild(mesh.facets_begin(), mesh.facets_end());
    for(; plane_begin != plane_end; ++plane_begin) {
      intersect_plane(*plane_begin, out);
    }
  }
};


template<class PK, class Polyhedron, class OutputIterator>
void intersection_of_plane_Polyhedra_3_using_AABB(Polyhedron& mesh, 
  const Point_3<PK>& center, 
  const Vector_3<PK>& base1,
  const Vector_3<PK>& base2,
  OutputIterator out
  )
{
  Plane_3<PK> plane(center, center + base1, center + base2);
  intersection_of_plane_Polyhedra_3_using_AABB<PK>(mesh, plane, out);
}

template<class PK, class Polyhedron, class OutputIterator>
void intersection_of_plane_Polyhedra_3_using_AABB(Polyhedron& mesh, 
  const Plane_3<PK>& plane,
  OutputIterator out
  )
{
  const Plane_3<PK>* plane_begin = &plane;
  const Plane_3<PK>* plane_end = plane_begin + 1;
  AABB_to_polyline_imp<Polyhedron, PK>()(mesh, plane_begin, plane_end, out);
}

template<class PK, class Polyhedron, class InputIterator, class OutputIterator>
void intersection_of_plane_Polyhedra_3_using_AABB(Polyhedron& mesh, 
  InputIterator plane_begin,
  InputIterator plane_end,
  OutputIterator out) 
{
  AABB_to_polyline_imp<Polyhedron, PK>()(mesh, plane_begin, plane_end, out);
}

///////////////////////////////////////////////////////////////////////////
// not good, since segment of constructed 'plane Polyhedron' can intersect
template<class Polyhedron, class PK, class OutputIterator>
Polyhedron use_intersection_of_Polyhedra_3_with_plane(Polyhedron& mesh, 
  const Point_3<PK>& center, 
  const Vector_3<PK>& base1,
  const Vector_3<PK>& base2,
  OutputIterator out
  )
{
  // Construct a polyhedron consisting two triangles to represent plane
  Polyhedron plane;
  const Point_3<PK>& v0 = center + (base1 + base2);
  const Point_3<PK>& v1 = center + (base1 - base2);
  const Point_3<PK>& v2 = center + (base2 - base1);
  const Point_3<PK>& v3 = center + (-base2 - base1);

  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

  Halfedge_handle he_v0 = plane.make_triangle(v0, v1, v2);

  Halfedge_handle h = he_v0->opposite();
  Halfedge_handle g = h->next();
  Halfedge_handle he_v3 = plane.add_vertex_and_facet_to_border(h, g);
  he_v3->vertex()->point() = v3;
  
  intersection_Polyhedron_3_Polyhedron_3(plane, mesh, out); 
  return plane;
}

template<class Polyhedron, class PK, class OutputIterator>
Polyhedron use_intersection_of_Polyhedra_3_with_triangle(Polyhedron& mesh, 
  const Point_3<PK>& center, 
  const Vector_3<PK>& base1,
  const Vector_3<PK>& base2,
  OutputIterator out
  )
{
  // Construct a polyhedron consisting two triangles to represent plane
  Polyhedron plane;
  const Point_3<PK>& v0 = center + (base1 + base2);
  const Point_3<PK>& v1 = center + (base1 - base2);
  const Point_3<PK>& v2 = center + (base2 - base1);

  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;

  plane.make_triangle(v0, v1, v2);

  intersection_Polyhedron_3_Polyhedron_3(plane, mesh, out); 
  return plane;
}
///////////////////////////////////////////////////////////////////////////

}// end of namespace CGAL
#endif //CGAL_INTERSECTION_OF_PLANE_POLYHEDRA_3_USING_AABB_H

//enum Intersection_type { NONE, POINT, SEGMENT, TRIANGLE };

//private:
// use halfedge which has smaller address
//struct Edge_comparator {
//  bool operator()(Halfedge_handle h1, Halfedge_handle h2) const {
//    return (std::min)(&*h1, &*(h1->opposite())) 
//           < (std::min)(&*h2, &*(h2->opposite()));
//  }
//};

//public:
//template<class OutputIterator>
//void imp_1(Polyhedron& mesh, const Plane& plane, OutputIterator out) 
//{
//  AABB_tree tree(mesh.facets_begin(), mesh.facets_end());
//  std::vector<Primitive_id> intersected_triangles;
//  tree.all_intersected_primitives(plane, std::back_inserter(intersected_triangles));
//  
//  typedef boost::tuple<Intersection_type, Object> intersection_info;
//  typedef std::map<Halfedge_handle, intersection_info, Edge_comparator> Halfedge_intersection_map;
//  typedef std::map<Facet_handle, Halfedge_intersection_map> Facet_intersection_map;

//  Facet_intersection_map facet_intersection_map;    

//  for(std::vector<Primitive_id>::iterator it = intersected_triangles.begin();
//    it != intersected_triangles.end(); ++it)
//  {
//    Halfedge_intersection_map& halfedge_intersections = Facet_intersection_map[*it];
//    // Object intersection = CGAL::intersection(plane, *it);
//    typename Polyhedron::Facet::Halfedge_around_facet_circulator edge_it = (*it)->facet_begin();
//    do {
//      std::pair<Halfedge_intersection_map::iterator, bool> inserted =
//      facet_intersection_map.insert(std::make_pair(edge_it, intersection_info(NONE, Object())));
//      if(inserted.second) {
//        const Point& s = edge_it->vertex()->point();
//        const Point& t = edge_it->opposite()->vertex()->point();
//        Object intersection = CGAL::intersection(plane, Segment(s,t));
//        
//        const Point* i_point; 
//        if(!(i_point = object_cast<Point>(&intersection))) { 
//          CGAL_warning(false && "degenarate case");
//          continue;
//        }

//        inserted.first->second.get<0>() = POINT;
//        inserted.first->second.get<1>() = intersection;
//        polyline.push_back(*i_point);
//      }
//    } while(++edge_it != (*it)->facet_begin());

//    for()
//    std::vector<Point> polyline;
//    *out++ = polyline;
//    const Point& a = m_halfedge_handle->vertex()->point();
//    const Point& b = m_halfedge_handle->opposite()->vertex()->point();

//    
//    if(!(i_point = object_cast<Point>(&object))) 
//    { continue; } // continue in case of segment.
//    facet_intersection_map
//  }

//Facet_intersection_map_iterator
//get_neighbor(Facet_handle fh) {
//  Facet_intersection_map_iterator ret;
//  typename Polyhedron::Facet::Halfedge_around_facet_circulator edge_it = fh->facet_begin();
//  do {
//    ret = facet_intersection_map[edge_it->opposite()->facet()];
//    if(ret != facet_intersection_map.end()) { return ret; }
//  } while(++edge_it != (*it)->facet_begin());

//  return ret;
//}
