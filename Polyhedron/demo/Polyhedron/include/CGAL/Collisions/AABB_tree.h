#ifndef _AABB_TREE_
#define _AABB_TREE_

#include <list>
#include <stack>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include "AABB_node.h"
#include "knn.h"

CGAL_BEGIN_NAMESPACE

template <class Kernel, class Input, class PSC>
class AABB_tree
{
public:

  typedef typename Kernel::FT FT;
  typedef typename CGAL::Bbox_3 Bbox;
  typedef typename Kernel::Ray_3 Ray;
  typedef typename Kernel::Line_3 Line;
  typedef typename Kernel::Plane_3 Plane;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef typename Kernel::Segment_3 Segment;
  typedef typename Kernel::Triangle_3 Triangle;
  typedef typename Kernel::Iso_cuboid_3 Iso_cuboid;

  typedef AABB_node<Kernel,Input,PSC> Node;
  typedef typename Node::Point_with_input Point_with_input;

  // types for K nearest neighbors search structure
  typedef CNeighbor_search<Kernel> Neighbor_search;
  typedef typename Node::BConverter BConverter;

private:

  // set of input primitives (halfedge or face handles)
  std::vector<Input> m_data;

  // set of nodes
  std::vector<Node> m_nodes;

  // single root node
  Node *m_root;

public:
  // life cycle
  AABB_tree() {m_root = NULL;}
  ~AABB_tree()
  {
    cleanup();
  }

  void cleanup()
  {
    int n = std::max<int>(m_data.size(), 2);
    m_data.clear();
    if(m_root != NULL)
      delete [] m_root;
  }

  // build tree when Input = face_handle
  bool build_faces(PSC& psc)
  {
    cleanup();
    set_face_data(psc);
    m_root = NULL;
    if(!empty())
    {
      m_root = new Node[m_data.size()-1]();
      m_root->expand(psc, m_data.begin(),m_data.end(), m_data.size());
      return true;
    }
    return false;
  }

  void set_face_data(PSC& psc)
  {
    unsigned int nbf = psc.size_of_facets();
    m_data.reserve(nbf);
    typename PSC::Facet_iterator f;
    for(f = psc.facets_begin();
      f != psc.facets_end();
      f++)
      m_data.push_back(f);
  }

  bool empty()
  {
    return m_data.size() < 2; // TODO: change this requirement to < 1
  }



  bool furthest_intersection(const Ray& ray,
    const Point& from,
    Point_with_input& furthest)
  {
    std::vector<Point_with_input> ps;
    m_root->list_intersections(ray, ps, m_data.size());

    if(ps.size() == 0)
      return false;

    furthest = furthest_point_from(ps,from);
    return true;
  }


  bool furthest_intersection(const Segment& seg,
    const Point& from,
    Point_with_input& furthest)
  {
    std::vector<Point_with_input> ps;
    m_root->list_intersections(seg, ps, m_data.size());

    if(ps.size() == 0)
      return false;

    furthest = furthest_point_from(ps,from);
    return true;
  }

  bool furthest_intersection(const Line& line,
    const Point& from,
    Point_with_input& closest)
  {
    std::vector<Point_with_input> ps;
    m_root->list_intersections(line, ps, m_data.size());

    if(ps.size() == 0)
      return false;

    closest = furthest_point_from(ps,from);
    return true;
  }


  // -------------- Any query ----------------------//


  template<class PointInputContainer, class Pt>
  Point_with_input furthest_point_from(const PointInputContainer& intersections,
    const Pt& from)
  {
    typename PointInputContainer::const_iterator it = intersections.begin();
    Point_with_input furthest_point = *it;
    FT max_sqd = CGAL::squared_distance(from,furthest_point.first);
    it++;
    for(;it != intersections.end();it++)
    {
      FT sqd = CGAL::squared_distance(from,(*it).first);
      if(sqd > max_sqd)
      {
	furthest_point = *it;
	max_sqd = sqd;
      }
    }
    return furthest_point;
  }


};

CGAL_END_NAMESPACE

#endif // _AABB_TREE_
