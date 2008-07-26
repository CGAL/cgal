#ifndef _AABB_TREE_
#define _AABB_TREE_

#include <list>
#include <stack>
#include <CGAL/basic.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include "meshing/PSC.h"
#include "AABB_node.h"
#include "knn.h"

namespace CGAL
{
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

    typedef typename AABB_node<Kernel,Input,PSC> Node;
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

    Neighbor_search m_knn;

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

    void build_knn_from_all_vertices(PSC& psc)
    {
      std::list<Point> points;
      PSC::Vertex_iterator v;
      for(v = psc.vertices_begin();
	v!= psc.vertices_end();
	v++)
	points.push_back(Point(v->point().x(),v->point().y(),v->point().z()));
      m_knn.init(points);
    }

    void build_knn_from_crease_vertices(PSC& psc)
    {
      std::list<Point> points;
      PSC::Vertex_iterator v;
      for(v = psc.vertices_begin();
	v!= psc.vertices_end();
	v++)
      {
	if(psc.is_feature(v))
	  points.push_back(Point(v->point().x(),v->point().y(),v->point().z()));
      }
      m_knn.init(points);
    }

    // build tree when Input = face_handle
    bool build_faces(PSC& psc)
    {
      cleanup();
      build_knn_from_all_vertices(psc);
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
      PSC::Facet_iterator f;
      for(f = psc.facets_begin();
	f != psc.facets_end();
	f++)
	m_data.push_back(f);
    }

    bool empty()
    {
      return m_data.size() < 2; // TODO: change this requirement to < 1
    }

    // build tree when Input = halfedge_handle
    bool build_creases(PSC& psc)
    {
      cleanup();
      build_knn_from_crease_vertices(psc);
      set_crease_data(psc);
      m_root = NULL;
      if(!empty())
      {
	m_root = new Node[m_data.size()-1]();
	m_root->expand(psc, m_data.begin(),m_data.end(), m_data.size());
	return true;
      }
      return false;
    }

    void set_crease_data(PSC& psc)
    {
      unsigned int nbc = psc.nb_crease_edges();
      m_data.reserve(nbc);

      std::list<PSC::Crease>& creases = psc.creases();
      std::list<PSC::Crease>::iterator cit;
      for(cit = creases.begin(); cit != creases.end(); cit++)
      {
	PSC::Crease::iterator heit;
	for(heit = cit->begin(); heit != cit->end(); heit++)
	  m_data.push_back(*heit);
      }
    }

    // --------------- RAY QUERIES ---------------------- // 

    bool closest_intersection(const Ray& ray,
      const Point& from,
      Point_with_input& closest)
    {
      std::vector<Point_with_input> ps;
      m_root->list_intersections(ray, ps, m_data.size());

      if(ps.size() == 0)
	return false;

      closest = closest_point_from(ps,from);
      return true;
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

    unsigned int nb_intersections(const Ray& ray)
    {
      return m_root->count_intersections(ray, m_data.size());		
    }

    unsigned int nb_intersections_plucker(const Ray& ray)
    {
      const Vector direction = ray.to_vector();
      int ray_type = 0;
      if(direction.x()> (FT)0.0)
	ray_type = ray_type+4;
      if(direction.y()> (FT)0.0)
	ray_type = ray_type+2;
      if(direction.z()> (FT)0.0)
	ray_type = ray_type+1;

      switch(ray_type)
      {
      case 0: return m_root->count_intersections(Node::Ray_plucker0(ray), m_data.size());
      case 1: return m_root->count_intersections(Node::Ray_plucker1(ray), m_data.size());
      case 2: return m_root->count_intersections(Node::Ray_plucker2(ray), m_data.size());
      case 3: return m_root->count_intersections(Node::Ray_plucker3(ray), m_data.size());
      case 4: return m_root->count_intersections(Node::Ray_plucker4(ray), m_data.size());
      case 5: return m_root->count_intersections(Node::Ray_plucker5(ray), m_data.size());
      case 6: return m_root->count_intersections(Node::Ray_plucker6(ray), m_data.size());
      default: return m_root->count_intersections(Node::Ray_plucker7(ray), m_data.size());
      }	
    }

    // --------------- SEGMENT QUERIES ---------------------- // 


    bool closest_intersection(const Segment& seg,
      const Point& from,
      Point_with_input& closest)
    {
      std::vector<Point_with_input> ps;
      m_root->list_intersections(seg, ps, m_data.size());

      if(ps.size() == 0)
	return false;

      closest = closest_point_from(ps,from);
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

    unsigned int nb_intersections(const Segment& s)
    {
      return m_root->count_intersections(s, m_data.size());
    }

    // --------------- LINE QUERIES ---------------------- // 

    unsigned int nb_intersections(const Line& line)
    {
      return m_root->count_intersections(line,m_data.size());
    }

    unsigned int nb_intersections_plucker(const Line& line)
    {
      const Vector direction = line.to_vector();
      int line_type = 0;
      if(direction.x()> (FT)0.0)
	line_type = line_type+4;
      if(direction.y()> (FT)0.0)
	line_type = line_type+2;
      if(direction.z()> (FT)0.0)
	line_type = line_type+1;
      return m_root->count_intersections(Node::Line_plucker(line,line_type),m_data.size());
    }

    bool closest_intersection(const Line& line,
      const Point& from,
      Point_with_input& closest)
    {
      std::vector<Point_with_input> ps;
      m_root->list_intersections(line, ps, m_data.size());

      if(ps.size() == 0)
	return false;

      closest = closest_point_from(ps,from);
      return true;
    }



    // -------------- PLANE QUERIES -------------------//
    bool closest_intersection_plane(const Plane& plane,
      const Point& from,
      Point_with_input& closest)
    {
      std::vector<Point_with_input> ps;
      m_root->list_intersections(plane, ps, m_data.size());

      if(ps.size() == 0)
	return false;

      closest = closest_point_from(ps, from);
      return true;
    }

    unsigned int nb_intersections(const Plane& plane)
    {
      return m_root->count_intersections(plane, m_data.size());
    }

    // -------------- PROJECTION QUERIES ---------------//

    Point_with_input project(const Point &point)
    {
      Point_with_input hint;
      hint.first = m_knn.nearest_point(point);
      //std::cout<<"First closest: ("<<hint.first.x()<<" , "<<hint.first.y()<<" , "<<hint.first.z()<<")"<<std::endl;
      return m_root->project(point, hint, m_data.size());
    }

    //--------------VORONOI SEGMENT QUERIES-------------//
    unsigned int nb_intersections_voronoi_segment(const Triangle& primal_facet,
      const Segment& voronoi_segment)
    {
      //if(empty()) return 0;
      //Bbox poly_bb = get_bb_polygon(polygon);
      return m_root->count_intersections(Node::WS_query(primal_facet,voronoi_segment), 
	m_data.size());
    }

    unsigned int nb_intersections_voronoi_ray(const Triangle& primal_facet,
      const Ray& voronoi_ray)
    {
      //if(empty()) return 0;
      //Bbox poly_bb = get_bb_polygon(polygon);
      return m_root->count_intersections(Node::WR_query(primal_facet,voronoi_ray), 
	m_data.size());
    }

    unsigned int nb_intersections_voronoi_line(const Triangle& primal_facet,
      const Line& voronoi_line)
    {
      //if(empty()) return 0;
      //Bbox poly_bb = get_bb_polygon(polygon);
      return m_root->count_intersections(Node::WL_query(primal_facet,voronoi_line), 
	m_data.size());
    }

    // ------------- POLYGON QUERIES ------------------ //

    // NEW VERSION (without triangulating the polygon)

    template<typename PointContainer>
    bool closest_intersection_polygon(const Segment& primal_segment,
      const PointContainer& polygon,
      const PointContainer& incidents,
      Point_with_input& closest)
    {
      Point from = CGAL::midpoint(primal_segment.source(), primal_segment.target());
      std::vector<Point_with_input> intersections;
      if(intersections_polygon(primal_segment, polygon, incidents, intersections))
      {
	closest = closest_point_from(intersections, from);
	return true;
      }
      else return false;
    }

    template<typename PointContainer>
    bool furthest_intersection_polygon(const Segment& primal_segment,
      const PointContainer& polygon,
      const PointContainer& incidents,
      Point_with_input& furthest)
    {
      Point from = CGAL::midpoint(primal_segment.source(), primal_segment.target());
      std::vector<Point_with_input> intersections;
      if(intersections_polygon(primal_segment, polygon, incidents, intersections))
      {
	furthest = furthest_point_from(intersections, from);
	return true;
      }
      else return false;
    }

    template<typename PointContainer, typename PointInputContainer>
    bool intersections_polygon(const Segment& primal_segment,
      const PointContainer& polygon,
      const PointContainer& incidents,
      PointInputContainer& intersections)
    {
      if(empty()) return false;
      Bbox poly_bb = get_bb_polygon(polygon);
      m_root->list_intersections(Node::WF_query(primal_segment, poly_bb, incidents,polygon),
	intersections, m_data.size());
      if(intersections.size() == 0)
	return false;
      else return true;
    }

    template <typename PointContainer>
    bool do_intersect_dual_face(const Segment& primal_edge,
      const PointContainer& incident_vertices,
      const PointContainer& voronoi_face)
    {
      if(empty()) return false;
      Bbox bbox = get_bb_polygon(voronoi_face);

      typedef typename CGAL::Quadruple<Segment,Bbox,PointContainer,PointContainer> WF_query;
      return m_root->do_intersect(WF_query(primal_edge,bbox,incident_vertices, voronoi_face),m_data.size()); 
    }

    template <typename PointContainer>
    unsigned int nb_intersect_dual_face(const Segment& primal_edge,
      const PointContainer& incident_vertices,
      const PointContainer& voronoi_face)
    {
      if(empty()) return 0;
      Bbox bbox = get_bb_polygon(voronoi_face);

      typedef typename CGAL::Quadruple<Segment,Bbox,PointContainer,PointContainer> WF_query;
      return m_root->count_intersections(WF_query(primal_edge,bbox,incident_vertices, voronoi_face),m_data.size()); 
    }

    template<typename PointContainer>
    bool do_intersect_finite_voronoi_face_2(const PointContainer& polygon)
    {
      if(empty()) return false;
      if(polygon.size() < 3) return false;
      return m_root->do_intersect(polygon, m_data.size());
    }

    template <typename PointContainer>
    static Bbox get_bb_polygon(const PointContainer& points)
    {
      PointContainer::const_iterator pit = points.begin();
      Bbox bbox = (*pit).bbox(); // bbox of one point exists
      for(; pit != points.end(); ++pit)
	bbox = bbox + (*pit).bbox();
      return bbox;
    }


    // -------------- Any query ----------------------//

    template<class PointInputContainer, class Pt>
    Point_with_input closest_point_from(const PointInputContainer& intersections,
      const Pt& from)
    {
      PointInputContainer::const_iterator it = intersections.begin();
      Point_with_input closest_point = *it;
      FT min_sqd = CGAL::squared_distance(from, closest_point.first);
      it++;
      for(;it != intersections.end();it++)
      {
	FT sqd = CGAL::squared_distance(from,(*it).first);
	if(sqd < min_sqd)
	{
	  closest_point = *it;
	  min_sqd = sqd;
	}
      }
      return closest_point;
    }

    template<class PointInputContainer, class Pt>
    Point_with_input furthest_point_from(const PointInputContainer& intersections,
      const Pt& from)
    {
      PointInputContainer::const_iterator it = intersections.begin();
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


    // RENDERING
    void gl_draw(PSC& psc)
    {
      m_root->gl_draw_recurse(psc, m_data.size());
    }
  };

}

#endif // _AABB_TREE_
