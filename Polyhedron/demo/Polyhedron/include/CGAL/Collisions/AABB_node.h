#ifndef _AABB_NODE_
#define _AABB_NODE_

#include <CGAL/basic.h>
#include <CGAL/Bbox_3.h>
#include "meshing/PSC.h"
#include <vector>

#include "Ray_3_Bbox_3_do_intersect.h"
#include "Bbox_3_Bbox_3_do_intersect.h"
#include "Segment_3_Bbox_3_do_intersect.h"
#include "Plane_3_Bbox_3_do_intersect.h"
#include "Triangle_3_Bbox_3_do_intersect.h"
#include "Line_3_Bbox_3_do_intersect.h"
#include "Plucker_ray_3_Bbox_3_do_intersect.h"
#include "Sphere_3_Bbox_do_intersect.h"

namespace CGAL
{

  template <class Kernel, class Input, class PSC>
  class AABB_node
  {
  public:

    // type with fixed (double) floating point arithmetic
    typedef CGAL::Bbox_3 Bbox;

    // basic kernel object types
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Ray_3 Ray;
    typedef typename Kernel::Line_3 Line;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Plane_3 Plane;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::Sphere_3 Sphere;
    typedef typename Kernel::Segment_3 Segment;
    typedef typename Kernel::Triangle_3 Triangle;

    typedef typename PSC PSC;
    typedef typename AABB_node<Kernel,Input,PSC> Node;
    typedef typename std::vector<Input>::iterator Iterator;
    typedef typename std::pair<Point, Input> Point_with_input;
    //typedef typename PSC::Point_with_F_handle PointPSC_with_F;

    typedef typename PSC::Traits PSC_kernel;
    typedef typename PSC_kernel::Point_3 PSC_Point;
    typedef typename CGAL::Cartesian_converter<PSC_kernel, Kernel > Converter;
    typedef typename CGAL::Cartesian_converter<Kernel, PSC_kernel > BConverter;

    // Voronoi face query
    // segment is the Delaunay edge (dual of Voronoi face)
    // bbox is aabb of Voronoi face
    // vector contains the set of incident vertices
    typedef typename CGAL::Quadruple<Segment,Bbox,std::vector<Point>,std::vector<Point> > WF_query; 

    //Voronoi segment query
    //we need the Voronoi egde (generally a segment)
    //we need something to traverse the hierarchy - aabb tree already built
    //query will consist of the voronoi edge de and the delaunay facet abc
    typedef typename std::pair<Triangle,Segment> WS_query;
    typedef typename std::pair<Triangle,Ray> WR_query;
    typedef typename std::pair<Triangle,Line> WL_query;

    typedef typename std::pair<Line,int> Line_plucker;

    //ray types for ray-plucker
    struct Ray_plucker0 : public Ray { Ray_plucker0(const Ray& ray) : Ray(ray) {} };
    struct Ray_plucker1 : public Ray { Ray_plucker1(const Ray& ray) : Ray(ray) {} };
    struct Ray_plucker2 : public Ray { Ray_plucker2(const Ray& ray) : Ray(ray) {} };
    struct Ray_plucker3 : public Ray { Ray_plucker3(const Ray& ray) : Ray(ray) {} };
    struct Ray_plucker4 : public Ray { Ray_plucker4(const Ray& ray) : Ray(ray) {} };
    struct Ray_plucker5 : public Ray { Ray_plucker5(const Ray& ray) : Ray(ray) {} };
    struct Ray_plucker6 : public Ray { Ray_plucker6(const Ray& ray) : Ray(ray) {} };
    struct Ray_plucker7 : public Ray { Ray_plucker7(const Ray& ray) : Ray(ray) {} };


  private:

    // bounding box
    Bbox m_bbox;

    // children nodes
    // either pointing towards children (if the children are not leaves)
    // or pointing toward Input primitives (if the children are leaves)
    void *m_left_child;
    void *m_right_child;

    //for pointer traversal comment the above 2 lines and uncomment the next 3 lines
    //interchange traversal and traversal2, cleanup_recurse and cleanup_recurse2
    //expand and expand2
    //use the other aabb_node

    /*void *m_left_child;
    Node *m_next_if_no;
    bool is_leaf;*/


  public:

    // life cycle
    AABB_node()
    {
      m_left_child = m_right_child = NULL;
    }

    /*AABB_node()
    {
    m_left_child = m_next_if_no = NULL;
    is_leaf = false;
    }*/

    ~AABB_node() {}

  public:

    // deletes the whole subtree rooted at this node (except this node, of course).
    // node = number of primitives contained in this node.
    void cleanup_recurse(int nb_primitives)
    {
      switch(nb_primitives)
      {
      case 2:
	break;
      case 3:
	delete static_cast<Node*>(m_right_child);
	break;
      default:
	static_cast<Node*>(m_left_child)->cleanup_recurse(nb_primitives/2);
	static_cast<Node*>(m_right_child)->cleanup_recurse(nb_primitives - nb_primitives/2);
	delete static_cast<Node*>(m_left_child);
	delete static_cast<Node*>(m_right_child);
      }
    }

    void cleanup_recurse2(int nb_primitives)
    {
      //not done yet
    }

  private:

    // compute bbox for iterator range of input primitives
    Bbox bbox(const PSC& psc, Iterator a, Iterator b)
    {
      Bbox bbox = psc.compute_bbox(*a);
      for(++a; a != b; ++a)
	bbox = bbox + psc.compute_bbox(*a);
      return bbox;
    }

  public:
    // builds the tree by recursive expansion.
    // [a,b[ is the range of primitives to be added to the tree.
    // node is the length of this range.
    // psc is needed to compute the bboxes of the leaves.
    void expand(const PSC& psc, Iterator a, Iterator b, int range)
    {
      m_bbox = bbox(psc, a, b);

      // sort primitives along longest axis aabb
      sort_primitives(a, b);

      switch(range)
      {
      case 2:
	m_left_child = &(*a);
	m_right_child = &(*(++a));
	break;
      case 3:
	m_left_child = &(*a);
	//m_right_child = new Node();
	m_right_child = static_cast<Node*>(this)+1;
	static_cast<Node*>(m_right_child)->expand(psc, a+1, b, 2);	
	break;
      default:
	//m_left_child = new Node();
	//m_right_child = new Node();
	m_left_child = static_cast<Node*>(this)+1;
	m_right_child = (static_cast<Node*>(this))+(range/2);
	static_cast<Node*>(m_left_child)->expand(psc, a, a+range/2, range/2);
	static_cast<Node*>(m_right_child)->expand(psc, a+range/2, b, range - range/2);				
      }
    }

    // this builds the tree with pointers
    // builds the tree by recursive expansion.
    // [a,b[ is the range of primitives to be added to the tree.
    // node is the length of this range.
    // psc is needed to compute the bboxes of the leaves.
    void expand2(const PSC& psc, Iterator a, Iterator b, int range, Node* next = NULL)
    {
      m_bbox = bbox(psc, a, b);

      Node* m_right_child;
      // sort primitives along longest axis aabb
      sort_primitives(a, b);

      switch(range)
      {
      case 2:
	m_left_child = new Node();
	m_right_child = new Node();

	static_cast<Node*>(m_left_child)->m_left_child = &(*a);
	static_cast<Node*>(m_left_child)->is_leaf = true;
	static_cast<Node*>(m_left_child)->m_next_if_no = m_right_child;

	m_right_child->m_left_child = &(*(++a));
	m_right_child->is_leaf = true;
	m_right_child->m_next_if_no = next;

	m_next_if_no = next;
	break;

      case 3:
	m_left_child = new Node();
	m_right_child = new Node();

	static_cast<Node*>(m_left_child)->m_left_child = &(*a);
	static_cast<Node*>(m_left_child)->is_leaf = true;
	static_cast<Node*>(m_left_child)->m_next_if_no = m_right_child;

	(m_right_child)->expand(psc, a+1, b, 2, next);

	m_next_if_no = next;
	break;

      default:
	m_left_child = new Node();
	m_right_child = new Node();

	static_cast<Node*>(m_left_child)->expand(psc, a, a+range/2, range/2, m_right_child);
	(m_right_child)->expand(psc, a+range/2, b, range - range/2, next);

	m_next_if_no = next;
      }
    }

  private:
    void sort_primitives(Iterator a, Iterator b)
    {
      switch(longest_axis())
      {
      case 0: // sort along x
	std::sort(a,b,PSC::lower_x<Input>);
	break;
      case 1: // sort along y
	std::sort(a,b,PSC::lower_y<Input>);
	break;
      default: // sort along z
	std::sort(a,b,PSC::lower_z<Input>);
      }
    }

    int longest_axis()
    {
      FT max_size = std::max(xsize(),std::max(ysize(),zsize()));
      if(max_size == xsize())
	return 0; // axis along x
      if(max_size == ysize())
	return 1; // axis along y
      return 2; // axis along z
    }

    // size of bounding box along each axis
    FT xsize() { return m_bbox.xmax() - m_bbox.xmin(); }
    FT ysize() { return m_bbox.ymax() - m_bbox.ymin(); }
    FT zsize() { return m_bbox.zmax() - m_bbox.zmin(); }

  public:
    // -----------------------------------------------------------//
    // -----------------------DRAWING-----------------------------//
    // -----------------------------------------------------------//

    // draws the subtree rooted at this node.
    // nb_primitives = number of primitives contained in this node.
    void gl_draw_recurse(const PSC& psc, int nb_primitives)
    {
      traversal<Drawing_traits>(0, psc, nb_primitives);
    }

    class Drawing_traits
    {
    private:
      const PSC& m_psc;
    public:
      Drawing_traits(const PSC& psc) : m_psc(psc) {}
    public:
      bool go_further() { return true; }
      bool intersection(const int&, const Input& i)
      {
	gl_draw(m_psc.compute_bbox(i));
	return true;
      }
      bool do_intersect(const int&, // unused
	const Node& node)
      {
	gl_draw(node.m_bbox);
	return true;
      }
    };

    // draw bbox
    static void gl_draw(const Bbox& bb)
    {
      ::glBegin(GL_LINES);
      gl_draw_edge(bb.xmin(), bb.ymin(), bb.zmin(),
	bb.xmax(), bb.ymin(), bb.zmin());
      gl_draw_edge(bb.xmin(), bb.ymin(), bb.zmin(),
	bb.xmin(), bb.ymax(), bb.zmin());
      gl_draw_edge(bb.xmin(), bb.ymin(), bb.zmin(),
	bb.xmin(), bb.ymin(), bb.zmax());

      gl_draw_edge(bb.xmax(), bb.ymin(), bb.zmin(),
	bb.xmax(), bb.ymax(), bb.zmin());
      gl_draw_edge(bb.xmax(), bb.ymin(), bb.zmin(),
	bb.xmax(), bb.ymin(), bb.zmax());

      gl_draw_edge(bb.xmin(), bb.ymax(), bb.zmin(),
	bb.xmax(), bb.ymax(), bb.zmin());
      gl_draw_edge(bb.xmin(), bb.ymax(), bb.zmin(),
	bb.xmin(), bb.ymax(), bb.zmax());

      gl_draw_edge(bb.xmin(), bb.ymin(), bb.zmax(),
	bb.xmax(), bb.ymin(), bb.zmax());
      gl_draw_edge(bb.xmin(), bb.ymin(), bb.zmax(),
	bb.xmin(), bb.ymax(), bb.zmax());

      gl_draw_edge(bb.xmax(), bb.ymax(), bb.zmax(),
	bb.xmin(), bb.ymax(), bb.zmax());
      gl_draw_edge(bb.xmax(), bb.ymax(), bb.zmax(),
	bb.xmax(), bb.ymin(), bb.zmax());
      gl_draw_edge(bb.xmax(), bb.ymax(), bb.zmax(),
	bb.xmax(), bb.ymax(), bb.zmin());
      ::glEnd();
    }

    static void gl_draw_edge(double px, double py, double pz,
      double qx, double qy, double qz)
    {
      ::glVertex3d(px,py,pz); 
      ::glVertex3d(qx,qy,qz);
    }

    // -----------------------------------------------------------//
    // --------------------LISTING QUERIES------------------------//
    // -----------------------------------------------------------//

    // puts the intersection points with the query q in the container result.
    // This function will work for any type QueryType such that the following
    // two oracles are defined:

    //		static bool intersection(const QueryType& q, const typename Input& f, Container::value_type& p)
    //		static bool do_intersect(const QueryType& q, const Iso_cuboid& bbox)
    // [do NOT use a template parameter for the QueryType when defining such oracles]

    // nb_primitives = number of primitives contained in this node.
    template<class QueryType, class Container>
    void list_intersections(const QueryType& q,
      Container& result, 
      int nb_primitives)
    {
      traversal<Listing_traits<QueryType, Container> >(q, result, nb_primitives);
    }

    template<class QueryType, class Container>
    class Listing_traits
    {
      typedef typename Container::value_type Pt;
    private:
      Container& r;
    public:
      bool go_further()
      {
	return true;
      }
      Listing_traits(Container& result) : r(result) {}
      bool intersection(const QueryType& q, const Input& i)
      {
	Pt p;
	if(Node::intersection(q, i, p))
	{
	  r.push_back(p);
	  return true;
	}
	return false;
      }
      bool do_intersect(const QueryType& q, const Node& node)
      {
	return Node::do_intersect(q, node);
      }
    };

    // -----------------------------------------------------------//
    // --------------------PROJECTING QUERIES---------------------//
    // -----------------------------------------------------------//

    // puts the intersection points with the query q in the container result.
    // This function will work for any type QueryType such that the following
    // two oracles are defined:

    //		static bool intersection(const QueryType& q, const typename Input& f, Container::value_type& p)
    //		static bool do_intersect(const QueryType& q, const Iso_cuboid& bbox)
    // [do NOT use a template parameter for the QueryType when defining such oracles]

    // nb_nodes = number of primitives contained from this node.
    template<class QueryType, class Result>
    Result project(const QueryType& query,
      const Result& hint,
      int nb_nodes)
    {
      Result result = hint;
      traversal<Projecting_traits<QueryType, Result, Node::Projection_sphere<QueryType, Result>::type > >(query, result, nb_nodes);
      return result;
    }

    template<class QueryType, class Pt, class Sphere>
    class Projecting_traits
    {
    private:
      Pt& r;

    public:
      bool go_further()
      {
	return true;
      }
      Projecting_traits(Pt& result) : r(result) {}
      bool intersection(const QueryType& q, const Input& primitive)
      {
	Pt p;
	Point projected;
	if(Node::closest_intersection(Node::Projection_sphere<QueryType, Pt>::construct_sphere(q, r.first), primitive, projected))
	{
	  r = Pt(projected,primitive);
	  return true;
	}
	return false;
      }
      bool do_intersect(const QueryType& q, const Node& node)
      {
	return Node::do_intersect(Node::Projection_sphere<QueryType, Pt>::construct_sphere(q, r.first), node);
      }
    };

    // -----------------------------------------------------------//
    // -------------------DETECTING QUERIES-----------------------//
    // -----------------------------------------------------------//

    // tells whether the query q intersects the data.
    // This function will work for any type QueryType such that the following
    // two oracles are defined:

    //		static bool do_intersect(const QueryType& q, const typename Input& f)
    //		static bool do_intersect(const QueryType& q, const Iso_cuboid& bbox)
    // [do NOT use a template parameter for the QueryType when defining such oracles]

    // nb_primitives = number of primitives contained in this node.
    template<class QueryType>
    bool do_intersect(const QueryType& q, 
      int nb_primitives)
    {
      bool result = false;
      traversal<Detecting_traits<QueryType> >(q, result, nb_primitives);
      return result;
    }

    template<class QueryType>
    class Detecting_traits
    {
    private:
      bool& r;
    public:
      bool go_further()
      {
	return ! r;
      }
      Detecting_traits(bool& result): r(result) {}
      bool intersection(const QueryType& q, const Input& i)
      {
	return r = Node::do_intersect(q, i);
      }
      bool do_intersect(const QueryType& q, const Node& node)
      {
	return Node::do_intersect(q, node);
      }
    };

    // -----------------------------------------------------------//
    // --------------------COUNTING QUERIES-----------------------//
    // -----------------------------------------------------------//

    // counts how many times the query q intersects the data.
    // This function will work for any type QueryType such that the following
    // two oracles are defined:

    //		static bool do_intersect(const QueryType& q, const typename Input& f)
    //		static bool do_intersect(const QueryType& q, const Bbox& bbox)
    // [do NOT use a template parameter for the QueryType when defining such oracles]

    // nb_primitives = number of primitives contained in this node.
    template<class QueryType>
    unsigned int count_intersections(const QueryType& q, int nb_primitives)
    {
      unsigned int result = 0;
      traversal<Counting_traits<QueryType> >(q, result, nb_primitives);
      return result;
    }

    template<class QueryType>
    class Counting_traits
    {
    private:
      unsigned int& r;
    public:
      bool go_further()
      {
	return true;
      }
      Counting_traits(unsigned int& result): r(result) {}
      bool intersection(const QueryType& q, const Input& i)
      {
	if(Node::do_intersect(q, i))
	{
	  ++r;
	  return true;
	}
	return false;
      }
      bool do_intersect(const QueryType& q, const Node& node)
      {
	return Node::do_intersect(q, node);
      }
    };

    // -----------------------------------------------------------//
    // -----------------------LINE ORACLES-------------------------//
    // -----------------------------------------------------------//

    static bool intersection(const Line& line, 
      const typename PSC::Facet_handle& f, 
      Point& p)
    {
      return PSC::intersection(line, f, p);
    }

    static bool intersection(const Line& line, 
      const typename PSC::Facet_handle& f, 
      Point_with_input& p)
    {
      Point p_alone;
      if(Node::intersection(line, f, p_alone))
      {
	p = Point_with_input(p_alone, f);
	return true;
      }
      return false;
    }

    static bool do_intersect(const Line& line, 
      const typename PSC::Facet_handle& f)
    {
      return PSC::do_intersect(line, f);
    }

    static bool do_intersect(const Line& line,
      const Node& node)
    {
      return CGAL::do_intersect(line, node.m_bbox);
    }

    static bool do_intersect(const Line_plucker& line_plucker, 
      const typename PSC::Facet_handle& f)
    {
      return PSC::do_intersect(line_plucker.first, f);
    }

    static bool do_intersect(const Line_plucker& line_plucker,
      const Node& node)
    {

      const Bbox& bbox = node.m_bbox;
      const Line& line = line_plucker.first;

      const Point source = line.point(0); 

      const FT xmin = bbox.xmin()-source.x();
      const FT xmax = bbox.xmax()-source.x();
      const FT ymin = bbox.ymin()-source.y();
      const FT ymax = bbox.ymax()-source.y();
      const FT zmin = bbox.zmin()-source.z();
      const FT zmax = bbox.zmax()-source.z();

      const Vector direction1 = line.to_vector();

      const FT dirx = direction1.x();
      const FT diry = direction1.y();
      const FT dirz = direction1.z();

      const FT ftzero = (FT)0.0;

      switch(line_plucker.second)
      {
      case 0: //MMM
	//if it lies on correct side of the silhouette(6 edges) then it intersects
	//MM*
	if((-dirx*ymin + diry*xmax) > ftzero)//DH
	  return false;
	if((-dirx*ymax + diry*xmin) < ftzero)//BF
	  return false;

	//*MM
	if((diry*zmax - dirz*ymin) > ftzero)//HE
	  return false;
	if((diry*zmin - dirz*ymax) < ftzero)//CB
	  return false;

	//M*M
	if((dirx*zmax - dirz*xmin) > ftzero)//EF
	  return false;
	if((dirx*zmin - dirz*xmax) < ftzero)//DC
	  return false;

	return true;
	break;

      case 1: //MMP

	//MM*
	if((-dirx*ymin + diry*xmax) > ftzero)//DH
	  return false;
	if((-dirx*ymax + diry*xmin) < ftzero)//BF
	  return false;

	//*MP
	if((diry*zmax - dirz*ymax) > ftzero)//GF
	  return false;
	if((diry*zmin - dirz*ymin) < ftzero)//DA
	  return false;

	//M*P
	if((dirx*zmax - dirz*xmax) > ftzero)//HG
	  return false;
	if((dirx*zmin - dirz*xmin) < ftzero)//AB
	  return false;

	return true;

	break;
      case 2: //MPM

	//MP*
	if((-dirx*ymax + diry*xmax) < ftzero)//CG
	  return false;
	if((-dirx*ymin + diry*xmin) > ftzero)//AE
	  return false;

	//*PM
	if((diry*zmax - dirz*ymax) < ftzero)//GF
	  return false;
	if((diry*zmin - dirz*ymin) > ftzero)//DA
	  return false;

	//M*M
	if((dirx*zmax - dirz*xmin) > ftzero)//EF
	  return false;
	if((dirx*zmin - dirz*xmax) < ftzero)//DC
	  return false;

	return true;

	break;
      case 3: //MPP

	//MP*
	if((-dirx*ymax + diry*xmax) < ftzero)//CG
	  return false;
	if((-dirx*ymin + diry*xmin) > ftzero)//AE
	  return false;

	//*PP
	if((diry*zmax - dirz*ymin) < ftzero)//HE
	  return false;
	if((diry*zmin - dirz*ymax) > ftzero)//CB
	  return false;

	//M*P
	if((dirx*zmax - dirz*xmax) > ftzero)//HG
	  return false;
	if((dirx*zmin - dirz*xmin) < ftzero)//AB
	  return false;

	return true;

	break;
      case 4: //PMM
	//PM*
	if((-dirx*ymax + diry*xmax) > ftzero)//CG
	  return false;
	if((-dirx*ymin + diry*xmin) < ftzero)//AE
	  return false;

	//*MM
	if((diry*zmax - dirz*ymin) > ftzero)//HE
	  return false;
	if((diry*zmin - dirz*ymax) < ftzero)//CB
	  return false;

	//P*M
	if((dirx*zmax - dirz*xmax) < ftzero)//HG
	  return false;
	if((dirx*zmin - dirz*xmin) > ftzero)//AB
	  return false;

	return true;

	break;
      case 5: //PMP

	//PM*
	if((-dirx*ymax + diry*xmax) > ftzero)//CG
	  return false;
	if((-dirx*ymin + diry*xmin) < ftzero)//AE
	  return false;

	//*MP
	if((diry*zmax - dirz*ymax) > ftzero)//GF
	  return false;
	if((diry*zmin - dirz*ymin) < ftzero)//DA
	  return false;

	//P*P;
	if((dirx*zmax - dirz*xmin) < ftzero)//EF
	  return false;
	if((dirx*zmin - dirz*xmax) > ftzero)//DC
	  return false;

	return true;

	break;
      case 6: //PPM

	//PP*
	if((-dirx*ymin + diry*xmax) < ftzero)//DH
	  return false;
	if((-dirx*ymax + diry*xmin) > ftzero)//BF
	  return false;

	//*PM
	if((diry*zmax - dirz*ymax) < ftzero)//GF
	  return false;
	if((diry*zmin - dirz*ymin) > ftzero)//DA
	  return false;

	//P*M
	if((dirx*zmax - dirz*xmax) < ftzero)//HG
	  return false;
	if((dirx*zmin - dirz*xmin) > ftzero)//AB
	  return false;

	return true;

	break;
      case 7: //PPP

	//PP*
	if((-dirx*ymin + diry*xmax) < ftzero)//DH
	  return false;
	if((-dirx*ymax + diry*xmin) > ftzero)//BF
	  return false;

	//*PP
	if((diry*zmax - dirz*ymin) < ftzero)//HE
	  return false;
	if((diry*zmin - dirz*ymax) > ftzero)//CB
	  return false;

	//P*P
	if((dirx*zmax - dirz*xmin) < ftzero)//EF
	  return false;
	if((dirx*zmin - dirz*xmax) > ftzero)//DC
	  return false;

	return true;

	break;
      default:
	return false;
      }
      return true;
    }


    // -----------------------------------------------------------//
    // -----------------------SEGMENT ORACLES-------------------------//
    // -----------------------------------------------------------//


    static bool intersection(const Segment& segment, 
      const typename PSC::Facet_handle& f, 
      Point& p)
    {
      return PSC::intersection(segment, f, p);
    }

    static bool intersection(const Segment& segment, 
      const typename PSC::Facet_handle& f, 
      Point_with_input& p)
    {
      Point p_alone;
      if(Node::intersection(segment, f, p_alone))
      {
	p = Point_with_input(p_alone, f);
	return true;
      }
      return false;
    }

    static bool do_intersect(const Segment& segment, 
      const typename PSC::Facet_handle& f)
    {
      return PSC::do_intersect(segment, f);
    }

    static bool do_intersect(const Segment& segment,
      const Node& node)
    {
      return CGAL::do_intersect(segment, node.m_bbox);
    }

    // -----------------------------------------------------------//
    // -----------------------RAY ORACLES-------------------------//
    // -----------------------------------------------------------//

    static bool intersection(const Ray& ray, 
      const typename PSC::Facet_handle& f, 
      Point& p)
    {
      return PSC::intersection(ray, f, p);
    }

    static bool intersection(const Ray& ray, 
      const typename PSC::Facet_handle& f, 
      Point_with_input& p)
    {
      Point p_alone;
      if(Node::intersection(ray, f, p_alone))
      {
	p = Point_with_input(p_alone, f);
	return true;
      }
      return false;
    }

    static bool do_intersect(const Ray& ray, 
      const typename PSC::Facet_handle& f)
    {
      return PSC::do_intersect(ray, f);
    }

    static bool do_intersect(const Ray& ray,
      const Node& node)
    {
      return CGAL::do_intersect(ray, node.m_bbox);
    }


    //Ray_plucker based on http://pages.cpsc.ucalgary.ca/~blob/ps/jgt04.pdf
    //ray classification
    //no division
    static bool do_intersect(const Ray_plucker0& ray, const typename PSC::Facet_handle& f)
    {	return PSC::do_intersect(ray, f); }

    static bool do_intersect(const Ray_plucker1& ray, const typename PSC::Facet_handle& f)
    {	return PSC::do_intersect(ray, f); }

    static bool do_intersect(const Ray_plucker2& ray, const typename PSC::Facet_handle& f)
    {	return PSC::do_intersect(ray, f); }

    static bool do_intersect(const Ray_plucker3& ray, const typename PSC::Facet_handle& f)
    {	return PSC::do_intersect(ray, f); }

    static bool do_intersect(const Ray_plucker4& ray, const typename PSC::Facet_handle& f)
    {	return PSC::do_intersect(ray, f); }

    static bool do_intersect(const Ray_plucker5& ray, const typename PSC::Facet_handle& f)
    {	return PSC::do_intersect(ray, f); }

    static bool do_intersect(const Ray_plucker6& ray, const typename PSC::Facet_handle& f)
    {	return PSC::do_intersect(ray, f); }

    static bool do_intersect(const Ray_plucker7& ray, const typename PSC::Facet_handle& f)
    {	return PSC::do_intersect(ray, f); }

    static bool do_intersect(const Ray_plucker0& ray, const Node& node)
    { return CGAL::do_intersect_type_0<Kernel>(ray, node.m_bbox); }

    static bool do_intersect(const Ray_plucker1& ray, const Node& node)
    { return CGAL::do_intersect_type_1<Kernel>(ray, node.m_bbox); }

    static bool do_intersect(const Ray_plucker2& ray, const Node& node)
    { return CGAL::do_intersect_type_2<Kernel>(ray, node.m_bbox); }

    static bool do_intersect(const Ray_plucker3& ray, const Node& node)
    { return CGAL::do_intersect_type_3<Kernel>(ray, node.m_bbox); }

    static bool do_intersect(const Ray_plucker4& ray, const Node& node)
    { return CGAL::do_intersect_type_4<Kernel>(ray, node.m_bbox); }

    static bool do_intersect(const Ray_plucker5& ray, const Node& node)
    { return CGAL::do_intersect_type_5<Kernel>(ray, node.m_bbox); }

    static bool do_intersect(const Ray_plucker6& ray, const Node& node)
    { return CGAL::do_intersect_type_6<Kernel>(ray, node.m_bbox); }

    static bool do_intersect(const Ray_plucker7& ray, const Node& node)
    { return CGAL::do_intersect_type_7<Kernel>(ray, node.m_bbox); }


    // -----------------------------------------------------------//
    // ----------------------PLANE ORACLES------------------------//
    // -----------------------------------------------------------//

    static bool intersection(const Plane& plane, 
      const typename PSC::Halfedge_handle& f, 
      Point& p)
    {
      return PSC::intersection(plane, f, p);
    }

    static bool intersection(const Plane& plane, 
      const typename PSC::Halfedge_handle& f, 
      Point_with_input& p)
    {
      Point p_alone;
      if(Node::intersection(plane, f, p_alone))
      {
	p = Point_with_input(p_alone, f);
	return true;
      }
      return false;
    }

    static bool do_intersect(const Plane& plane, const Node& node)
    {
      return CGAL::do_intersect(plane, node.m_bbox);
    }

    // -----------------------------------------------------------//
    // -------------------PROJECTION ORACLES----------------------//
    // -----------------------------------------------------------//

    // Generic template, which should be specialized for each use.
    template <class QueryType, class ProjectionType>
    class Projection_sphere {};

    template<>
    class Projection_sphere<Point, Point_with_input> 
    {
    public:
      typedef typename Kernel::Sphere_3	type;
      static type construct_sphere(const Point& p, const Point& q)
      {
	//Converter c;
	//return type(p, CGAL::squared_distance(p, c(q)));
	return type(p, CGAL::squared_distance(p,q));
      }
    };

    static bool closest_intersection(const Sphere& sphere, 
      Input input,
      Point& projected)
    {
      Converter c;
      BConverter bc;
      PSC_Point psc_point = bc(projected);
      bool success = PSC::closest_intersection_from_sphere_center(bc(sphere), input, psc_point);
      projected = c(psc_point);
      return success;
      //return PSC::closest_intersection_from_sphere_center(bc(sphere), input, bc(projected));
    }

    static bool do_intersect(const typename Kernel::Sphere_3& sphere, const Node& node)
    {
      return CGAL::do_intersect(sphere, node.m_bbox);
    }

    // -----------------------------------------------------------//
    // ------------------FINITE VORONOI SEGMENT ORACLES--------------//
    // -----------------------------------------------------------//


    // intersect bbox of query with bbox of node
    static bool do_intersect(const WS_query& query,
      const Node& node)
    {
      //return CGAL::do_intersect(query.second, node.m_bbox);
      return true;
    }

    static bool do_intersect(const WS_query& query,
      const typename PSC::Facet_handle f)
    {

      const Triangle& abc = query.first;
      const Point& a = abc.vertex(0);
      const Point& b = abc.vertex(1);
      const Point& c = abc.vertex(2);

      Converter convert;

      const Point& p = convert(f->halfedge()->vertex()->point());
      const Point& q = convert(f->halfedge()->next()->vertex()->point());
      const Point& r = convert(f->halfedge()->next()->next()->vertex()->point());

      int pos[3][3];
      int exored[3];
      int summed[3];
      int orpos[3];
      int andofexor,orofexor;

      //can be -1, 0, 1
      pos[0][0] = CGAL::compare_distance_to_point(p,a,b);  //(ap < bp);
      pos[1][0] = CGAL::compare_distance_to_point(q,a,b);  //(aq < bq);
      pos[2][0] = CGAL::compare_distance_to_point(r,a,b);  //(ar < br);

      //exit if all three point are on one side
      summed[0] = abs(pos[0][0] + pos[1][0] + pos[2][0]);
      if(summed[0] == 3)
	return false;

      pos[0][1] = CGAL::compare_distance_to_point(p,b,c);  //(bp < cp);
      pos[1][1] = CGAL::compare_distance_to_point(q,b,c);  //(bq < cq);
      pos[2][1] = CGAL::compare_distance_to_point(r,b,c);  //(br < cr);
      summed[1] = abs(pos[0][1] + pos[1][1] + pos[2][1]);
      if(summed[1] == 3)
	return false;

      pos[0][2] = CGAL::compare_distance_to_point(p,c,a);  //(cp < ap);
      pos[1][2] = CGAL::compare_distance_to_point(q,c,a);  //(cq < aq);
      pos[2][2] = CGAL::compare_distance_to_point(r,c,a);  //(cr < ar);
      summed[2] = abs(pos[0][2] + pos[1][2] + pos[2][2]);
      if(summed[2] == 3)
	return false;

      const Segment& de = query.second;

      const Point& d = de.source();
      const Point& e = de.target();


      CGAL::Orientation pqrd = CGAL::orientation(p,q,r,d);
      if(pqrd == CGAL::orientation(p,q,r,e) && pqrd != CGAL::COLLINEAR)
	return false;

      orpos[0] = pos[0][0] | pos[0][1];
      orpos[1] = pos[1][0] | pos[1][1];
      orpos[2] = pos[2][0] | pos[2][1];
      if((orpos[0] & orpos[1] & orpos[2])==0) //orpos can be 0 or -1, hence if "and of all the orpos" is zero then we have a zero
	return true;
      if((summed[0] | summed[1] | summed[2])>=2) // atleast one of them is 2 - they can either be 0, 1 or 2
      {
	return false;
      }

      exored[0] = pos[0][0] ^ pos[1][0] ^ pos[2][0];
      exored[1] = pos[0][1] ^ pos[1][1] ^ pos[2][1];
      exored[2] = pos[0][2] ^ pos[1][2] ^ pos[2][2];
      andofexor = exored[0] & exored[1] & exored[2];
      orofexor = exored[0] | exored[1] | exored[2];
      //if exoreds are 1 1 1 then return true
      //if exoreds are -1 -1 -1 then return true
      if(andofexor == -1 || (andofexor == 1 && orofexor == 1))
	return true;

      //if exoreds are -2 -2 -2 then return true - dont think this is worth it


      return PSC::do_intersect(de,f);

    }

    // intersect bbox of query with bbox of node
    static bool do_intersect(const WR_query& query,
      const Node& node)
    {
      //return CGAL::do_intersect(query.second, node.m_bbox);
      return true;
    }


    static bool do_intersect(const WR_query& query,
      const typename PSC::Facet_handle f)
    {

      const Triangle& abc = query.first;
      const Point& a = abc.vertex(0);
      const Point& b = abc.vertex(1);
      const Point& c = abc.vertex(2);

      Converter convert;

      const Point& p = convert(f->halfedge()->vertex()->point());
      const Point& q = convert(f->halfedge()->next()->vertex()->point());
      const Point& r = convert(f->halfedge()->next()->next()->vertex()->point());


      int pos[3][3];
      int exored[3];
      int summed[3];
      int orpos[3];
      int andofexor,orofexor;

      //can be -1, 0, 1
      pos[0][0] = CGAL::compare_distance_to_point(p,a,b);  //(ap < bp);
      pos[1][0] = CGAL::compare_distance_to_point(q,a,b);  //(aq < bq);
      pos[2][0] = CGAL::compare_distance_to_point(r,a,b);  //(ar < br);

      //exit if all three point are on one side
      summed[0] = abs(pos[0][0] + pos[1][0] + pos[2][0]);
      if(summed[0] == 3)
	return false;

      pos[0][1] = CGAL::compare_distance_to_point(p,b,c);  //(bp < cp);
      pos[1][1] = CGAL::compare_distance_to_point(q,b,c);  //(bq < cq);
      pos[2][1] = CGAL::compare_distance_to_point(r,b,c);  //(br < cr);
      summed[1] = abs(pos[0][1] + pos[1][1] + pos[2][1]);
      if(summed[1] == 3)
	return false;

      pos[0][2] = CGAL::compare_distance_to_point(p,c,a);  //(cp < ap);
      pos[1][2] = CGAL::compare_distance_to_point(q,c,a);  //(cq < aq);
      pos[2][2] = CGAL::compare_distance_to_point(r,c,a);  //(cr < ar);
      summed[2] = abs(pos[0][2] + pos[1][2] + pos[2][2]);
      if(summed[2] == 3)
	return false;

      const Ray& ray = query.second;


      if(!CGAL::do_intersect(ray,Plane(p,q,r)))
	return false;

      orpos[0] = pos[0][0] | pos[0][1];
      orpos[1] = pos[1][0] | pos[1][1];
      orpos[2] = pos[2][0] | pos[2][1];
      if((orpos[0] & orpos[1] & orpos[2])==0) //orpos can be 0 -1 or 1, hence if and of all the orpos is zero then we have a zero
	return true;
      if((summed[0] | summed[1] | summed[2])>=2) // atleast one of them is 2 - they can either be 0, 1 or 2
      {
	return false;
      }

      exored[0] = pos[0][0] ^ pos[1][0] ^ pos[2][0];
      exored[1] = pos[0][1] ^ pos[1][1] ^ pos[2][1];
      exored[2] = pos[0][2] ^ pos[1][2] ^ pos[2][2];
      andofexor = exored[0] & exored[1] & exored[2];
      orofexor = exored[0] | exored[1] | exored[2];
      //if exoreds are 1 1 1 then return true
      //if exoreds are -1 -1 -1 then return true
      if(andofexor == -1 || (andofexor == 1 && orofexor == 1))
	return true;

      //if exoreds are -2 -2 -2 then return true - dont think this is worth it

      return PSC::do_intersect(ray,f);


    }

    // intersect bbox of query with bbox of node
    static bool do_intersect(const WL_query& query,
      const Node& node)
    {
      //return CGAL::do_intersect(query.second, node.m_bbox);
      return true;
    }

    static bool do_intersect(const WL_query& query,
      const typename PSC::Facet_handle f)
    {

      const Triangle& abc = query.first;
      const Point& a = abc.vertex(0);
      const Point& b = abc.vertex(1);
      const Point& c = abc.vertex(2);

      Converter convert;

      const Point& p = convert(f->halfedge()->vertex()->point());
      const Point& q = convert(f->halfedge()->next()->vertex()->point());
      const Point& r = convert(f->halfedge()->next()->next()->vertex()->point());

      int pos[3][3];
      int exored[3];
      int summed[3];
      int orpos[3];
      int andofexor,orofexor;


      //can be -1, 0, 1
      pos[0][0] = CGAL::compare_distance_to_point(p,a,b);  //(ap < bp);
      pos[1][0] = CGAL::compare_distance_to_point(q,a,b);  //(aq < bq);
      pos[2][0] = CGAL::compare_distance_to_point(r,a,b);  //(ar < br);

      //exit if all three point are on one side
      summed[0] = abs(pos[0][0] + pos[1][0] + pos[2][0]);
      if(summed[0] == 3)
	return false;

      pos[0][1] = CGAL::compare_distance_to_point(p,b,c);  //(bp < cp);
      pos[1][1] = CGAL::compare_distance_to_point(q,b,c);  //(bq < cq);
      pos[2][1] = CGAL::compare_distance_to_point(r,b,c);  //(br < cr);
      summed[1] = abs(pos[0][1] + pos[1][1] + pos[2][1]);
      if(summed[1] == 3)
	return false;

      pos[0][2] = CGAL::compare_distance_to_point(p,c,a);  //(cp < ap);
      pos[1][2] = CGAL::compare_distance_to_point(q,c,a);  //(cq < aq);
      pos[2][2] = CGAL::compare_distance_to_point(r,c,a);  //(cr < ar);
      summed[2] = abs(pos[0][2] + pos[1][2] + pos[2][2]);
      if(summed[2] == 3)
	return false;

      Line line = query.second;



      if(!CGAL::do_intersect(line,Plane(p,q,r)))
	return false;

      orpos[0] = pos[0][0] | pos[0][1];
      orpos[1] = pos[1][0] | pos[1][1];
      orpos[2] = pos[2][0] | pos[2][1];
      if((orpos[0] & orpos[1] & orpos[2])==0) //orpos can be 0 -1 or 1, hence if and of all the orpos is zero then we have a zero
	return true;
      if((summed[0] | summed[1] | summed[2])>=2) // atleast one of them is 2 - they can either be 0, 1 or 2
      {
	return false;
      }

      exored[0] = pos[0][0] ^ pos[1][0] ^ pos[2][0];
      exored[1] = pos[0][1] ^ pos[1][1] ^ pos[2][1];
      exored[2] = pos[0][2] ^ pos[1][2] ^ pos[2][2];
      andofexor = exored[0] & exored[1] & exored[2];
      orofexor = exored[0] | exored[1] | exored[2];
      //if exoreds are 1 1 1 then return true
      //if exoreds are -1 -1 -1 then return true
      if(andofexor == -1 || (andofexor == 1 && orofexor == 1))
	return true;

      //if exoreds are -2 -2 -2 then return true - dont think this is worth it

      return PSC::do_intersect(line,f);

    }

    // -----------------------------------------------------------//
    // ------------------FINITE VORONOI FACE ORACLES--------------//
    // -----------------------------------------------------------//

    static bool intersection(const WF_query& query,
      const typename PSC::Halfedge_handle& he,
      Point& z)
    {
      if(! Node::do_intersect(query,he))
	return false;

      // primal edge ab
      const Segment& primal_edge = query.first;
      Point m = CGAL::midpoint(primal_edge.source(), primal_edge.target());
      Plane plane(m, Vector(primal_edge.source(), primal_edge.target()));

      // halfedge pq
      Converter convert;
      Point p = convert(he->vertex()->point());
      Point q = convert(he->opposite()->vertex()->point());

      CGAL::Object o = CGAL::intersection(plane, Segment(p,q));
      if(CGAL::assign(z, o))
	return true;
      return false;
    }

    static bool intersection(const WF_query& query,
      const typename PSC::Halfedge_handle& he,
      Point_with_input& z)
    {
      Point p_alone;
      if(Node::intersection(query, he, p_alone))
      {
	z = Point_with_input(p_alone, he);
	return true;
      }
      return false;
    }

    // intersect bbox of query with bbox of node
    static bool do_intersect(const WF_query& query,
      const Node& node)
    {
      return CGAL::do_overlap(query.second,node.m_bbox);
    }

    //original version
    static bool do_intersect2(const WF_query& query,
      const typename PSC::Halfedge_handle he)
    {
      // primal edge ab
      const Segment& primal_edge = query.first;
      const Point& a = primal_edge.source();
      const Point& b = primal_edge.target();

      // halfedge pq
      Converter convert;
      Point p = convert(he->vertex()->point());
      Point q = convert(he->opposite()->vertex()->point());

      Vector pq = q - p;
      Vector ab = b - a;
      FT pqab = pq * ab;
      // case when the crease of the PSC and the 
      // edge of the mesh are orthogonal:
      if(pqab == (FT)0.0) 
	return false;

      // ratio is such that m (= p + ratio * pq) is the intersection
      // point of the mediator plane of [a b] with (p q)
      FT ratio = (((a - p) + (b - p))*ab) / ((FT)2.0 * pqab);
      // case when m does not belong to the segment [p q]
      if(ratio > (FT)1.0 || ratio < (FT)0.0)
	return false;

      Point m = p + ratio * pq;

      const std::vector<Point>& incidents = query.third;
      std::vector<Point>::const_iterator it;
      for(it = incidents.begin(); 
	it != incidents.end();
	++it)
      {
	// we have d(m,a) == d(m,b) by construction of m.
	// we now check that m belongs to the cells of a and b rather
	// than to the ones of their neighbors.
	if(CGAL::has_smaller_distance_to_point(m,*it,a))
	  return false;
      }
      return true;
    }

    //new version
    static bool do_intersect(const WF_query& query,
      const typename PSC::Halfedge_handle he)
    {

      // primal edge ab
      const Segment& primal_edge = query.first;
      const Point& a = primal_edge.source();
      const Point& b = primal_edge.target();

      // halfedge pq
      Converter convert;
      Point p = convert(he->vertex()->point());
      Point q = convert(he->opposite()->vertex()->point());

      const std::vector<Point>& polygon = query.fourth;
      unsigned int m = 0;

      //assuming atleast 3 points
      //assuming closed convex polygon
      //assuming planarity of points of polygon
      Vector pq = q - p;
      Vector ab = b - a;
      FT pqab = pq * ab;

      FT ratio = (((a - p) + (b - p))*ab) / ((FT)2.0 * pqab);
      if(ratio > (FT)1.0 || ratio < (FT)0.0)
	return false;

      CGAL::Orientation o = CGAL::orientation(polygon[polygon.size()-1],polygon[0],p,q);
      if(o==CGAL::COPLANAR)
	return true;
      CGAL::Orientation o2;
      unsigned int i;
      unsigned int sizem1 = polygon.size() - 1;
      for(i=0;i<sizem1;i++)
      {
	o2 = CGAL::orientation(polygon[i],polygon[i+1],p,q);
	if(o2==o)
	  continue;
	if(o2==CGAL::COPLANAR)
	  return true;
	return false;
      }

      return true;
    }


    // -----------------------------------------------------------//
    // ----------------------GENERAL QUERY------------------------//
    // -----------------------------------------------------------//

    // general traversal query, the traits class allows to use it for the various 
    // traversal methods we need: listing, counting, detecting intersections, drawing the boxes.

    // TODO: document the traits class requirements...
    template<class Traits, class QueryType, class ResultType>
    void traversal(const QueryType& query,
      ResultType& result,
      int nb_primitives)
    {
      Traits traits(result);
      bool left_test;
      switch(nb_primitives)
      {
      case 2: // comment
	left_test = traits.intersection(query, *static_cast<Input*>(m_left_child));
	if((! left_test) || (traits.go_further()))
	  traits.intersection(query, *static_cast<Input*>(m_right_child));
	break;
      case 3: // comment
	left_test = traits.intersection(query, *static_cast<Input*>(m_left_child));
	if((! left_test) || (traits.go_further()))
	  if(traits.do_intersect(query, *static_cast<Node*>(m_right_child)))
	    static_cast<Node*>(m_right_child)->traversal<Traits>(query, result, 2);
	break;
      default: // comment
	if(traits.do_intersect(query, *static_cast<Node*>(m_left_child)))
	{
	  static_cast<Node*>(m_left_child)->traversal<Traits>(query, result, nb_primitives/2);
	  if(traits.go_further())
	    if(traits.do_intersect(query, *static_cast<Node*>(m_right_child)))
	      static_cast<Node*>(m_right_child)->traversal<Traits>(query, result, nb_primitives - nb_primitives/2);
	}
	else
	  if(traits.do_intersect(query, *static_cast<Node*>(m_right_child)))
	    static_cast<Node*>(m_right_child)->traversal<Traits>(query, result, nb_primitives - nb_primitives/2);
      }
    }


    //traversal using skip pointers - prevents recursion overhead
    template<class Traits, class QueryType, class ResultType>
    void traversal2(const QueryType& query,
      ResultType& result,
      int nb_primitives)
    {
      Traits traits(result);
      const bool go_further = traits.go_further();
      Node *curr = this;
      while(curr!=NULL)
      {
	if(curr->is_leaf)
	{
	  if(traits.intersection(query, *static_cast<Input*>(curr->m_left_child)))
	  {
	    if(go_further)
	    {
	      curr = (curr->m_next_if_no);
	    }
	    else
	      curr = NULL;
	  }
	  else
	  {
	    curr = (curr->m_next_if_no);
	  }
	}
	else
	{
	  if(traits.do_intersect(query, *curr))
	  {
	    curr = static_cast<Node*>(curr->m_left_child);
	  }
	  else
	    curr = (curr->m_next_if_no);
	}
      }
    }
  };

} // end namespace VMESH

#endif // _AABB_NODE_

