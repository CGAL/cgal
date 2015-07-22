/***********************************************************************

Takes a list of points and returns a list of segments corresponding to
the weighted Alpha Shape.

************************************************************************/

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

#include <CGAL/Weighted_point.h>
#include <CGAL/Weighted_alpha_shape_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>


typedef double coord_type;

typedef CGAL::Simple_cartesian<coord_type>  SC;
typedef CGAL::Filtered_kernel<SC> K;
typedef K::Point_2 Point_base;
typedef CGAL::Weighted_point<Point_base,coord_type>  Point;
typedef K::Segment_2  Segment;
typedef K::Line_2  Line;
typedef K::Triangle_2  Triangle;

typedef CGAL::Weighted_alpha_shape_euclidean_traits_2<K> Gt;
typedef CGAL::Regular_triangulation_vertex_base_2<Gt> Rvb;
typedef CGAL::Regular_triangulation_face_base_2<Gt> Rf;

//ExactComparisonTag is Tag_false
typedef CGAL::Alpha_shape_vertex_base_2<Gt,Rvb> Vb;
typedef CGAL::Alpha_shape_face_base_2<Gt, Rf>  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Regular_triangulation_2<Gt,Tds> Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;

//ExactComparisonTag is Tag_true
typedef CGAL::Alpha_shape_vertex_base_2<Gt,Rvb,CGAL::Tag_true> Vb_TT;
typedef CGAL::Alpha_shape_face_base_2<Gt, Rf,CGAL::Tag_true>  Fb_TT;
typedef CGAL::Triangulation_data_structure_2<Vb_TT,Fb_TT> Tds_TT;
typedef CGAL::Regular_triangulation_2<Gt,Tds_TT> Triangulation_2_TT;
typedef CGAL::Alpha_shape_2<Triangulation_2_TT,CGAL::Tag_true>  Alpha_shape_2_TT;

typedef Alpha_shape_2::Face  Face;
typedef Alpha_shape_2::Vertex Vertex;
typedef Alpha_shape_2::Edge Edge;
typedef Alpha_shape_2::Face_handle  Face_handle;
typedef Alpha_shape_2::Vertex_handle Vertex_handle;

typedef Alpha_shape_2::Face_circulator  Face_circulator;
typedef Alpha_shape_2::Vertex_circulator  Vertex_circulator;

typedef Alpha_shape_2::Locate_type Locate_type;

typedef Alpha_shape_2::Face_iterator  Face_iterator;
typedef Alpha_shape_2::Vertex_iterator  Vertex_iterator;
typedef Alpha_shape_2::Edge_iterator  Edge_iterator;
typedef Alpha_shape_2::Edge_circulator  Edge_circulator;

typedef Alpha_shape_2::Alpha_iterator Alpha_iterator;
//---------------------------------------------------------------------

template <class Alpha_shape,class InputIterator, class OutputIterator>
void
alpha_edges(InputIterator begin, InputIterator end,
	    const typename Alpha_shape::FT& Alpha,
	    bool mode,
	    OutputIterator out)
  // Generate Alpha Shape
{
  typedef typename Alpha_shape::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;
  std::vector<Segment> V_seg;
  Alpha_shape A(begin,end);
  
  if (mode) 
    { A.set_mode(Alpha_shape::GENERAL); } 
  else 
    { A.set_mode(Alpha_shape::REGULARIZED); };
  A.set_alpha(Alpha);

  for(Alpha_shape_edges_iterator it =  A.alpha_shape_edges_begin();
      it != A.alpha_shape_edges_end();
      ++it){
    *out++ = A.segment(*it);
  }
}

//---------------------------------------------------------------------
bool
file_input(std::list<Point>& L)
{

  std::ifstream is("./data/fin", std::ios::in);

  if(is.fail())
    {
      std::cerr << "unable to open file for input" << std::endl;
      return false;
    }

  CGAL::set_ascii_mode(is);

  int n;
  is >> n;
  std::cout << "Reading " << n << " points" << std::endl;
  Point_base p(1,2,3);
  for( ; n>0 ; n--)
    {
      is >> p;
      L.push_back(Point (p,coord_type(10)));
    }
  std::cout << "Points inserted" << std::endl;
  return true;
}
    
//------------------ main -------------------------------------------

int main()
{
  std::list<Point> points;
  file_input(points);
  //ExactComparisonTag is Tag_false
  {
    std::vector<Segment> segments;
    alpha_edges<Alpha_shape_2>(points.begin(), points.end(),
          10000.,Alpha_shape_2::GENERAL,
          std::back_inserter(segments));
    std::cout << segments.size() << " alpha shape edges." << std::endl;
  }
  //ExactComparisonTag is Tag_true
  {
    std::vector<Segment> segments;
    alpha_edges<Alpha_shape_2_TT>(points.begin(), points.end(),
          10000.,Alpha_shape_2_TT::GENERAL,
          std::back_inserter(segments));
    std::cout << segments.size() << " alpha shape edges." << std::endl;
  }
  
  return 0;
}
