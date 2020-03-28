/***********************************************************************

Takes a list of weighted points and returns a list of segments
corresponding to the weighted Alpha Shape.

************************************************************************/

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <list>
#include <vector>

typedef double                                          coord_type;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                      Point;
typedef K::Weighted_point_2                             Weighted_point;
typedef K::Segment_2                                    Segment;
typedef K::Line_2                                       Line;
typedef K::Triangle_2                                   Triangle;

typedef K                                               Gt;
typedef CGAL::Regular_triangulation_vertex_base_2<Gt>   Rvb;
typedef CGAL::Regular_triangulation_face_base_2<Gt>     Rf;

//ExactComparisonTag is Tag_false
typedef CGAL::Alpha_shape_vertex_base_2<Gt,Rvb>         Vb;
typedef CGAL::Alpha_shape_face_base_2<Gt, Rf>           Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>     Tds;
typedef CGAL::Regular_triangulation_2<Gt,Tds>           Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>            Alpha_shape_2;

//ExactComparisonTag is Tag_true
typedef CGAL::Alpha_shape_vertex_base_2<Gt,Rvb,
                                        CGAL::Tag_true /* exact */,
                                        CGAL::Tag_true /* weighted */> Vb_TT;
typedef CGAL::Alpha_shape_face_base_2<Gt, Rf,
                                      CGAL::Tag_true /* exact */,
                                      CGAL::Tag_true /* weighted */>   Fb_TT;
typedef CGAL::Triangulation_data_structure_2<Vb_TT,Fb_TT>              Tds_TT;
typedef CGAL::Regular_triangulation_2<Gt,Tds_TT>                       Triangulation_2_TT;
typedef CGAL::Alpha_shape_2<Triangulation_2_TT,CGAL::Tag_true>         Alpha_shape_2_TT;

typedef Alpha_shape_2::Face               Face;
typedef Alpha_shape_2::Vertex             Vertex;
typedef Alpha_shape_2::Edge               Edge;
typedef Alpha_shape_2::Face_handle        Face_handle;
typedef Alpha_shape_2::Vertex_handle      Vertex_handle;

typedef Alpha_shape_2::Face_circulator    Face_circulator;
typedef Alpha_shape_2::Vertex_circulator  Vertex_circulator;

typedef Alpha_shape_2::Locate_type        Locate_type;

typedef Alpha_shape_2::Face_iterator      Face_iterator;
typedef Alpha_shape_2::Vertex_iterator    Vertex_iterator;
typedef Alpha_shape_2::Edge_iterator      Edge_iterator;
typedef Alpha_shape_2::Edge_circulator    Edge_circulator;

typedef Alpha_shape_2::Alpha_iterator     Alpha_iterator;

template <class Alpha_shape,class InputIterator, class OutputIterator>
void alpha_edges(InputIterator begin, InputIterator end,
                 const typename Alpha_shape::FT& Alpha,
                 bool mode,
                 OutputIterator out)
{
  typedef typename Alpha_shape::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;

  // Generate Alpha Shape
  Alpha_shape A(begin, end);

  if (mode) {
    A.set_mode(Alpha_shape::GENERAL);
  } else {
    A.set_mode(Alpha_shape::REGULARIZED);
  }

  A.set_alpha(Alpha);

  Alpha_shape_edges_iterator eit = A.alpha_shape_edges_begin(),
                             eend = A.alpha_shape_edges_end();
  for( ; eit!=eend; ++eit)
    *out++ = A.segment(*eit);
}

bool file_input(std::list<Weighted_point>& L)
{
  std::ifstream is("./data/fin_weighted", std::ios::in);

  if(is.fail())
  {
    std::cerr << "unable to open file for input" << std::endl;
    return false;
  }

  int n;
  is >> n;
  std::cout << "Reading " << n << " points" << std::endl;
  Weighted_point wp;
  for( ; n>0 ; n--)
  {
    is >> wp;
    L.push_back(wp);
  }

  std::cout << "Points inserted" << std::endl;
  return true;
}

//------------------ main -------------------------------------------
int main()
{
  std::list<Weighted_point> wpoints;
  if(!file_input(wpoints))
    return -1;


  //ExactComparisonTag is Tag_false
  {
    std::vector<Segment> segments;
    alpha_edges<Alpha_shape_2>(wpoints.begin(), wpoints.end(),
          10000.,Alpha_shape_2::GENERAL,
          std::back_inserter(segments));
    std::cout << segments.size() << " alpha shape edges." << std::endl;
  }

  //ExactComparisonTag is Tag_true
  {
    std::vector<Segment> segments;
    alpha_edges<Alpha_shape_2_TT>(wpoints.begin(), wpoints.end(),
          10000.,Alpha_shape_2_TT::GENERAL,
          std::back_inserter(segments));
    std::cout << segments.size() << " alpha shape edges." << std::endl;
  }

  return 0;
}
