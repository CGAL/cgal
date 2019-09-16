#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <fstream>
#include <iostream>
#include <list>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::FT                                               FT;
typedef K::Weighted_point_2                                 Weighted_point;
typedef K::Segment_2                                        Segment;

typedef CGAL::Regular_triangulation_vertex_base_2<K>        Rvb;
typedef CGAL::Alpha_shape_vertex_base_2<K,Rvb>              Vb;
typedef CGAL::Regular_triangulation_face_base_2<K>          Rf;
typedef CGAL::Alpha_shape_face_base_2<K,Rf>                 Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>         Tds;
typedef CGAL::Regular_triangulation_2<K,Tds>                Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>                Alpha_shape_2;

typedef Alpha_shape_2::Alpha_shape_edges_iterator           Alpha_shape_edges_iterator;

template <class OutputIterator>
void alpha_edges(const Alpha_shape_2& A, OutputIterator out)
{
  Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(),
                             end = A.alpha_shape_edges_end();
  for( ; it!=end; ++it)
    *out++ = A.segment(*it);
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
  for( ; n>0; n--)
  {
    Weighted_point wp;
    is >> wp;
    L.push_back(wp);
  }

  return true;
}

// Reads a list of points and returns a list of segments corresponding to
// the weighted Alpha Shape.
int main()
{
  std::list<Weighted_point> wpoints;
  if(!file_input(wpoints))
    return -1;

  Alpha_shape_2 A(wpoints.begin(), wpoints.end(),
                  FT(10000),
                  Alpha_shape_2::GENERAL);

  std::vector<Segment> segments;
  alpha_edges(A, std::back_inserter(segments));

  std::cout << "Alpha Shape computed" << std::endl;
  std::cout << segments.size() << " alpha shape edges" << std::endl;
  std::cout << "Optimal alpha: " << *A.find_optimal_alpha(1)<<std::endl;
  return 0;
}
