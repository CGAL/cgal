#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_xy_3.h>

#include <CGAL/algorithm.h>

#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>

typedef double                                coord_type;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;
typedef CGAL::Projection_traits_xy_3<Epic>      K;

typedef K::Point_2                            Point;
typedef K::Segment_2                          Segment;

// -------------------------------------------------------------------
// Since K::Point_2 is here in fact CGAL::Point_3<FK>, the basic Cartesian_converter
// cannot be used (and thus ExactAlphaComparisonTag cannot be set to 'true') because
// it does not know how to convert from CGAL::Point_3<FK> to CGAL::Point_2<EK>.

// Thus, we must provide a specialization of Cartesian_converter to be able
// to set ExactAlphaComparisonTag to 'true'

namespace CGAL {

template < class K2, class C >
class Cartesian_converter<Epic, K2, C>
{
public:
  typedef CGAL::Projection_traits_xy_3<Epic> Source_kernel;
  typedef K2                                 Target_kernel;
  typedef C                                  Number_type_converter;

  typedef typename Source_kernel::Point_2    SP2;
  typedef typename Target_kernel::Point_2    TP2;

  TP2 operator()(const SP2& p) const
  {
    return TP2(c(p.x()), c(p.y()));
  }

private:
  C c;
};

} // namespace CGAL

// The partial specialization must be defined before Alpha Shapes-related headers
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>

// ExactAlphaComparisonTag is false
typedef K Gt;
typedef CGAL::Alpha_shape_vertex_base_2<Gt>             Vb;
typedef CGAL::Alpha_shape_face_base_2<Gt>               Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>     Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds>          Triangulation_2;
typedef CGAL::Alpha_shape_2<Triangulation_2>            Alpha_shape_2;

// ExactAlphaComparisonTag is true
typedef K Gt;
typedef CGAL::Alpha_shape_vertex_base_2<Gt,CGAL::Default,CGAL::Tag_true> Vb_TT;
typedef CGAL::Alpha_shape_face_base_2<Gt,CGAL::Default,CGAL::Tag_true>   Fb_TT;
typedef CGAL::Triangulation_data_structure_2<Vb_TT,Fb_TT>                Tds_TT;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds_TT>                        Triangulation_2_TT;
typedef CGAL::Alpha_shape_2<Triangulation_2_TT,CGAL::Tag_true>           Alpha_shape_2_TT;

template <class Alpha_shape,class InputIterator, class OutputIterator>
void alpha_edges(InputIterator begin, InputIterator end,
                 const typename Alpha_shape::FT& Alpha,
                 bool mode,
                 OutputIterator out)
{
  typedef typename Alpha_shape::Alpha_shape_edges_iterator Alpha_shape_edges_iterator;
  Alpha_shape A(begin,end);

  if (mode) { A.set_mode(Alpha_shape::GENERAL); }
  else { A.set_mode(Alpha_shape::REGULARIZED); }

  A.set_alpha(Alpha);

  for(Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin();
      it != A.alpha_shape_edges_end(); ++it) {
    *out++ = A.segment(*it);
  }
}

template <class OutputIterator>
bool file_input(OutputIterator out)
{
  std::ifstream is("./data/fin3", std::ios::in);

  if(is.fail()) {
    std::cerr << "unable to open file for input" << std::endl;
    return false;
  }

  int n;
  is >> n;
  std::cout << "Reading " << n << " points from file" << std::endl;
  std::copy_n(std::istream_iterator<Point>(is), n, out);

  return true;
}

int main()
{
  std::list<Point> points;
  if(! file_input(std::back_inserter(points)))
    return -1;

  // ExactAlphaComparisonTag is False
  {
    std::vector<Segment> segments;
    alpha_edges<Alpha_shape_2>(points.begin(), points.end(),
                               10000.,Alpha_shape_2::GENERAL,
                               std::back_inserter(segments));

    std::cout << "Alpha Shape computed with ExactAlphaComparisonTag = false" << std::endl;
    std::cout << segments.size() << " alpha shape edges" << std::endl;
  }

  // ExactAlphaComparisonTag is True
  {
    std::vector<Segment> segments;
    alpha_edges<Alpha_shape_2_TT>(points.begin(), points.end(),
                                  10000.,Alpha_shape_2_TT::GENERAL,
                                  std::back_inserter(segments));

    std::cout << "Alpha Shape computed with ExactAlphaComparisonTag = true" << std::endl;
    std::cout << segments.size() << " alpha shape edges" << std::endl;
  }

  return 0;
}
