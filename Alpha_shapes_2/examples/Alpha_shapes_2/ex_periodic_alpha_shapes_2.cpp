#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>

#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>

#include <fstream>
#include <iostream>

// Traits
typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K>     Gt;

// Vertex type
typedef CGAL::Periodic_2_triangulation_vertex_base_2<Gt>        Vb;
typedef CGAL::Alpha_shape_vertex_base_2<Gt, Vb>                 AsVb;
// Cell type
typedef CGAL::Periodic_2_triangulation_face_base_2<Gt>          Cb;
typedef CGAL::Alpha_shape_face_base_2<Gt, Cb>                   AsCb;

typedef CGAL::Triangulation_data_structure_2<AsVb, AsCb>        Tds;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<Gt, Tds>      P2DT2;
typedef CGAL::Alpha_shape_2<P2DT2>                              Alpha_shape_2;

typedef Gt::Point_2                                             Point;
typedef Gt::Segment_2                                           Segment;
typedef Alpha_shape_2::Alpha_shape_edges_iterator               Alpha_shape_edges_iterator;

template <class OutputIterator>
void alpha_edges( const Alpha_shape_2& A, OutputIterator out)
{
  Alpha_shape_edges_iterator it = A.alpha_shape_edges_begin(),
                             end = A.alpha_shape_edges_end();
  for( ; it!=end; ++it)
    *out++ = A.segment(*it);
}

template <class OutputIterator>
bool file_input(OutputIterator out)
{
  std::ifstream is("./data/fin", std::ios::in);

  if(is.fail())
  {
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

  // Define the periodic square
  P2DT2 pdt(Gt::Iso_rectangle_2(-10,-10, 700,700));

  // Heuristic for inserting large point sets (if pts is reasonably large)
  pdt.insert(points.begin(), points.end(), true);

  // As pdt won't be modified anymore switch to 1-sheeted cover if possible
  if(pdt.is_triangulation_in_1_sheet())
    pdt.convert_to_1_sheeted_covering();
  std::cout << "Periodic Delaunay computed." << std::endl;

  // compute alpha shape
  Alpha_shape_2 as(pdt);
  std::cout << "Alpha shape computed in REGULARIZED mode by default." << std::endl;

   // find optimal alpha values
  Alpha_shape_2::NT alpha_solid = as.find_alpha_solid();
  Alpha_shape_2::Alpha_iterator opt = as.find_optimal_alpha(1);
  std::cout << "Smallest alpha value to get a solid through data points is " << alpha_solid << std::endl;
  std::cout << "Optimal alpha value to get one connected component is " <<  *opt    << std::endl;

  as.set_alpha(*opt);
  assert(as.number_of_solid_components() == 1);

  as.set_alpha(10000);
  std::vector<Segment> segments;
  alpha_edges(as, std::back_inserter(segments));
  std::cout << segments.size() << " alpha shape edges" << std::endl;

  return 0;
}
