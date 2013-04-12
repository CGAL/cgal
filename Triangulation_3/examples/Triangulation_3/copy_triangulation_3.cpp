#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/point_generators_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC;
typedef CGAL::Exact_predicates_exact_constructions_kernel EPEC;
typedef CGAL::Delaunay_triangulation_3< EPIC > DT3_epic;
typedef CGAL::Delaunay_triangulation_3< EPEC > DT3_epec;

typedef CGAL::Cartesian_converter<
  CGAL::Exact_predicates_inexact_constructions_kernel,
  CGAL::Exact_predicates_exact_constructions_kernel > Converter;
typedef CGAL::Creator_uniform_3<double,EPIC::Point_3>  Creator;

struct Convert_vertex{
  mutable bool first_vertex;
  Convert_vertex():first_vertex(true) {}
  DT3_epec::Vertex operator()(const DT3_epic::Vertex&) const { return DT3_epec::Vertex(); }
  void operator()(const DT3_epic::Vertex& src,DT3_epec::Vertex& tgt) const
  {
    if (!first_vertex)
      tgt.point() = Converter()( src.point() );
    else
      first_vertex=false;
  }
};

struct Convert_cell{
  DT3_epec::Cell operator()(const DT3_epic::Cell&) const { return DT3_epec::Cell(); }
  void operator()(const DT3_epic::Cell&,DT3_epec::Cell&) const {}  
};

int main()
{
  std::vector< EPIC::Point_3> points;
  CGAL::Random_points_in_sphere_3<EPIC::Point_3,Creator> g(1.0);
  CGAL::cpp11::copy_n( g, 600, std::back_inserter(points) );

  DT3_epic dt3_epic(points.begin(), points.end());
  DT3_epec dt3_epec;
  dt3_epec.set_infinite_vertex( 
    dt3_epec.tds().copy_tds( dt3_epic.tds(),dt3_epic.infinite_vertex(), Convert_vertex(), Convert_cell() ) );

  CGAL_assertion( dt3_epec.is_valid() );
}
