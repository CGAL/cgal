#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/point_generators_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC;
typedef CGAL::Exact_predicates_exact_constructions_kernel EPEC;
typedef CGAL::Delaunay_triangulation_2< EPIC > DT2_epic;
typedef CGAL::Delaunay_triangulation_2< EPEC > DT2_epec;

typedef CGAL::Cartesian_converter<
  CGAL::Exact_predicates_inexact_constructions_kernel,
  CGAL::Exact_predicates_exact_constructions_kernel > Converter;
typedef CGAL::Creator_uniform_2<double,EPIC::Point_2>  Creator;

struct Convert_vertex{
  mutable bool first_vertex;
  Convert_vertex():first_vertex(true) {}
  DT2_epec::Vertex operator()(const DT2_epic::Vertex&) const { return DT2_epec::Vertex(); }
  void operator()(const DT2_epic::Vertex& src,DT2_epec::Vertex& tgt) const
  {
    if (!first_vertex)
      tgt.point() = Converter()( src.point() );
    else
      first_vertex=false;
  }
};

struct Convert_face{
  DT2_epec::Face operator()(const DT2_epic::Face&) const { return DT2_epec::Face(); }
  void operator()(const DT2_epic::Face&,DT2_epec::Face&) const {}  
};

int main()
{
  std::vector< EPIC::Point_2> points;
  CGAL::Random_points_in_disc_2<EPIC::Point_2,Creator> g(1.0);
  CGAL::cpp11::copy_n( g, 600, std::back_inserter(points) );

  DT2_epic dt2_epic;
  dt2_epic.insert(points.begin(), points.end());
  DT2_epec dt2_epec;
  dt2_epec.set_infinite_vertex( 
    dt2_epec.tds().copy_tds( dt2_epic.tds(),dt2_epic.infinite_vertex(), Convert_vertex(), Convert_face() ) );

  CGAL_assertion( dt2_epec.is_valid() );
}
