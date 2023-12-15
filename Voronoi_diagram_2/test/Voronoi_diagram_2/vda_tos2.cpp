#include <CGAL/Voronoi_diagram_2.h>
#include "vda_test.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_on_sphere_adaptation_traits_2.h>
#include <CGAL/Projection_on_sphere_traits_3.h>
#include <CGAL/Delaunay_triangulation_on_sphere_adaptation_policies_2.h>
#include <CGAL/Delaunay_triangulation_on_sphere_2.h>

#include <iostream>
#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point_3;
typedef K::Segment_3 Segment_3;

typedef CGAL::Projection_on_sphere_traits_3<K> Traits;
typedef CGAL::Triangulation_on_sphere_face_base_2<Traits> Fb;
typedef CGAL::Triangulation_on_sphere_vertex_base_2<Traits> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Delaunay_triangulation_on_sphere_2<Traits, Tds> DToS2;

typedef Traits::SK SK;
typedef Traits::Point_on_sphere_2 Point_on_sphere_2;
typedef Traits::Arc_on_sphere_2 Arc_on_sphere_2;
typedef DToS2::Face_handle Face_handle;

typedef CGAL::Delaunay_triangulation_on_sphere_adaptation_traits_2<DToS2> AT;
typedef CGAL::Delaunay_triangulation_on_sphere_caching_degeneracy_removal_policy_2<DToS2> AP;
typedef CGAL::Voronoi_diagram_2<DToS2, AT, AP> VD;

int main(int argc, char** argv)
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/radar.xyz");
  const double radius = (argc > 2) ? std::stod(argv[2]) : 1;
  const double step = (argc > 3) ? std::stod(argv[3]) : 0.01;

  std::vector<Point_3> points;
  double x, y, z;

  std::ifstream in(filename);
  if(!in)
  {
    std::cerr << "Invalid input file: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  while(in >> x >> y >> z)
    points.emplace_back(x, y, z);

  std::cout << points.size() << " points" << std::endl;

  Traits traits(Point_3(1,1,1), radius); // radius is 1 by default
  DToS2 dtos(traits);

  Traits::Construct_point_on_sphere_2 cst = traits.construct_point_on_sphere_2_object();

  for(const auto& pt : points)
  {
//    std::cout << "----- Inserting (" << pt
//              << ") at squared distance " << CGAL::squared_distance(pt, traits.center())
//              << " from the center of the sphere" << std::endl;

    // shenanigans just to get a OFF with pts on the sphere (it should just be 'cst(pt)')
    dtos.insert(cst(pt).get_projection(traits.center(), traits.radius()));

//    std::cout << "The triangulation now has dimension: " << dtos.dimension() << " and\n";
//    std::cout << dtos.number_of_vertices() << " vertices" << std::endl;
//    std::cout << dtos.number_of_edges() << " edges" << std::endl;
//    std::cout << dtos.number_of_faces() << " solid faces" << std::endl;
//    std::cout << dtos.number_of_ghost_faces() << " ghost faces" << std::endl;
  }

#ifdef CGAL_VD2_TEST_TOS2_OUTPUT
  CGAL::write_OFF("result.off", dtos, CGAL::parameters::stream_precision(17));

  std::ofstream out_primal("edges_primal.polylines.cgal");
  out_primal.precision(17);
  std::ofstream out_dual("edges_dual.polylines.cgal");
  out_dual.precision(17);
#endif

  for(typename DToS2::All_faces_iterator fit = dtos.all_faces_begin(); fit!=dtos.all_faces_end(); ++fit)
  {
    Face_handle fh = fit;

    const Point_3 c = dtos.dual(fh);
    CGAL_USE(c);

    const Point_on_sphere_2 cs = dtos.dual_on_sphere(fh);
    const FT r = dtos.geom_traits().radius();
    assert(CGAL::abs(CGAL::squared_distance(dtos.construct_point(cs),
                                            dtos.geom_traits().center()) - r*r) <= 1e-10);

    for(int i=0; i<3; ++i)
    {
      Face_handle nfh = fh->neighbor(i);
      if(/*!dtos.is_ghost(nfh) &&*/ fh > nfh)
        continue;

      typename DToS2::Edge e(fh, i);
      const bool diametral_edge = CGAL::collinear(dtos.construct_point(dtos.point(fh, (i+1)%3)),
                                                  dtos.construct_point(dtos.point(fh, (i+2)%3)),
                                                  dtos.geom_traits().center());

      // primal
      if(!diametral_edge)
      {
        const Arc_on_sphere_2 as = dtos.segment_on_sphere(e);
        std::vector<typename SK::Point_3> discretization_points;
        CGAL::Triangulations_on_sphere_2::internal::subsample_arc_on_sphere_2<SK>(as, std::back_inserter(discretization_points), step);
        assert(discretization_points.size() >= 2);

#ifdef CGAL_VD2_TEST_TOS2_OUTPUT
        for(std::size_t i=0; i<discretization_points.size()-1; ++i)
          out_primal << "2 " << discretization_points[i] << " " << discretization_points[i+1] << "\n";
#endif
      }
      else
      {
#ifdef CGAL_VD2_TEST_TOS2_OUTPUT
        const Segment_3 s = dtos.segment(e);
        out_primal << "2 " << s.source() << " " << s.target() << "\n";
#endif
      }

      // Dual
//      if(dtos.is_contour(e))
//        continue;

      const Point_on_sphere_2 c1 = dtos.circumcenter_on_sphere(e.first);
      const Point_on_sphere_2 c2 = dtos.circumcenter_on_sphere(dtos.mirror_edge(e).first);

      // That should never be possible, but with constructions...
      const bool diametral_dual = CGAL::collinear(dtos.construct_point(c1),
                                                  dtos.construct_point(c2),
                                                  dtos.geom_traits().center());
      if(!diametral_dual)
      {
        const Arc_on_sphere_2 ad = dtos.dual_on_sphere(e);
        std::vector<typename SK::Point_3> discretization_points;
        CGAL::Triangulations_on_sphere_2::internal::subsample_arc_on_sphere_2<SK>(ad, std::back_inserter(discretization_points), step);
        assert(discretization_points.size() >= 2);

#ifdef CGAL_VD2_TEST_TOS2_OUTPUT
        for(std::size_t i=0; i<discretization_points.size()-1; ++i)
          out_dual << "2 " << discretization_points[i] << " " << discretization_points[i+1] << "\n";
#endif
      }
      else
      {
#ifdef CGAL_VD2_TEST_TOS2_OUTPUT
        const Segment_3 d = dtos.dual(e);
        out_dual << "2 " << d.source() << " " << d.target() << "\n";
#endif
      };
    }
  }

  std::cout << "Voronoi diagram:" << std::endl;
  AT adaptation_traits(dtos.geom_traits());

  VD vd(dtos, false /*don't give the diagram ownership of the dt*/, adaptation_traits);

  std::cout << vd.number_of_vertices() << " vertices" << std::endl;
  std::cout << vd.number_of_halfedges() << " halfedges" << std::endl;
  std::cout << vd.number_of_faces() << " faces" << std::endl;
  std::cout << "dimension = " << vd.dual().dimension() << std::endl;

  // redirect std::cout to cout_output
  std::stringstream cout_output;
  std::streambuf* old_cout_buf = std::cout.rdbuf(cout_output.rdbuf());

  VD::Vertex_iterator vit = vd.vertices_begin(), vend = vd.vertices_end();
  for(; vit!=vend; ++vit)
  {
    std::cout << "Voronoi vertex at " << vit->point() << std::endl;
    std::cout << "valence: " << vit->degree() << std::endl;
    DToS2::Face_handle dual_face = vit->dual();
    CGAL_USE(dual_face);

    VD::Halfedge_around_vertex_circulator cir = vd.incident_halfedges(vit), cend = cir;
    do {
      std::cout << "  Incident Halfedge with source: " << cir->source()->point() << std::endl;
      VD::Face_handle incident_Voronoi_face = cir->face();
      DToS2::Vertex_handle Voronoi_face_seed = incident_Voronoi_face->dual();
      CGAL_USE(Voronoi_face_seed);
    } while(cir != cend);
  }

  VD::Halfedge_iterator hit = vd.halfedges_begin(), hend = vd.halfedges_end();
  for(; hit!=hend; ++hit)
  {
    std::cout << "Halfedge between " << hit->source()->point() << " and " << hit->target()->point() << std::endl;
    VD::Halfedge_handle next_h = hit->next();
    VD::Halfedge_handle opposite_h = hit->opposite();
    CGAL_USE(next_h);
    CGAL_USE(opposite_h);
  }

  // now restore std::cout and display the output
  std::cout.rdbuf(old_cout_buf);
  const auto output = cout_output.str();
  const auto str_size = output.size();
  const auto str_begin = output.data();
  const auto str_end = str_begin + str_size;
  constexpr auto nb = static_cast<std::remove_cv_t<decltype(str_size)>>(10000);
  auto pos1 = str_begin + (std::min)(nb, str_size);
  assert(pos1 <= str_end);
  const auto pos2 = str_end - (std::min)(nb, str_size);
  assert(pos2 >= str_begin);
  if (pos2 <= pos1) {
    pos1 = str_end;
  }
  std::cout << "NOW THE FIRST AND LAST 10k CHARACTERS OF THE COUT OUTPUT:\n";
  std::cout << "-----\n" << std::string(str_begin, pos1) << "\n-----\n";
  if (pos1 != str_end) {
    std::cout << "[...]\n-----\n" << std::string(pos2, str_end) << "\n-----\n";
  }
  const auto file_name = "vda_tos2_test_output.txt";
  std::ofstream file_output(file_name);
  file_output << output;
  std::cout << "Full log is output to " << file_name << "\n";

  return EXIT_SUCCESS;
}
