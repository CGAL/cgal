#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Geographical_coordinates_traits_2.h>
#include <CGAL/Delaunay_triangulation_on_sphere_traits_2.h>
#include <CGAL/Projection_on_sphere_traits_3.h>
#include <CGAL/Delaunay_triangulation_on_sphere_2.h>

#include <CGAL/enum.h>
#include <CGAL/point_generators_3.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          EPICK;

template <typename K, typename PointRange>
void test(const PointRange& points)
{
//  typedef CGAL::Geographical_coordinates_traits_2<K>                  Gt;
  typedef CGAL::Projection_on_sphere_traits_3<K>                      Gt;
  typedef CGAL::Delaunay_triangulation_on_sphere_2<Gt>                Dt;

  typedef typename Gt::SK                                             SK;

  typedef typename Dt::FT                                             FT;
  typedef typename Dt::Point_3                                        Point_3;
  typedef typename Dt::Point                                          Point_on_sphere_2;
  typedef typename Dt::Segment_3                                      Segment_3;
  typedef typename Dt::Arc_on_sphere_2                                Arc_on_sphere_2;
  typedef typename Dt::Face_handle                                    Face_handle;

  Dt dt;
  dt.insert(points.begin(), points.end());

  // The triangulation, using straight edges
  CGAL::IO::write_OFF("result.off", dt, CGAL::parameters::stream_precision(17));

  std::ofstream out_primal("edges_primal.polylines.cgal");
  out_primal.precision(17);
  std::ofstream out_dual("edges_dual.polylines.cgal");
  out_dual.precision(17);

  for(typename Dt::Solid_faces_iterator fit = dt.solid_faces_begin(); fit!=dt.solid_faces_end(); ++fit)
  {
    Face_handle fh = fit;

    const Point_3 c = dt.dual(fh);
    CGAL_USE(c);

    const Point_on_sphere_2 cs = dt.dual_on_sphere(fh);
    const FT r = dt.geom_traits().radius();
    assert(CGAL::abs(CGAL::squared_distance(dt.construct_point(cs),
                                            dt.geom_traits().center()) - r*r) <= 1e-10);

    for(int i=0; i<3; ++i)
    {
      Face_handle nfh = fh->neighbor(i);
      if(!dt.is_ghost(nfh) && fh < nfh)
        continue;

      typename Dt::Edge e(fh, i);
      const bool diametral_edge = CGAL::collinear(dt.construct_point(dt.point(fh, (i+1)%3)),
                                                  dt.construct_point(dt.point(fh, (i+2)%3)),
                                                  dt.geom_traits().center());

      // primal
      if(!diametral_edge)
      {
        Arc_on_sphere_2 as = dt.segment_on_sphere(e);
        std::vector<typename SK::Point_3> discretization_points;
        CGAL::Triangulations_on_sphere_2::internal::subsample_arc_on_sphere_2<SK>(as, std::back_inserter(discretization_points));
        assert(discretization_points.size() >= 2);

        for(std::size_t i=0; i<discretization_points.size()-1; ++i)
          out_primal << "2 " << discretization_points[i] << " " << discretization_points[i+1] << "\n";
      }
      else
      {
        Segment_3 s = dt.segment(e);
        out_primal << "2 " << s.source() << " " << s.target() << "\n";
      }

      // Dual
      if(dt.is_contour(e))
        continue;

      const Point_on_sphere_2 c1 = dt.circumcenter_on_sphere(e.first);
      const Point_on_sphere_2 c2 = dt.circumcenter_on_sphere(dt.mirror_edge(e).first);

      // That should never be possible, but with constructions...
      const bool diametral_dual = CGAL::collinear(dt.construct_point(c1),
                                                  dt.construct_point(c2),
                                                  dt.geom_traits().center());

      if(!diametral_dual)
      {
        Arc_on_sphere_2 ad = dt.dual_on_sphere(e);
        std::vector<typename SK::Point_3> discretization_points;
        CGAL::Triangulations_on_sphere_2::internal::subsample_arc_on_sphere_2<SK>(ad, std::back_inserter(discretization_points));
        assert(discretization_points.size() >= 2);

        for(std::size_t i=0; i<discretization_points.size()-1; ++i)
          out_dual << "2 " << discretization_points[i] << " " << discretization_points[i+1] << "\n";
      }
      else
      {
        Segment_3 d = dt.dual(e);
        out_dual << "2 " << d.source() << " " << d.target() << "\n";
      };
    }
  }
}

template <typename K>
void test_with_ghost_faces()
{
  std::vector<typename K::Point_3> points;
  points.emplace_back(   1,    0,    0);
  points.emplace_back(   0,    1,    0);
  points.emplace_back(   1,    1,    1);
  points.emplace_back(  -1,   -1,   -1);
  points.emplace_back( 0.5,  0.5,    0);
  points.emplace_back(   1,    1,    0);
  points.emplace_back(0.25, 0.25,    0);
//  points.emplace_back(   0,    0,    1);
//  points.emplace_back(   2,    0,    0);
//  points.emplace_back( 0.5,    0,    0);
//  points.emplace_back( 0.5,    0,  0.5);
//  points.emplace_back( 0.5,    0, -0.5);
//  points.emplace_back(   0,    0,  0.5);
//  points.emplace_back(  -1,    0,    0);

  return test<K>(points);
}

template <typename K>
void test_random_data(const std::size_t num_of_pts)
{
  typedef typename K::Point_3                                         Point_3;
  typedef CGAL::Creator_uniform_3<double, Point_3>                    Creator;

  CGAL::Random r;
  std::cout << "Seed is " << r.get_seed() << std::endl;
  CGAL::Random_points_on_sphere_3<Point_3, Creator> on_sphere(1 /*radius*/, r);

  std::vector<Point_3> points;
  points.reserve(num_of_pts);

  for(std::size_t c=0; c<num_of_pts; ++c)
    points.push_back(*on_sphere++);

  return test<K>(points);
}

int main(int argc, char** argv)
{
  const std::size_t num_pts = (argc > 1) ? std::atoi(argv[1]) : 100;

  test_with_ghost_faces<EPICK>();
  test_random_data<EPICK>(num_pts);

  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
