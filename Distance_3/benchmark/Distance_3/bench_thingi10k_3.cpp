#include <CGAL/Simple_cartesian.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>

#include <CGAL/squared_distance_3.h>
#include <CGAL/Real_timer.h>

#include <CGAL/box_intersection_d.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <boost/container/small_vector.hpp>

#include <cassert>
#include <iostream>

template <typename K>
class Old_compare_squared_distance_3
{
  typedef typename K::FT                 FT;
public:
  typedef typename K::Comparison_result  result_type;

  template <class T1, class T2>
  CGAL::Needs_FT<result_type>
  operator()(const T1& p, const T2& q, const FT& d2) const
  {
    return CGAL::compare(CGAL::internal::squared_distance(p, q, K()), d2);
  }
};

template <typename K>
struct Test
{
  typedef typename K::RT              RT;
  typedef typename K::FT              FT;
  typedef typename K::Comparison_result Comparison_result;
  typedef typename K::Point_3         P;
  typedef typename K::Segment_3       S;
  typedef typename K::Vector_3        V;
  typedef typename K::Ray_3           R;
  typedef typename K::Line_3          L;
  typedef typename K::Triangle_3      T;
  typedef typename K::Plane_3         Pl;
  typedef typename K::Tetrahedron_3   Tet;
  typedef typename K::Iso_cuboid_3    Cub;

  typedef std::vector<boost::container::small_vector<std::size_t, 3> >::iterator Iterator;
  typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator> CBox;

  size_t nb_closed_pairs;

public:
   Test() : nb_closed_pairs(0){ }

  template<bool optimized_version>
  void close_triangles(std::vector<P> &points, std::vector<boost::container::small_vector<std::size_t, 3> >& triangles, FT d2){
    std::vector< CBox > boxes;
    auto extend_bbox3=[&](const Iterator it, FT& d2){
      CGAL::Bbox_3 bb = points[(*it)[0]].bbox()+points[(*it)[1]].bbox()+points[(*it)[2]].bbox();
      return CGAL::Bbox_3(bb.xmin(),bb.ymin(),bb.zmin(),bb.xmax()+CGAL::to_double(d2),bb.ymax()+CGAL::to_double(d2),bb.zmax()+CGAL::to_double(d2));
    };

    auto callback=[&](const CBox &ba, const CBox &bb){
      boost::container::small_vector<std::size_t, 3> &a = *(ba.handle());
      boost::container::small_vector<std::size_t, 3> &b = *(bb.handle());

      std::sort(a.begin(), a.end());
      std::sort(b.begin(), b.end());
      std::vector<int> v;
      std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),std::back_inserter(v));

      if(v.size()!=0) //they have common vertices
        return;

      bool comp;
      T tr1(points[a[0]], points[a[1]], points[a[2]]);
      T tr2(points[b[0]], points[b[1]], points[b[2]]);
      if constexpr(optimized_version)
        comp = K().compare_squared_distance_3_object()(tr1, tr2, d2)!=CGAL::LARGER;
      else
        comp = Old_compare_squared_distance_3<K>()(tr1, tr2, d2)!=CGAL::LARGER;

      if(comp)
      {
        nb_closed_pairs++;
      }
    };

    for(Iterator it=triangles.begin(); it!=triangles.end(); ++it)
      boxes.emplace_back(extend_bbox3(it, d2), it);

    CGAL::box_self_intersection_d(boxes.begin(), boxes.end(), callback);
  }

  void run(std::string filename, FT d2)
  {
    nb_closed_pairs=0;
    std::vector<P> input_points;
    std::vector<boost::container::small_vector<std::size_t, 3>> input_triangles;

    if (!CGAL::IO::read_polygon_soup(filename, input_points, input_triangles))
    {
      std::cerr << "Cannot read " << filename << "\n";
      return;
    }
    CGAL::Real_timer t;
    t.start();
    close_triangles<true>(input_points, input_triangles, d2);
    t.stop();
    std::cout << "New version: #points = " << input_points.size() << " and #triangles = " << input_triangles.size() << " has " << nb_closed_pairs << " pairs at squared distance " << d2 << " in " << t.time() << " sec." << std::endl;
    nb_closed_pairs=0;
    t.reset();
    t.start();
    close_triangles<false>(input_points, input_triangles, d2);
    t.stop();
    std::cout << "Old version: #points = " << input_points.size() << " and #triangles = " << input_triangles.size() << " has " << nb_closed_pairs << " pairs at squared distance " << d2 << " in " << t.time() << " sec." << std::endl;
  }
};

int main(int argc, char** argv)
{
  const std::string filename = argc == 1 ? CGAL::data_file_path("meshes/elephant.off")
                                         : std::string(argv[1]);

  // const std::string out_file = argc <= 2 ? "rounded_soup.off"
  //                                        : std::string(argv[2]);

  std::cout.precision(17);
  std::cerr.precision(17);

  std::cout << "3D Distance bench" << std::endl;

  std::vector<CGAL::Simple_cartesian<double>::Point_3> input_points;
  std::vector<boost::container::small_vector<std::size_t, 3>> input_triangles;

  if (!CGAL::IO::read_polygon_soup(filename, input_points, input_triangles))
  {
    std::cerr << "Cannot read " << filename << "\n";
    return 1;
  }
  CGAL::Bbox_3 bb = CGAL::bbox_3(input_points.begin(), input_points.end());
  double max= (std::max)((std::max)(bb.xmax()-bb.xmin(),bb.ymax()-bb.ymin()),bb.zmax()-bb.zmin());

//  Test<CGAL::Simple_cartesian<double> >().run(filename);
//  Test<CGAL::Simple_homogeneous<double> >().run(filename);
//  Test<CGAL::Simple_cartesian<CGAL::Interval_nt<true> > >(r).run();

  // Test<CGAL::Homogeneous<CGAL::Exact_integer> >(r).run();

  // Test<CGAL::Exact_predicates_inexact_constructions_kernel>().run(filename, max*max);
  // Test<CGAL::Exact_predicates_inexact_constructions_kernel>().run(filename, 10*max*max/input_points.size());

  double average_length=0;
  double min_sq_length=CGAL::squared_distance(input_points[input_triangles[0][0]],input_points[input_triangles[0][1]]);
  for(auto &tr: input_triangles){
    for(int i=0; i<3; ++i){
      double l=CGAL::squared_distance(input_points[tr[i]], input_points[tr[(i+1)%3]]);
      min_sq_length=(std::min)(min_sq_length, l);
      average_length+=std::sqrt(l);
    }
  }
  average_length/=(3*input_triangles.size());

  // Equivalent to EPECK since there are only predicates
  // Test<CGAL::Exact_predicates_inexact_constructions_kernel>().run(filename, 100*max/input_points.size());
  // std::cout << "EPICK" << std::endl;
  // Test<CGAL::Exact_predicates_inexact_constructions_kernel>().run(filename, average_length*average_length/256);
  // Test<CGAL::Exact_predicates_inexact_constructions_kernel>().run(filename, min_sq_length*4);

  std::cout << "EPECK" << std::endl;
  Test<CGAL::Exact_predicates_exact_constructions_kernel>().run(filename, average_length*average_length/256);
  Test<CGAL::Exact_predicates_exact_constructions_kernel>().run(filename, min_sq_length*4);

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
