#include <iostream>
#include <fstream>
#include <algorithm>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Shape_detection_3.h>
#include <CGAL/structure_point_set.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/IO/read_xyz_points.h>

#include <boost/lexical_cast.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3  Point;
typedef std::pair<Kernel::Point_3, Kernel::Vector_3>         Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

// Efficient RANSAC types
typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits
  <Kernel, Pwn_vector, Point_map, Normal_map>                Traits;
typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits>    Efficient_ransac;
typedef CGAL::Shape_detection_3::Plane<Traits>               Plane;

typedef CGAL::internal::Point_set_structuring<Traits>        Structuring;

typedef CGAL::Advancing_front_surface_reconstruction_vertex_base_3<Kernel> LVb;
typedef CGAL::Advancing_front_surface_reconstruction_cell_base_3<Kernel> LCb;

typedef CGAL::Triangulation_data_structure_3<LVb,LCb> Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel,Tds> Triangulation_3;
typedef typename Triangulation_3::Vertex_handle Vertex_handle;

typedef CGAL::cpp11::array<std::size_t,3> Facet;


struct On_the_fly_pair{
  const Pwn_vector& points;
  typedef std::pair<Point, std::size_t> result_type;

  On_the_fly_pair(const Pwn_vector& points) : points(points) {}
  
  result_type
  operator()(std::size_t i) const
  {
    return result_type(points[i].first,i);
  }
};


template <typename Structuring>
struct Structure_coherence {

  Structuring& structuring;

  Structure_coherence(Structuring& structuring)
    : structuring (structuring)
  {}

  template <typename Vertex_handle>
  bool operator()(const Vertex_handle& p, const Vertex_handle& q, const Vertex_handle& r) const
  {
    Facet f = {{ p->info(), q->info(), r->info() }};
    return (structuring.facet_coherence (f) == 0);
  }
};

typedef CGAL::Advancing_front_surface_reconstruction<Triangulation_3,
                                                     Structure_coherence<Structuring> > Reconstruction;



int main (int argc, char* argv[])
{
  // Points with normals.
  Pwn_vector points;

  // Loading point set from a file. 
  std::ifstream stream(argc>1 ? argv[1] : "data/cube.pwn");

  if (!stream || 
    !CGAL::read_xyz_points_and_normals(stream,
      std::back_inserter(points),
      Point_map(),
      Normal_map()))
  {
      std::cerr << "Error: cannot read file" << std::endl;
      return EXIT_FAILURE;
  }

  std::cerr << "Shape detection... ";
  // Shape detection
  Efficient_ransac ransac;
  ransac.set_input(points);
  ransac.add_shape_factory<Plane>();

  Efficient_ransac::Parameters op;
  op.probability = 0.05;
  op.min_points = 100;
  op.epsilon = (argc>2 ? boost::lexical_cast<double>(argv[2]) : 0.002);
  op.cluster_epsilon = (argc>3 ? boost::lexical_cast<double>(argv[3]) : 0.02);
  op.normal_threshold = 0.7;

  ransac.detect(op);
  std::cerr << "done\nPoint set structuring... ";
  Pwn_vector structured_pts;
  
  Structuring pss (points.begin (), points.end (), ransac);

  pss.run (op.cluster_epsilon);

  pss.get_output (std::back_inserter (structured_pts));
  std::cerr << "done\nAdvancing front... ";
  std::vector<std::size_t> point_indices(boost::counting_iterator<std::size_t>(0),
                                         boost::counting_iterator<std::size_t>(structured_pts.size()));

  Triangulation_3 dt (boost::make_transform_iterator(point_indices.begin(), On_the_fly_pair(structured_pts)),
                      boost::make_transform_iterator(point_indices.end(), On_the_fly_pair(structured_pts)));

  Structure_coherence<Structuring> filter (pss);
  
  Reconstruction R(dt, filter);
  R.run (5., 0.52);
  std::cerr << "done\nWriting result... ";
  std::vector<Facet> output;
  CGAL::write_triple_indices (std::back_inserter (output), R);

  std::ofstream f ("out.off");
  f << "OFF\n" << structured_pts.size () << " " << output.size() << " 0\n";
  for (std::size_t i = 0; i < structured_pts.size (); ++ i)
    f << structured_pts[i].first << std::endl;
  for (std::size_t i = 0; i < output.size (); ++ i)
    f << "3 "
      << output[i][0] << " "
      << output[i][1] << " "
      << output[i][2] << std::endl;
  std::cerr << "all done\n" << std::endl;
  
  return 0;
}
