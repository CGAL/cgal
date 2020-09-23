#include <iostream>
#include <fstream>
#include <algorithm>
#include <array>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_detection/Efficient_RANSAC.h>
#include <CGAL/structure_point_set.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/disable_warnings.h>

#include <boost/lexical_cast.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3  Point;
typedef std::pair<Kernel::Point_3, Kernel::Vector_3>         Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;

// Efficient RANSAC types
typedef CGAL::Shape_detection::Efficient_RANSAC_traits
  <Kernel, Pwn_vector, Point_map, Normal_map>              Traits;
typedef CGAL::Shape_detection::Efficient_RANSAC<Traits>    Efficient_ransac;
typedef CGAL::Shape_detection::Plane<Traits>               Plane;

// Point set structuring type
typedef CGAL::Point_set_with_structure<Kernel>               Structure;

// Advancing front types
typedef CGAL::Advancing_front_surface_reconstruction_vertex_base_3<Kernel> LVb;
typedef CGAL::Advancing_front_surface_reconstruction_cell_base_3<Kernel> LCb;
typedef CGAL::Triangulation_data_structure_3<LVb,LCb> Tds;
typedef CGAL::Delaunay_triangulation_3<Kernel,Tds> Triangulation_3;
typedef Triangulation_3::Vertex_handle Vertex_handle;

typedef std::array<std::size_t,3> Facet;


// Functor to init the advancing front algorithm with indexed points
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

// Specialized priority functor that favor structure coherence
template <typename Structure>
struct Priority_with_structure_coherence {

  Structure& structure;
  double bound;

  Priority_with_structure_coherence(Structure& structure,
                                    double bound)
    : structure (structure), bound (bound)
  {}

  template <typename AdvancingFront, typename Cell_handle>
  double operator() (AdvancingFront& adv, Cell_handle& c,
                     const int& index) const
  {
    // If perimeter > bound, return infinity so that facet is not used
    if (bound != 0)
      {
        double d  = 0;
        d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                                  c->vertex((index+2)%4)->point()));
        if(d>bound) return adv.infinity();
        d += sqrt(squared_distance(c->vertex((index+2)%4)->point(),
                                   c->vertex((index+3)%4)->point()));
        if(d>bound) return adv.infinity();
        d += sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                                   c->vertex((index+3)%4)->point()));
        if(d>bound) return adv.infinity();
      }

    Facet f = {{ c->vertex ((index + 1) % 4)->info (),
                 c->vertex ((index + 2) % 4)->info (),
                 c->vertex ((index + 3) % 4)->info () }};

    // facet_coherence takes values between -1 and 3, 3 being the most
    // coherent and -1 being incoherent. Smaller weight means higher
    // priority.
    double weight = 100. * (5 - structure.facet_coherence (f));

    return weight * adv.smallest_radius_delaunay_sphere (c, index);
  }

};

// Advancing front type
typedef CGAL::Advancing_front_surface_reconstruction
         <Triangulation_3,
          Priority_with_structure_coherence<Structure> >
        Reconstruction;


int main (int argc, char* argv[])
{
  // Points with normals.
  Pwn_vector points;

  const char* fname = (argc>1) ? argv[1] : "data/cube.pwn";
  // Loading point set from a file.
  std::ifstream stream(fname);

  if (!stream ||
    !CGAL::read_xyz_points(stream,
      std::back_inserter(points),
      CGAL::parameters::point_map(Point_map()).
      normal_map(Normal_map())))
  {
      std::cerr << "Error: cannot read file" << std::endl;
      return EXIT_FAILURE;
  }

  std::cerr << "Shape detection... ";

  Efficient_ransac ransac;
  ransac.set_input(points);
  ransac.add_shape_factory<Plane>(); // Only planes are useful for stucturing

  // Default RANSAC parameters
  Efficient_ransac::Parameters op;
  op.probability = 0.05;
  op.min_points = 100;
  op.epsilon = (argc>2 ? boost::lexical_cast<double>(argv[2]) : 0.002);
  op.cluster_epsilon = (argc>3 ? boost::lexical_cast<double>(argv[3]) : 0.02);
  op.normal_threshold = 0.7;

  ransac.detect(op); // Plane detection

  Efficient_ransac::Plane_range planes = ransac.planes();

  std::cerr << "done\nPoint set structuring... ";

  Pwn_vector structured_pts;
  Structure pss (points,
                 planes,
                 op.cluster_epsilon,  // Same parameter as RANSAC
                 CGAL::parameters::point_map (Point_map()).
                 normal_map (Normal_map()).
                 plane_map (CGAL::Shape_detection::Plane_map<Traits>()).
                 plane_index_map(CGAL::Shape_detection::Point_to_shape_index_map<Traits>(points, planes)));


  for (std::size_t i = 0; i < pss.size(); ++ i)
    structured_pts.push_back (pss[i]);

  std::cerr << "done\nAdvancing front... ";

  std::vector<std::size_t> point_indices(boost::counting_iterator<std::size_t>(0),
                                         boost::counting_iterator<std::size_t>(structured_pts.size()));

  Triangulation_3 dt (boost::make_transform_iterator(point_indices.begin(), On_the_fly_pair(structured_pts)),
                      boost::make_transform_iterator(point_indices.end(), On_the_fly_pair(structured_pts)));


  Priority_with_structure_coherence<Structure> priority (pss,
                                                         1000. * op.cluster_epsilon); // Avoid too large facets
  Reconstruction R(dt, priority);
  R.run ();

  std::cerr << "done\nWriting result... ";

  std::vector<Facet> output;
  const Reconstruction::TDS_2& tds = R.triangulation_data_structure_2();

  for(Reconstruction::TDS_2::Face_iterator fit = tds.faces_begin(); fit != tds.faces_end(); ++fit)
    if(fit->is_on_surface())
      output.push_back (CGAL::make_array(fit->vertex(0)->vertex_3()->id(),
                                         fit->vertex(1)->vertex_3()->id(),
                                         fit->vertex(2)->vertex_3()->id()));

  std::ofstream f ("out.off");
  f << "OFF\n" << structured_pts.size () << " " << output.size() << " 0\n"; // Header
  for (std::size_t i = 0; i < structured_pts.size (); ++ i)
    f << structured_pts[i].first << std::endl;
  for (std::size_t i = 0; i < output.size (); ++ i)
    f << "3 "
      << output[i][0] << " "
      << output[i][1] << " "
      << output[i][2] << std::endl;
  std::cerr << "all done\n" << std::endl;

  f.close();

  return 0;
}
