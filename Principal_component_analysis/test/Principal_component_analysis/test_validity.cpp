/*
  This test checks that the outputs of all variants of PCA are
  correct.

  To do so:
  - N Point_3 objects are generated on a random line
  - N*N Point_3 objects are generated on a random plane

  For each variant of the PCA, we generate objects "centered" on the
  reference points:
  - Sphere_3 centered on each point with random radius
  - Equilateral Triangle_3 centered on each point with random length
  - Point_2 using first 2 coordinates of Point_3
  - etc.

  PCA is then computed for Line_2 in 2D or for Line_3 and Plane_3 in
  3D. Then, the mean distance between the original points (or 2D
  variants if 2D case) and the estimated support is computed: if it
  higher than 0.01% of the diameter of the point cloud, an assertion
  is raised.

  This only tests an obvious case where the PCA is computed on a set
  of exactly aligned objects, but it should help avoiding the
  introduction of bugs.
 */

#include <CGAL/Simple_cartesian.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Random.h>

#include <fstream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Segment_3 Segment_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Sphere_3 Sphere_3;
typedef Kernel::Triangle_3 Triangle_3;
typedef Kernel::Tetrahedron_3 Tetrahedron_3;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Line_3 Line_3;

typedef Kernel::Point_2 Point_2;
typedef Kernel::Vector_2 Vector_2;
typedef Kernel::Circle_2 Circle_2;
typedef Kernel::Iso_rectangle_2 Iso_rectangle_2;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Triangle_2 Triangle_2;
typedef Kernel::Line_2 Line_2;

CGAL::Random rnd;

/*
  Functions used to generate objects "centered" on points (objects
  whose centers of mass are on the points).
 */

void generate_objects_centered_on_points (const std::vector<Point_3>& source,
                                          std::vector<Point_3>& target)
{
  target.reserve (source.size());
  std::copy (source.begin(), source.end(), std::back_inserter (target));
}

void generate_objects_centered_on_points (const std::vector<Point_3>& source,
                                          std::vector<Point_2>& target)
{
  target.reserve (source.size());
  for (std::size_t i = 0; i < source.size(); ++ i)
    target.push_back (Point_2 (source[i].x(), source[i].y()));
}

void generate_objects_centered_on_points (const std::vector<Point_3>& source,
                                          std::vector<Segment_3>& target)
{
  Vector_3 diff (rnd.get_double(0, 0.01),
                 rnd.get_double(0, 0.01),
                 rnd.get_double(0, 0.01));

  target.reserve (source.size());
  for (std::size_t i = 0; i < source.size(); ++ i)
    target.push_back (Segment_3 (source[i] + diff, source[i] - diff));
}

void generate_objects_centered_on_points (const std::vector<Point_3>& source,
                                          std::vector<Segment_2>& target)
{
  Vector_2 diff (rnd.get_double(0, 0.01),
                 rnd.get_double(0, 0.01));

  target.reserve (source.size());
  for (std::size_t i = 0; i < source.size(); ++ i)
  {
    Point_2 p (source[i].x(), source[i].y());
    target.push_back (Segment_2 (p + diff, p - diff));
  }
}

void generate_objects_centered_on_points (const std::vector<Point_3>& source,
                                          std::vector<Sphere_3>& target)
{
  double radius = rnd.get_double(0, 0.01);
  target.reserve (source.size());
  for (std::size_t i = 0; i < source.size(); ++ i)
    target.push_back (Sphere_3 (source[i], radius * radius));
}

void generate_objects_centered_on_points (const std::vector<Point_3>& source,
                                          std::vector<Circle_2>& target)
{
  double radius = rnd.get_double(0, 0.01);
  target.reserve (source.size());
  for (std::size_t i = 0; i < source.size(); ++ i)
  {
    Point_2 p (source[i].x(), source[i].y());
    target.push_back (Circle_2 (p, radius * radius));
  }
}

void generate_objects_centered_on_points (const std::vector<Point_3>& source,
                                          std::vector<Triangle_3>& target)
{
  // General equilateral triangles
  double length = rnd.get_double(0, 0.01);
  double height = length * std::sqrt(3.) / 2.;

  target.reserve (source.size());
  for (std::size_t i = 0; i < source.size(); ++ i)
    target.push_back
      (Triangle_3 (source[i] + Vector_3 (0, 0, 2. * height / 3.),
                   source[i] + Vector_3 (0, -length/2., -height/3.),
                   source[i] + Vector_3 (0, length/2., -height/3.)));
}

void generate_objects_centered_on_points (const std::vector<Point_3>& source,
                                          std::vector<Triangle_2>& target)
{
  // General equilateral triangles
  double length = rnd.get_double(0, 0.01);
  double height = length * std::sqrt(3.) / 2.;

  target.reserve (source.size());
  for (std::size_t i = 0; i < source.size(); ++ i)
  {
    Point_2 p (source[i].x(), source[i].y());
    target.push_back
      (Triangle_2 (p + Vector_2 (0, 2. * height / 3.),
                   p + Vector_2 (-length/2., -height/3.),
                   p + Vector_2 (length/2., -height/3.)));
  }
}

void generate_objects_centered_on_points (const std::vector<Point_3>& source,
                                          std::vector<Tetrahedron_3>& target)
{
  // General regular tetrahedra
  double length = rnd.get_double(0, 0.01);
  double height = length * std::sqrt(2. / 3.);

  double height_tri = length * std::sqrt(3.) / 2.;

  target.reserve (source.size());
  for (std::size_t i = 0; i < source.size(); ++ i)
    target.push_back
      (Tetrahedron_3 (source[i] + Vector_3 (0, 0, 3. * height / 4.),
                      source[i] + Vector_3 (0, 2. * height_tri / 3., -height/4.),
                      source[i] + Vector_3 (-length/2., -height_tri / 3., -height/4.),
                      source[i] + Vector_3 (length/2., -height_tri / 3., -height/4.)));
}

void generate_objects_centered_on_points (const std::vector<Point_3>& source,
                                          std::vector<Iso_cuboid_3>& target)
{
  Vector_3 diff (rnd.get_double(0, 0.01),
                 rnd.get_double(0, 0.01),
                 rnd.get_double(0, 0.01));

  target.reserve (source.size());
  for (std::size_t i = 0; i < source.size(); ++ i)
    target.push_back (Iso_cuboid_3 (source[i] + diff, source[i] - diff));
}

void generate_objects_centered_on_points (const std::vector<Point_3>& source,
                                          std::vector<Iso_rectangle_2>& target)
{
  Vector_2 diff (rnd.get_double(0, 0.01),
                 rnd.get_double(0, 0.01));

  target.reserve (source.size());
  for (std::size_t i = 0; i < source.size(); ++ i)
  {
    Point_2 p (source[i].x(), source[i].y());
    target.push_back (Iso_rectangle_2 (p + diff, p - diff));
  }
}

/*
  Compute distances between original points to Line_3/Plane_3/Line_2
 */

double distance (const Point_3& p, const Line_3& l)
{ return std::sqrt (CGAL::squared_distance (p, l)); }
double distance (const Point_3& p, const Plane_3& pl)
{ return std::sqrt (CGAL::squared_distance (p, pl)); }
double distance (const Point_3& p, const Line_2& l)
{ return std::sqrt (CGAL::squared_distance (Point_2 (p.x(), p.y()), l)); }

/*
  Test quality of fit and raise assertion if too low.
 */

template <typename Fitted>
void assert_quality (const std::vector<Point_3>& points, const Fitted& fitted)
{
  double mean_dist = 0;
  for (std::size_t i = 0; i < points.size(); ++ i)
  {
    double dist = distance (points[i], fitted);
    mean_dist += dist;
  }
  mean_dist /= points.size();

  std::cerr << "mean distance = " << mean_dist << std::endl;

  CGAL_assertion_code
    (double limit = 1e-5 * std::sqrt (CGAL::squared_distance (points.front(), points.back())));
  CGAL_assertion (mean_dist < limit);
}

/*
  Reference test functions, both in 2D and 3D.
 */

template <typename Fitted, typename Object, int dim>
void test (const std::vector<Object>& objects,
           const std::vector<Point_3>& points,
           const CGAL::Dimension_tag<2>& /* ambient dimension */)
{
  std::cerr << " CGAL::Dimension_tag<" << dim << ">: ";
  Fitted fitted;
  CGAL::linear_least_squares_fitting_2
    (objects.begin(), objects.end(), fitted, CGAL::Dimension_tag<dim>());
  assert_quality (points, fitted);
}

template <typename Fitted, typename Object, int dim>
void test (const std::vector<Object>& objects,
           const std::vector<Point_3>& points,
           const CGAL::Dimension_tag<3>& /* ambient dimension */)
{
  std::cerr << " CGAL::Dimension_tag<" << dim << ">: ";
  Fitted fitted;
  CGAL::linear_least_squares_fitting_3
    (objects.begin(), objects.end(), fitted, CGAL::Dimension_tag<dim>());
  assert_quality (points, fitted);
}

template <typename Object, typename Fitted, int dim>
void test (const std::vector<Point_3>& points)
{
  std::vector<Object> objects;
  generate_objects_centered_on_points (points, objects);
  test<Fitted, Object, dim>
    (objects, points, typename CGAL::Ambient_dimension<Object, Kernel>::type());
}




int main()
{
  std::cerr << "Validity test with seed " << rnd.get_seed() << std::endl;

  Point_3 origin (rnd.get_double(), rnd.get_double(), rnd.get_double());
  Vector_3 base1 (rnd.get_double(), rnd.get_double(), rnd.get_double());
  Vector_3 base2 (rnd.get_double(), rnd.get_double(), rnd.get_double());

  std::cerr << "Origin = " << origin << std::endl
            << "Base1 = " << base1 << std::endl
            << "Base2 = " << base2 << std::endl;

  std::size_t nb_points_on_line = 100;
  std::vector<Point_3> points_on_line;
  points_on_line.reserve (nb_points_on_line);
  for (std::size_t i = 0; i < nb_points_on_line; ++ i)
    points_on_line.push_back (Point_3 (origin + double(i) * base1));

  std::vector<Point_3> points_on_plane;
  points_on_plane.reserve (nb_points_on_line * nb_points_on_line);
  for (std::size_t i = 0; i < nb_points_on_line; ++ i)
    for (std::size_t j = 0; j < nb_points_on_line; ++ j)
      points_on_plane.push_back (Point_3 (origin + double(i) * base1 + double(j) * base2));

  std::cerr << std::endl << "=== 2D ===" << std::endl << std::endl;

  std::cerr << "[Testing line fitting on Point_2 objects]" << std::endl;
  test<Point_2, Line_2, 0> (points_on_line);

  std::cerr << "[Testing line fitting on Segment_2 objects]" << std::endl;
  test<Segment_2, Line_2, 1> (points_on_line);
  test<Segment_2, Line_2, 0> (points_on_line);

  std::cerr << "[Testing line fitting on Circle_2 objects]" << std::endl;
  test<Circle_2, Line_2, 2> (points_on_line);
  test<Circle_2, Line_2, 1> (points_on_line);

  std::cerr << "[Testing line fitting on Triangle_2 objects]" << std::endl;
  test<Triangle_2, Line_2, 2> (points_on_line);
  test<Triangle_2, Line_2, 1> (points_on_line);
  test<Triangle_2, Line_2, 0> (points_on_line);

  std::cerr << "[Testing line fitting on Iso_rectangle_2 objects]" << std::endl;
  test<Iso_rectangle_2, Line_2, 2> (points_on_line);
  test<Iso_rectangle_2, Line_2, 1> (points_on_line);
  test<Iso_rectangle_2, Line_2, 0> (points_on_line);

  std::cerr << std::endl << "=== 3D ===" << std::endl << std::endl;

  std::cerr << "[Testing line fitting on Point_3 objects]" << std::endl;
  test<Point_3, Line_3, 0> (points_on_line);

  std::cerr << "[Testing plane fitting on Point_3 objects]" << std::endl;
  test<Point_3, Plane_3, 0> (points_on_plane);

  std::cerr << "[Testing line fitting on Segment_3 objects]" << std::endl;
  test<Segment_3, Line_3, 1> (points_on_line);
  test<Segment_3, Line_3, 0> (points_on_line);

  std::cerr << "[Testing plane fitting on Segment_3 objects]" << std::endl;
  test<Segment_3, Plane_3, 1> (points_on_plane);
  test<Segment_3, Plane_3, 0> (points_on_plane);

  std::cerr << "[Testing line fitting on Sphere_3 objects]" << std::endl;
  test<Sphere_3, Line_3, 3> (points_on_line);
  test<Sphere_3, Line_3, 2> (points_on_line);

  std::cerr << "[Testing plane fitting on Sphere_3 objects]" << std::endl;
  test<Sphere_3, Plane_3, 3> (points_on_plane);
  test<Sphere_3, Plane_3, 2> (points_on_plane);

  std::cerr << "[Testing line fitting on Triangle_3 objects]" << std::endl;
  test<Triangle_3, Line_3, 2> (points_on_line);
  test<Triangle_3, Line_3, 1> (points_on_line);
  test<Triangle_3, Line_3, 0> (points_on_line);

  std::cerr << "[Testing plane fitting on Triangle_3 objects]" << std::endl;
  test<Triangle_3, Plane_3, 2> (points_on_plane);
  test<Triangle_3, Plane_3, 1> (points_on_plane);
  test<Triangle_3, Plane_3, 0> (points_on_plane);

  std::cerr << "[Testing line fitting on Tetrahedron_3 objects]" << std::endl;
  test<Tetrahedron_3, Line_3, 3> (points_on_line);
  test<Tetrahedron_3, Line_3, 2> (points_on_line);
  test<Tetrahedron_3, Line_3, 1> (points_on_line);
  test<Tetrahedron_3, Line_3, 0> (points_on_line);

  std::cerr << "[Testing plane fitting on Tetrahedron_3 objects]" << std::endl;
  test<Tetrahedron_3, Plane_3, 3> (points_on_plane);
  test<Tetrahedron_3, Plane_3, 2> (points_on_plane);
  test<Tetrahedron_3, Plane_3, 1> (points_on_plane);
  test<Tetrahedron_3, Plane_3, 0> (points_on_plane);

  std::cerr << "[Testing line fitting on Iso_cuboid_3 objects]" << std::endl;
  test<Iso_cuboid_3, Line_3, 3> (points_on_line);
  test<Iso_cuboid_3, Line_3, 2> (points_on_line);
  test<Iso_cuboid_3, Line_3, 1> (points_on_line);
  test<Iso_cuboid_3, Line_3, 0> (points_on_line);

  std::cerr << "[Testing plane fitting on Iso_cuboid_3 objects]" << std::endl;
  test<Iso_cuboid_3, Plane_3, 3> (points_on_plane);
  test<Iso_cuboid_3, Plane_3, 2> (points_on_plane);
  test<Iso_cuboid_3, Plane_3, 1> (points_on_plane);
  test<Iso_cuboid_3, Plane_3, 0> (points_on_plane);

  return EXIT_SUCCESS;
}
