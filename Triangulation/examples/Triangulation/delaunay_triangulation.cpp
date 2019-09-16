#include <CGAL/config.h>
#if defined(BOOST_GCC) && (__GNUC__ <= 4) && (__GNUC_MINOR__ < 4)

#include <iostream>
int main()
{
  std::cerr << "NOTICE: This test requires G++ >= 4.4, and will not be compiled." << std::endl;
}

#else

#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>

int main()
{
  double pointsIn[][7] = {
    { 42.89, 0, 60.55, 30.72, 0, 0, 0 },
    { 45.65, 50.83, 50.37, 16.13, 0, 0, 0 },
    { 79.06, 57.84, 61.59, 2.52, 0, 0, 0 },
    { 44.47, 39.46, 39.53, 28.72, 0, 0, 0 },
    { 0, 100, 0, 0, 100, 0, 53.47 },
    { 66.95, 100, 33.6, 0, 0, 0, 0 },
    { 42.89, 0, 0, 30.72, 100, 0, 53.47 },
    { 100, 100, 100, 100, 100, 100, 100 }
  };

  typedef CGAL::Delaunay_triangulation<CGAL::Epick_d< CGAL::Dimension_tag<7> > >      T;
  T dt(7);

  std::vector<T::Point> points;
  points.reserve(8);
  for (int j = 0; j < 8; ++j) {
    T::Point p(&pointsIn[j][0], &pointsIn[j][7]);
    points.push_back(p);
  }

  T::Vertex_handle hint;
  int i = 0;
  for (std::vector<T::Point>::iterator it = points.begin(); it != points.end(); ++it) {
    if (T::Vertex_handle() != hint) {
      hint = dt.insert(*it, hint);
  }
    else {
      hint = dt.insert(*it);
    }
    printf("Processing: %d/%d\n", ++i, (int)points.size());
  }
  return 0;
}

#endif
