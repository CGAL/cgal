#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Union_of_balls_3.h>

int main(int, char *[])
{
  double pts[][3]={{3.3874, 3.3577, 2.86547},
                   {4.20832, 3.04325, 3.05838},
                   {3.63033, 2.62921, 2.50657},
                   {4.3492, 2.80494, 1.99437},
                   {5.24092, 2.81322, 2.11588},
                   {6.00076, 3.29489, 2.1409},
                   {5.53583, 3.6421, 1.45294},
                   {5.97231, 2.95352, 1.07171},
                   {5.29922, 3.54395, 0.980338},
                   {5.46575, 3.92853, 0.183865}};
  typedef CGAL::Exact_predicates_inexact_constructions_kernel IKernel;
  typedef IKernel::Point_3                                     Bare_point;
  typedef IKernel::Weighted_point_3        Weighted_point;
  size_t size=sizeof(pts)/(3*sizeof(double));
  std::vector<Weighted_point> l(size);
  for (size_t i=0; i< size; ++i) {
    l[i]= Weighted_point(Bare_point(pts[i][0], pts[i][1], pts[i][2]),
        .9*.9);
    std::cout << ".color " << i << std::endl;
    std::cout << ".sphere " << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << " "
              << .9 << std::endl;
  }

  CGAL::Polyhedron_3<IKernel> p;
  CGAL::Union_of_balls_3<CGAL::Skin_surface_traits_3<IKernel> >
      skin_surface(l.begin(), l.end());

  return 0;
}
