
// Test program for the kernel checker.

#include <CGAL/basic.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Kernel_checker.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/Point_3.h>

#include <CGAL/MP_Float.h>
#include <CGAL/Arithmetic_filter.h>
#include <CGAL/Random.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

typedef CGAL::Filtered_exact<double, CGAL::MP_Float> NT1;
// typedef double NT1;

class conv
{
public:
  CGAL::MP_Float operator()(const NT1 &n) const
  {
    return CGAL::to_double(n);
  }
};

typedef CGAL::Cartesian<NT1>            K1;
typedef CGAL::Cartesian<CGAL::MP_Float> K2;
typedef CGAL::Kernel_checker<K1, K2, CGAL::Cartesian_converter<K1, K2, conv> >
        Rep;

typedef Rep::RT  NT;
typedef Rep Gt3d;
typedef CGAL::Triangulation_vertex_base_3<Gt3d> Vb3d;
typedef CGAL::Triangulation_cell_base_3<Gt3d> Ce3d;
typedef CGAL::Triangulation_data_structure_3<Vb3d,Ce3d> Tds3d;
typedef CGAL::Delaunay_triangulation_3<Gt3d, Tds3d> Delaunay3d;

int my_rand()
{
  return int(CGAL::default_random.get_double()*(1<<31));
}

int main()
{
  Delaunay3d D;
  Rep K;
  for (int i=0; i<100; i++)
    D.insert(K.construct_point_3_object()(NT(my_rand()),
                                          NT(my_rand()),
                                          NT(my_rand())));
  return 0;
}
