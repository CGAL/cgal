
// Test program for the new filtering wrapper for kernel traits.
// Currently, it's only a compilation test.

#include <CGAL/basic.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <CGAL/Random.h>
#include <CGAL/MP_Float.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

typedef CGAL::Filtered_kernel<CGAL::Cartesian<double>,
                              CGAL::Simple_cartesian<CGAL::MP_Float> > Rep;

typedef Rep::FT  NT;
typedef CGAL::Triangulation_geom_traits_3<Rep> Gt3d;
typedef CGAL::Triangulation_vertex_base_3<Gt3d> Vb3d;
// typedef CGAL::Triangulation_vertex_base_pointer_3<Gt3d> Vb3d;
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
  int loops = 100;
  for (int i=0; i<loops; i++)
    D.insert(K.construct_point_3_object()(NT(my_rand()),
                                          NT(my_rand()),
                                          NT(my_rand())));
  return 0;
}
