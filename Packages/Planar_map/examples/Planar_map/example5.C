// examples/Planar_map/example5.C
// ------------------------------
#include <CGAL/Homogeneous.h>
#include <CGAL/Pm_segment_exact_traits.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

#include <CGAL/Pm_naive_point_location.h>

typedef CGAL::Homogeneous<long>                  Coord_t;
typedef CGAL::Pm_segment_exact_traits<Coord_t>   Pmtraits;
typedef Pmtraits::Point                          Point;
typedef Pmtraits::X_curve                        Curve;
typedef CGAL::Pm_default_dcel<Pmtraits>          Pmdcel;
typedef CGAL::Planar_map_2<Pmdcel,Pmtraits>      Pmap;

int main()
{
  // creating an instance of Planar_map_2<Pmdcel,Pmtraits>
  //with a naive point location strategy
  CGAL::Pm_naive_point_location<Pmap> naive_pl;
  Pmap pm(&naive_pl);

  Curve cv[4];
  int i;

  CGAL::set_ascii_mode(std::cout);

  Point a1(1, 1), a2(1, 0), a3(0, 0), a4(0, 1);

  cv[0] = Curve(a1, a2);
  cv[1] = Curve(a2, a3);
  cv[2] = Curve(a3, a4);
  cv[3] = Curve(a4, a1);


  // inserting the curves to the map
  Pmap::Halfedge_handle e[4];  

  e[0]=pm.insert_in_face_interior(cv[0],pm.unbounded_face());

  for (i = 1; i < 3; i++)
  {
    e[i]=pm.insert_from_vertex(cv[i],e[i-1]->target(), true);
  }

  e[3]=pm.insert_at_vertices(cv[3],e[0]->source(),e[2]->target() );

  // check the validity of the map
  CGAL_assertion(pm.is_valid());

  return 0;  
}
