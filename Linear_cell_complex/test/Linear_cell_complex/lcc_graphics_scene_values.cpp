#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <cassert>
#include <cstdlib>
#include <iostream>

// Face and volume values passed by the LCC drawer to the graphics scene.

using LCC=CGAL::Linear_cell_complex_for_combinatorial_map<3>;
using Point=LCC::Point;
using Dart=LCC::Dart_const_descriptor;
using GSO=CGAL::Graphics_scene_options<LCC, Dart, Dart, Dart, Dart>;

void make_hexa(LCC& lcc, double x)
{
  lcc.make_hexahedron(Point(x,0,0), Point(x+1,0,0), Point(x+1,1,0), Point(x,1,0),
                      Point(x,1,1), Point(x,0,1), Point(x+1,0,1), Point(x+1,1,1));
}

int main()
{
  LCC lcc;
  make_hexa(lcc, 0);
  make_hexa(lcc, 10);

  // Default options: no value at all, so the scene is as before.
  CGAL::Graphics_scene gs;
  CGAL::add_to_graphics_scene(lcc, gs, GSO());
  assert(gs.get_volume_faces().size()==2);
  assert(!gs.has_volume_values() && !gs.has_face_values());
  assert(!gs.is_volume_valued(0) && !gs.is_volume_valued(1));

  // A value on the first hexahedron only.
  GSO gso;
  gso.is_volume_valued=[](const LCC& l, Dart d) { return l.point(d).x()<5; };
  gso.volume_value=[](const LCC&, Dart) { return 0.25f; };
  gso.volume_value_name="quality";
  gso.is_face_valued=[](const LCC&, Dart) { return true; };
  gso.face_value=[](const LCC& l, Dart d) { return float(l.point(d).x()); };

  CGAL::Graphics_scene gs2;
  CGAL::add_to_graphics_scene(lcc, gs2, gso);
  assert(gs2.has_volume_values());
  assert(gs2.get_volume_values().size()==gs2.get_volume_faces().size());
  assert(gs2.is_volume_valued(0)!=gs2.is_volume_valued(1));
  const std::size_t v=gs2.is_volume_valued(0) ? 0 : 1;
  assert(gs2.get_volume_values()[v]==0.25f);
  assert(gs2.volume_value_min()==0.25f && gs2.volume_value_max()==0.25f);
  assert(gs2.volume_value_name()=="quality");
  assert(gs2.has_face_values());
  assert(gs2.get_face_values().size()==gs2.number_of_faces());
  assert(gs2.face_value_min()<gs2.face_value_max());

  // Outside a volume, a volume value is ignored.
  CGAL::Graphics_scene gs3;
  gs3.set_volume_value(1.f);
  assert(!gs3.has_volume_values());

  // clear() drops the values with the faces and volumes.
  gs2.clear();
  assert(!gs2.has_volume_values() && !gs2.has_face_values());
  assert(gs2.get_volume_values().empty() && gs2.get_face_values().empty());

  std::cout<<"lcc_graphics_scene_values: OK"<<std::endl;
  return EXIT_SUCCESS;
}
