#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <cstdlib>
#include <iostream>
#include <vector>

// A face shared by two volumes is stored once in the graphics scene. The scene
// records the volume of the stored copy (its normal points out of this volume)
// and whether its color is the volume one, so that the viewer can draw it
// correctly for the other volume when the first one is clipped.

using LCC=CGAL::Linear_cell_complex_for_combinatorial_map<3>;
using Point=LCC::Point;
using Dart=LCC::Dart_const_descriptor;
using GSO=CGAL::Graphics_scene_options<LCC, Dart, Dart, Dart, Dart>;
using GS=CGAL::Graphics_scene;

void check(bool b, const char* msg)
{
  if(!b)
  {
    std::cerr<<"lcc_graphics_scene_shared_faces FAILED: "<<msg<<std::endl;
    std::exit(EXIT_FAILURE);
  }
}

// Whether the normal of face f points to the center of volume v.
bool side(const GS& gs, int v, unsigned int f)
{
  const std::vector<float>& pos=gs.get_array_of_index(GS::POS_FACES);
  const std::vector<float>& normals=gs.get_array_of_index(GS::FLAT_NORMAL_FACES);
  const std::pair<unsigned int, unsigned int>& r=gs.face_range(f);
  const CGAL::Bbox_3& b=gs.get_volume_bboxes()[v];
  const double c[3]={(b.xmin()+b.xmax())/2., (b.ymin()+b.ymax())/2., (b.zmin()+b.zmax())/2.};
  double d=0.;
  for(int k=0; k<3; ++k)
  { d+=normals[3*r.first+k]*(c[k]-pos[3*r.first+k]); }
  return d>0.;
}

void make_hexa(LCC& lcc, double x)
{
  lcc.make_hexahedron(Point(x,0,0), Point(x+1,0,0), Point(x+1,1,0), Point(x,1,0),
                      Point(x,1,1), Point(x,0,1), Point(x+1,0,1), Point(x+1,1,1));
}

int main()
{
  // Two unit cubes sharing the face x=1.
  LCC lcc;
  make_hexa(lcc, 0);
  make_hexa(lcc, 1);
  check(lcc.sew3_same_facets()==1, "the two cubes share one face");

  GSO gso;
  gso.colored_volume=[](const LCC&, Dart) { return true; };
  gso.volume_color=[](const LCC&, Dart) { return CGAL::IO::Color(200,50,50); };
  GS gs;
  CGAL::add_to_graphics_scene(lcc, gs, gso);

  const std::vector<std::vector<unsigned int>>& vols=gs.get_volume_faces();
  check(vols.size()==2, "two volumes");
  check(gs.number_of_faces()==11, "the shared face is stored once");

  // The shared face is in both volumes.
  std::vector<unsigned int> shared;
  for(unsigned int f: vols[0])
    for(unsigned int g: vols[1])
      if(f==g) { shared.push_back(f); }
  check(shared.size()==1, "one face in both volumes");
  const unsigned int fs=shared[0];

  // Every other face belongs to its own volume.
  for(int v=0; v<2; ++v)
    for(unsigned int f: vols[v])
      if(f!=fs) { check(gs.face_volume(f)==v, "face of its own volume"); }

  // The shared face is the copy of one of the two volumes: it is oriented like the
  // other faces of this volume, and the opposite way of those of the other volume.
  const int owner=gs.face_volume(fs);
  check(owner==0 || owner==1, "shared face belongs to one of the volumes");
  for(int v=0; v<2; ++v)
    for(unsigned int f: vols[v])
      if(f!=fs)
      {
        if(v==owner)
        { check(side(gs, v, f)==side(gs, v, fs), "oriented like its volume"); }
        else
        { check(side(gs, v, f)!=side(gs, v, fs), "opposite to the other volume"); }
      }
  check(gs.face_has_volume_color(fs), "shared face has the volume color");

  // A face with its own color does not have the volume color.
  GSO gso2=gso;
  gso2.colored_face=[](const LCC&, Dart) { return true; };
  gso2.face_color=[](const LCC&, Dart) { return CGAL::IO::Color(50,50,200); };
  GS gs2;
  CGAL::add_to_graphics_scene(lcc, gs2, gso2);
  for(unsigned int f=0; f<gs2.number_of_faces(); ++f)
  { check(!gs2.face_has_volume_color(f), "face with its own color"); }

  // A face outside any volume has no volume.
  GS gs3;
  gs3.face_begin();
  gs3.add_point_in_face(Point(0,0,0));
  gs3.add_point_in_face(Point(1,0,0));
  gs3.add_point_in_face(Point(0,1,0));
  gs3.face_end();
  check(gs3.number_of_faces()==1, "one face outside volumes");
  check(gs3.face_volume(0)==-1 && !gs3.face_has_volume_color(0),
        "face outside volumes has no volume");

  gs.clear();
  check(gs.number_of_faces()==0, "clear");

  std::cout<<"lcc_graphics_scene_shared_faces: OK"<<std::endl;
  return EXIT_SUCCESS;
}
