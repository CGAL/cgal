// Copyright (c) 2022 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of 3d-query-replace.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
// Prog that uses the query/replace method to generate hexaedral meshes.

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>

#include <map>
#include <sys/types.h>
#include <unistd.h>

#include "cmap_copy.h"
#include "cmap_query_replace.h"
#include "Compute_stats.h"
#include "init_to_preserve_for_query_replace.h"
#include "lcc_read_depending_extension.h"
#include "lcc_save_load_mesh.h"
#include "Print_txt.h"

////////////////////////////////////////////////////////////////////////////////
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT       FT;
typedef Kernel::Point_3  Point;
typedef Kernel::Vector_3 Vector;
////////////////////////////////////////////////////////////////////////////////
// Class My_vertex used to compute marching cube
template<typename Refs>
struct My_vertex: public CGAL::Cell_attribute_with_point<Refs>
{
  using Base=CGAL::Cell_attribute_with_point<Refs>;

  My_vertex():
    m_outside(false),
    m_new_vertex(false)
  {}

  My_vertex(const Point& p):
    Base(p),
    m_outside(false),
    m_new_vertex(false)
  {}

  bool m_outside; 
  bool m_new_vertex;
};
////////////////////////////////////////////////////////////////////////////////
class Myitems
{
public:
  template < class Refs >
  struct Dart_wrapper
  {
    typedef My_vertex<Refs> Vertex_attrib;
    typedef CGAL::Cell_attribute<Refs, unsigned int> Volume_attrib;
    typedef std::tuple<Vertex_attrib, void, void, Volume_attrib> Attributes;
  };
};
typedef CGAL::Linear_cell_complex_traits<3,Kernel> Mytraits;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3,Mytraits,Myitems> LCC3;

typedef typename LCC3::Dart_handle Dart_handle;
typedef typename LCC3::Vertex_attribute_handle Vertex_handle;
typedef typename LCC3::Attribute_handle<3>::type Volume_handle;
typedef typename LCC3::size_type   size_type;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef typename Kernel::Triangle_3 Triangle;
typedef typename Kernel::Segment_3 Segment;
typedef CGAL::Side_of_triangle_mesh<Polyhedron, Kernel> Side_of_mesh;
////////////////////////////////////////////////////////////////////////////////
[[ noreturn ]] void usage(int /*argc*/, char** argv)
{
  // Name
  std::cout<<"Name"<<std::endl;
  std::cout<<"        "<<argv[0]<<" - subdivides the given off file in hexahedra.";
  std::cout<<std::endl<<std::endl;
  // Synopsis
  std::cout<<"SYNOPSIS"<<std::endl;
  std::cout<<"        "<<argv[0]<<" [--help|-h|-?] "
           <<"[-create-transitions] [-draw] [-init X] [-lmax L] [-marching-cubes]"
           <<"[-no-remove-outside] [-no-signature] [-save] [-smooth] "
           <<"[-truncated-cubes] filename"
           <<std::endl<<std::endl;
  // Description
  std::cout<<"DESCRIPTION"<<std::endl;
  std::cout<<"        "<<" subdivides the given off file in hexahedra."
           <<std::endl<<std::endl;
  // Options
  std::cout<<"        --help, -h, -?"<<std::endl
           <<"                display this help and exit."
           <<std::endl<<std::endl;
  std::cout<<"        -create-transitions"<<std::endl
           <<"                use the 325 patterns to create transitions between hexahedra with 1 level of difference."
           <<std::endl<<std::endl;
  std::cout<<"        -draw"<<std::endl
           <<"                draw the final lcc."
           <<std::endl<<std::endl;
  std::cout<<"        -init X"<<std::endl
           <<"                fix the number of hexa in x axis to X (1 by default). Number of hexa in Y and Z are automatically computed to try to keep the ratio of the initial bounding box."
           <<std::endl<<std::endl;
  std::cout<<"        -lmax L"<<std::endl
           <<"                gives the maximum subdvision level (4 by default)."
           <<std::endl<<std::endl;
  std::cout<<"        -marching-cubes"<<std::endl
           <<"                create an approximation of the surface using marching cube method (imply !create-transition !no-remove-outside and !smooth)."
           <<std::endl<<std::endl;
  std::cout<<"        -no-remove-outside"<<std::endl
           <<"                do not remove outside hexa."
           <<std::endl<<std::endl;
  std::cout<<"        -no-signature"<<std::endl
           <<"                do not use signature to test isomorphism (but classical test)."
           <<std::endl<<std::endl;
  std::cout<<"        -save"<<std::endl
           <<"                save the lcc resulting of the subdvision (format mesh)."
           <<std::endl<<std::endl;
  std::cout<<"        -smooth"<<std::endl
           <<"                smooth the final volumic mesh. Contrary to marching-cubes, keep all the volumes instead of creating a surface."
           <<std::endl<<std::endl;
  std::cout<<"        -truncated-cubes"<<std::endl
           <<"                create a mesh with truncated cubes instead of cubes (imply all hexa with same level and !marching-cubes !smooth)."
           <<std::endl<<std::endl;
  exit(EXIT_FAILURE);
}
////////////////////////////////////////////////////////////////////////////////
[[ noreturn ]] void error_command_line(int argc, char** argv, const char* msg)
{
  std::cout<<"ERROR: "<<msg<<std::endl;
  usage(argc, argv);
}
////////////////////////////////////////////////////////////////////////////////
void process_command_line(int argc, char** argv,
                          std::string& filename,
                          bool& create_transitions,
                          bool& draw,
                          unsigned int& initX,
                          unsigned int& lmax,
                          bool& marching_cubes,
                          bool& no_remove_outside,
                          bool& no_signature,
                          bool& save,
                          bool& smooth,
                          bool& truncated_cubes
                          )
{
  filename="";
  create_transitions=false;
  draw=false;
  initX=1;
  lmax=4;
  marching_cubes=false;
  no_remove_outside=false;
  no_signature=false;
  save=false;
  smooth=false;
  truncated_cubes=false;

  bool helprequired=false;
  std::string arg;
  for (int i=1; i<argc; ++i)
  {
    arg=std::string(argv[i]);
    if(arg==std::string("-h") || arg==std::string("--help") || arg==std::string("-?"))
    { helprequired=true; }
    else if(arg=="-create-transitions")
    { create_transitions=true; }
    else if(arg=="-draw")
    { draw=true; }
    else if(arg=="-init")
    {
      if (argc-1-i<1)
      { error_command_line(argc, argv, "no numbers after -init"); }
      initX=std::stoi(std::string(argv[++i]));
    }
    else if(arg=="-lmax")
    {
      if (argc-1-i<1)
      { error_command_line(argc, argv, "no number after -lmax"); }
      lmax=std::stoi(std::string(argv[++i]));
    }
    else if(arg=="-marching-cubes")
    { marching_cubes=true; }
    else if(arg=="-no-remove-outside")
    { no_remove_outside=true; }
    else if(arg=="-no-signature")
    { no_signature=true; }
    else if(arg=="-save")
    { save=true; }
    else if(arg=="-smooth")
    { smooth=true; }
    else if(arg=="-truncated-cubes")
    { truncated_cubes=true; }
    else if(arg[0]=='-')
    { std::cout<<"Unknown option "<<arg<<", ignored."<<std::endl; }
    else { filename=arg; }
  }
  if (helprequired || filename.empty()) { usage(argc, argv); }

  if(truncated_cubes) { create_transitions=false; marching_cubes=false; smooth=false; }
  if(marching_cubes) { create_transitions=false; no_remove_outside=false; smooth=false; }
}
///////////////////////////////////////////////////////////////////////////////
/* To check the intersection between a particular square (given by its four
 * points) and the aabbtree t
*/
bool is_intersect(double x1, double y1, double z1,
                  double x2, double y2, double z2,
                  double x3, double y3, double z3,
                  double x4, double y4, double z4,
                  const Tree& t)
{
  Kernel::Point_3 p1(x1,y1,z1);
  Kernel::Point_3 p2(x2,y2,z2);
  Kernel::Point_3 p3(x3,y3,z3);
  Kernel::Point_3 p4(x4,y4,z4);

  // And compute the two triangles
  Triangle t1(p1, p2, p3);
  if(t.do_intersect(t1))
  { return true; }

  t1=Triangle(p1, p3, p4);
  if(t.do_intersect(t1))
  { return true; }

  return false;
}
///////////////////////////////////////////////////////////////////////////////
/* Test the intersection between a particular voxel (given by its two
 * extremal points) and the aabbtree t
*/
bool is_intersect(double x1, double y1, double z1,
                  double x2, double y2, double z2,
                  const Tree& t)
{
  return
      is_intersect(x1,y1,z1, x2,y1,z1, x2,y1,z2, x1,y1,z2, t) || // f1 y1
      is_intersect(x2,y2,z1, x2,y1,z1, x2,y1,z2, x2,y2,z2, t) || // f2 x2
      is_intersect(x1,y2,z1, x2,y2,z1, x2,y2,z2, x1,y2,z2, t) || // f3 y2
      is_intersect(x1,y1,z1, x1,y1,z2, x1,y2,z2, x1,y2,z1, t) || // f4 x1
      is_intersect(x1,y1,z1, x2,y1,z1, x2,y2,z1, x1,y2,z1, t) || // f5 z1
      is_intersect(x1,y1,z2, x2,y1,z2, x2,y2,z2, x1,y2,z2, t);   // f6 z2
}
///////////////////////////////////////////////////////////////////////////////
/* Test if a particular point is outside of the object (Tree), knowing there is
 * no intersection between its voxel and the tree.
 */
bool is_outside_knowing_no_intersect(const Kernel::Point_3& p,
                                     const Tree& t)
{
  Side_of_mesh s(t);
  CGAL::Bounded_side res=s(p);
  return res!=CGAL::ON_BOUNDED_SIDE; // && !=CGAL::ON_BOUNDARY ?
}
///////////////////////////////////////////////////////////////////////////////
/* Test if a particular voxel (given by its two extremal points) is outside
   * of the object (Tree).
   */
bool is_outside(double x1, double y1, double z1,
                double x2, double y2, double z2,
                const Tree& t)
{
  if (is_intersect(x1, y1, z1, x2, y2, z2, t)) { return false; }
  return is_outside_knowing_no_intersect(Kernel::Point_3(x1, y1, z1), t);
}
///////////////////////////////////////////////////////////////////////////////
/* Test the intersection between a particular voxel of the LCC3
 * and the Tree (variable t - the object - created from the load of a surface mesh file).
 * (Rq: the variable lcc represents a set of voxels - a block of voxels of the grid).
 */
bool is_intersect(LCC3& lcc, Dart_handle dh, const Tree& t)
{
  CGAL::Bbox_3 bbox=lcc.point(dh).bbox();
  // For each vertex of the volume
  for(auto it=lcc.one_dart_per_incident_cell<0,3>(dh).begin(),
      itend=lcc.one_dart_per_incident_cell<0,3>(dh).end(); it!=itend; ++it)
  { bbox+=lcc.point(it).bbox(); }

  return is_intersect(bbox.xmin(), bbox.ymin(), bbox.zmin(),
                      bbox.xmax(), bbox.ymax(), bbox.zmax(), t);
}
///////////////////////////////////////////////////////////////////////////////
/* Test if a particular voxel of the LCC3 (lcc) is outside of the object (Tree).
 * (Rq: the variable lcc represents a set of voxels - a block of voxels of the grid).
 */
bool is_outside(LCC3& lcc, Dart_handle dh, const Tree& t)
{
  CGAL::Bbox_3 bbox=lcc.point(dh).bbox();
  // For each vertex of the volume
  for(auto it=lcc.one_dart_per_incident_cell<0,3>(dh).begin(),
      itend=lcc.one_dart_per_incident_cell<0,3>(dh).end(); it!=itend; ++it)
  { bbox+=lcc.point(it).bbox(); }

  return is_outside(bbox.xmin(), bbox.ymin(), bbox.zmin(),
                    bbox.xmax(), bbox.ymax(), bbox.zmax(), t);
}
////////////////////////////////////////////////////////////////////////////////
void count_subdivisions(LCC3& lcc, std::map<std::size_t, std::size_t>& levels)
{
  levels.clear();
  unsigned int mylevel;
  for(auto it=lcc.attributes<3>().begin(), itend=lcc.attributes<3>().end();
      it!=itend; ++it)
  {
    mylevel=lcc.info_of_attribute<3>(it);
    auto l=levels.find(mylevel);
    if(l==levels.end())
    { levels[mylevel]=1; }
    else
    { ++(l->second); }
  }
  for(auto it: levels)
  { std::cout<<it.first<<"->"<<it.second<<"  "; }
  std::cout<<std::endl;
}
void count_subdivisions(LCC3& lcc)
{
  std::map<std::size_t, std::size_t> levels;
  count_subdivisions(lcc, levels);
}
///////////////////////////////////////////////////////////////////////////////
LCC3::Dart_handle make_one_hexa(LCC3& lcc, double x1, double y1, double z1,
                                double x2, double y2, double z2)
{
  LCC3::Dart_handle dh = lcc.make_hexahedron
    (LCC3::Point(x1, y1, z1), LCC3::Point(x2, y1, z1),
     LCC3::Point(x2, y2, z1),
     LCC3::Point(x1, y2, z1),
     LCC3::Point(x1, y2, z2),
     LCC3::Point(x1, y1, z2),
     LCC3::Point(x2, y1, z2),
     LCC3::Point(x2, y2, z2));
  return dh;
}
///////////////////////////////////////////////////////////////////////////////
/// Create a block of (nbx x nby x nbz) hexa. External hexa are not created.
/// First hexa starts at position (startx, starty, starty) and each hexa has size
/// (sizex, sizey, sizez)
void create_hexa_that_intersect(LCC3& lcc,
                                const Tree& t,
                                unsigned int nbx, unsigned int nby, unsigned nbz,
                                double startx, double starty, double startz,
                                double sizex, double sizey, double sizez)
{
  Dart_handle dh;
  double x1, y1, z1, x2, y2, z2;
  for (unsigned int x=0; x<nbx; ++x)
  {
    for (unsigned int y=0; y<nby; ++y)
    {
      for (unsigned int z=0; z<nbz; ++z)
      {
        x1=startx+x*sizex; y1=starty+y*sizey; z1=startz+z*sizez;
        x2=startx+(x+1)*sizex; y2=starty+(y+1)*sizey; z2=startz+(z+1)*sizez;
        if (!is_outside(x1, y1, z1, x2, y2, z2, t))
        {
          dh=make_one_hexa(lcc,x1, y1, z1, x2, y2, z2);
          lcc.set_attribute<3>(dh, lcc.create_attribute<3>(0));
        }
      }
    }
  }
  lcc.sew3_same_facets();
}
///////////////////////////////////////////////////////////////////////////////
void compute_size_and_init(unsigned int init, const Tree& aabb_tree,
                           int& longestAxis,
                           double& sx, double& sy, double& sz,
                           unsigned int& initX,
                           unsigned int& initY,
                           unsigned int& initZ)
{
  sx=(aabb_tree.bbox().xmax()-aabb_tree.bbox().xmin());
  sy=(aabb_tree.bbox().ymax()-aabb_tree.bbox().ymin());
  sz=(aabb_tree.bbox().zmax()-aabb_tree.bbox().zmin());

  if (sx>=sy && sx>=sz)     longestAxis = 0;
  else if(sy>=sx && sy>=sz) longestAxis = 1;
  else                      longestAxis = 2;

  if(longestAxis == 0)
  {
    initX = init;
    initY = std::ceil(init * (sy / sx));
    initZ = std::ceil(init * (sz / sx));
  }
  else if(longestAxis == 1)
  {
    initX = std::ceil(init * (sx / sy));
    initY = init;
    initZ = std::ceil(init * (sz / sy));
  }
  else
  {
    initX = std::ceil(init * (sx / sz));
    initY = std::ceil(init * (sy / sz));
    initZ = init;
  }

  sx/=initX;
  sy/=initY;
  sz/=initZ;
}
///////////////////////////////////////////////////////////////////////////////
bool create_initial_voxels(const std::string& file,
                           LCC3& lcc,
                           unsigned int init,
                           CGAL::Polyhedron_3<Kernel>& surface,
                           Tree& aabb_tree)
{
  // 1) Load surface mesh
  std::ifstream off_file(file);
  if(!off_file.good()) return false; // File open error

  std::chrono::system_clock::time_point start2=std::chrono::system_clock::now();

  surface.clear();
  off_file>>surface;
  CGAL::Polygon_mesh_processing::triangulate_faces(surface);

  std::chrono::duration<double> diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Load off: ", diff, "s.");

  // 2) Compute AABB tree
  start2 = std::chrono::system_clock::now();

  aabb_tree.insert(faces(surface).first, faces(surface).second, surface);
  aabb_tree.accelerate_distance_queries();
  aabb_tree.bbox();

  diff = std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Compute AABB tree: ", diff, "s.");

  // std::cout<<"AABB tree bbox="<<aabb_tree.bbox()<<std::endl;

  double sx, sy, sz;
  int longestAxis;    // 0 = x, 1 = y, 2 = z
  unsigned int initX, initY, initZ;
  compute_size_and_init(init, aabb_tree, longestAxis, sx, sy, sz,
                        initX, initY, initZ);

  start2 = std::chrono::system_clock::now();

  // 3) Create a block of hexahedra
  create_hexa_that_intersect(lcc, aabb_tree,
                             initX, initY, initZ,
                             aabb_tree.bbox().xmin(),
                             aabb_tree.bbox().ymin(),
                             aabb_tree.bbox().zmin(),
                             sx, sy, sz);

  diff = std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Create block of hexa: ", diff, "s.");

  return true;
}
////////////////////////////////////////////////////////////////////////////////
void subdivide_edges(LCC3& lcc,
                     const std::vector<Dart_handle>& hexa_tosubdivide,
                     size_type corner_mark)
{
  size_type amark=lcc.get_new_mark();
  std::vector<Dart_handle> edges_to_subdivide;
  edges_to_subdivide.reserve(hexa_tosubdivide.size()); // a la louche ;)

  for(Dart_handle its: hexa_tosubdivide)
  {
    for(auto itd=lcc.darts_of_cell<3>(its).begin(),
        itdend=lcc.darts_of_cell<3>(its).end(); itd!=itdend; ++itd)
    {
      if(!lcc.is_marked(itd, amark) && lcc.is_marked(itd, corner_mark) &&
         lcc.is_marked(lcc.beta<2>(itd), corner_mark))
      {
        lcc.mark_cell<1>(itd, amark);
        edges_to_subdivide.push_back(itd);
      }
    }
  }

  for(Dart_handle itd: edges_to_subdivide)
  {
    lcc.unmark_cell<1>(itd, amark);
    lcc.insert_barycenter_in_cell<1>(itd);
  }

  assert(lcc.is_whole_map_unmarked(amark));
  lcc.free_mark(amark);
}
////////////////////////////////////////////////////////////////////////////////
void subdivide_faces(LCC3& lcc,
                     const std::vector<Dart_handle>& hexa_tosubdivide,
                     size_type corner_mark,
                     Pattern_substituer<LCC3>& ps,
                     std::size_t& nb_replaced)
{
  size_type amark=lcc.get_new_mark();
  std::vector<Dart_handle> faces_to_subdivide;
  faces_to_subdivide.reserve(hexa_tosubdivide.size()); // a la louche ;)
  for(Dart_handle its: hexa_tosubdivide)
  {
    for(auto itd=lcc.darts_of_cell<3>(its).begin(),
        itdend=lcc.darts_of_cell<3>(its).end(); itd!=itdend; ++itd)
    {
      if(!lcc.is_marked(itd, amark) && lcc.is_marked(itd, corner_mark) &&
         lcc.beta<1,1,1,1>(itd)!=itd)
      {
        lcc.mark_cell<2>(itd, amark);
        faces_to_subdivide.push_back(itd);
      }
    }
  }

  for(Dart_handle itd: faces_to_subdivide)
  {
    lcc.unmark_cell<2>(itd, amark);   
    /*if(ps.query_replace_one_face(lcc, itd, corner_mark)==
       std::numeric_limits<std::size_t>::max())
    { std::cout<<"[ERROR] in subdivide_faces: one query/replace failed."
               <<std::endl;
      / *LCC3 lcc_error;
      copy_cell<3>(lcc, itd, lcc_error);
      for(auto itnewv=lcc.one_dart_per_incident_cell<2,3>(itd).begin(),
          itnewvend=lcc.one_dart_per_incident_cell<2,3>(itd).end(); itnewv!=itnewvend; ++itnewv)
      {
        if(!lcc.is_free<3>(itnewv))
        { copy_cell<3>(lcc, lcc.beta<3>(itnewv), lcc_error); }
      }
      CGAL::draw(lcc_error);
      * /
    }*/
    ps.replace_one_face_from_dart(lcc, itd, ps.m_fpatterns[0],
                                  ps.m_fsignatures.begin()->second.first);
    ++nb_replaced;
  }

  assert(lcc.is_whole_map_unmarked(amark));
  lcc.free_mark(amark);
}
////////////////////////////////////////////////////////////////////////////////
void subdivide_volumes(LCC3& lcc,
                       const std::vector<Dart_handle>& hexa_tosubdivide,
                       size_type corner_mark,
                       Pattern_substituer<LCC3>& ps,
                       std::size_t& nb_replaced)
{
  std::vector<Dart_handle> vol_darts;
  vol_darts.reserve(96);
  Dart_handle corner_dart=nullptr;

  for(Dart_handle its: hexa_tosubdivide)
  {
    vol_darts.clear();
    corner_dart=nullptr;
    for(auto it=lcc.darts_of_cell<3>(its).begin(),
        itend=lcc.darts_of_cell<3>(its).end(); it!=itend; ++it)
    {
      if(!lcc.is_marked(it, corner_mark)) { vol_darts.push_back(it); }
      else if(corner_dart==nullptr) { corner_dart=it; }
    }

    lcc.info<3>(its)=1+lcc.info<3>(its);
    /* if(ps.query_replace_one_volume(lcc, its, corner_mark)==
       std::numeric_limits<std::size_t>::max())
    { std::cout<<"[ERROR] in subdivide_volumes: one query/replace failed."
              <<std::endl; } */
    ps.replace_one_volume_from_dart(lcc, corner_dart, ps.m_vpatterns[0],
                                    ps.m_vsignatures.begin()->second.first);
    ++nb_replaced;

    // Mark all old darts of the 8 volume as corners
    for(Dart_handle itv: vol_darts)
    {
      if(!lcc.is_marked(itv, corner_mark))
      {
        for(auto itnewv=lcc.darts_of_cell<3>(itv).begin(),
            itnewvend=lcc.darts_of_cell<3>(itv).end(); itnewv!=itnewvend; ++itnewv)
        { lcc.mark(itnewv, corner_mark); }
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void subdivide_hexa(LCC3& lcc, Dart_handle dh, LCC3::size_type corner_mark,
                    Pattern_substituer<LCC3>& ps)
{
  size_type mark_to_preserve_v=lcc.get_new_mark();
  mark_volume_corners<LCC3>(lcc, dh, mark_to_preserve_v);
  std::vector<Dart_handle> tounmark;
  //size_type mark_to_preserve_s=lcc.get_new_mark();
  std::unordered_set<Vertex_handle> vertices;
  std::vector<Dart_handle> one_dart_per_corner;
  one_dart_per_corner.reserve(8);

  lcc.info<3>(dh)=1+lcc.info<3>(dh);

  // 1) Subdivide each edge which is not yet subdivide
  std::vector<Dart_handle> tosubdivide;
  for(auto it=lcc.darts_of_cell<3>(dh).begin(),
      itend=lcc.darts_of_cell<3>(dh).end(); it!=itend; ++it)
  {
    if(it<lcc.beta<2>(it)) // considere only one dart per edge of the hexa
    {
      if(lcc.is_marked(it, corner_mark) && lcc.is_marked(lcc.beta<2>(it), corner_mark))
      { tosubdivide.push_back(it); }
    }
    if(lcc.is_marked(it, corner_mark) && vertices.count(lcc.vertex_attribute(it))==0)
    {
      one_dart_per_corner.push_back(it);
      vertices.insert(lcc.vertex_attribute(it));
    }
    tounmark.push_back(it);
  }
  assert(one_dart_per_corner.size()==8);

  for(auto it: tosubdivide)
  { lcc.insert_barycenter_in_cell<1>(it); }

  // 2) Subdivide faces.
  tosubdivide.clear();
  for(auto it=lcc.one_dart_per_incident_cell<2,3>(dh).begin(),
      itend=lcc.one_dart_per_incident_cell<2,3>(dh).end(); it!=itend; ++it)
  { tosubdivide.push_back(it); }

  for(auto it: tosubdivide)
  { ps.query_replace_one_face(lcc, it, mark_to_preserve_v); }

  // 3) Subdivide the hexa itself
  ps.query_replace_one_volume(lcc, dh, mark_to_preserve_v);

  // 4) Mark all darts of the 8 volume as corners
  for(auto itv: one_dart_per_corner)
  {
    std::size_t nb=0;
    for(auto it=lcc.darts_of_cell<3>(itv).begin(),
        itend=lcc.darts_of_cell<3>(itv).end(); it!=itend; ++it)
    { lcc.mark(it, corner_mark); ++nb; }
    assert(nb==24);
  }

  for(auto it: tounmark)
  { lcc.unmark(it, mark_to_preserve_v); }
  assert(lcc.is_whole_map_unmarked(mark_to_preserve_v));
  lcc.free_mark(mark_to_preserve_v);

  // assert(lcc.is_valid());
}
////////////////////////////////////////////////////////////////////////////////
std::size_t compute_hexa_to_subdivide_for_truncated_cubes
(LCC3& lcc,
 const Tree& aabb_tree,
 bool no_remove_outside,
 std::vector<std::vector<Dart_handle>>& hexa_tosubdivide)
{
  std::size_t res=0;
  for(auto& it: hexa_tosubdivide)
  { it.clear(); }

  std::queue<Volume_handle> totreat;
  unsigned int mylevel;
  for (auto it=lcc.attributes<3>().begin(); it!=lcc.attributes<3>().end(); ++it)
  {
    mylevel=lcc.info_of_attribute<3>(it);
    if(is_intersect(lcc, lcc.dart_of_attribute<3>(it), aabb_tree))
    {
      totreat.push(it);
      hexa_tosubdivide[mylevel].push_back(lcc.dart_of_attribute<3>(it));
      ++res;
    }
    else if(!no_remove_outside &&
            is_outside_knowing_no_intersect
            (lcc.point(lcc.dart_of_attribute<3>(it)), aabb_tree))
    { lcc.remove_cell<3>(lcc.dart_of_attribute<3>(it)); }
    else // Hexa inside
    {
      totreat.push(it);
      hexa_tosubdivide[mylevel].push_back(lcc.dart_of_attribute<3>(it));
      ++res;
    }
  }
  return res;
}
////////////////////////////////////////////////////////////////////////////////
std::size_t compute_hexa_to_subdivide(LCC3& lcc,
                                      const Tree& aabb_tree,
                                      unsigned int curlevel,
                                      bool no_remove_outside,
                                      bool marching_cubes,
                                      std::vector<std::vector<Dart_handle>>&
                                      hexa_tosubdivide)
{
  std::size_t res=0;

  for(auto& it: hexa_tosubdivide)
  { it.clear(); }

  std::queue<Volume_handle> totreat;
  unsigned int mylevel;
  size_type treated=lcc.get_new_mark();
  for (auto it=lcc.attributes<3>().begin(); it!=lcc.attributes<3>().end(); ++it)
  {
    mylevel=lcc.info_of_attribute<3>(it);
    if (mylevel==curlevel)
    {
      if(is_intersect(lcc, lcc.dart_of_attribute<3>(it), aabb_tree))
      {
        totreat.push(it);
        assert(lcc.dart_of_attribute<3>(it)!=nullptr);
        hexa_tosubdivide[mylevel].push_back(lcc.dart_of_attribute<3>(it));
        ++res;
        lcc.mark_cell<3>(lcc.dart_of_attribute<3>(it), treated);
      }
      else if(marching_cubes ||
              (!no_remove_outside &&
               is_outside_knowing_no_intersect
               (lcc.point(lcc.dart_of_attribute<3>(it)), aabb_tree)))
      { lcc.remove_cell<3>(lcc.dart_of_attribute<3>(it)); }
    }
  }
  if(!marching_cubes)
  {
    Volume_handle curvh;
    while(!totreat.empty())
    {
      curvh=totreat.front();
      totreat.pop();
      mylevel=lcc.info_of_attribute<3>(curvh);
      // Test all neighboors of curvh
      for(auto it=lcc.darts_of_cell<3>(lcc.dart_of_attribute<3>(curvh)).begin(),
          itend=lcc.darts_of_cell<3>(lcc.dart_of_attribute<3>(curvh)).end();
          it!=itend; ++it)
      {
        // Adjacency through faces AND edges !
        for(auto ite=lcc.darts_of_cell<1>(it).begin(),
            iteend=lcc.darts_of_cell<1>(it).end(); ite!=iteend; ++ite)
        {
          if(!lcc.is_marked(ite, treated) && mylevel==1+lcc.info<3>(ite))
          {
            totreat.push(lcc.attribute<3>(ite));
            hexa_tosubdivide[lcc.info<3>(ite)].push_back(ite);
            ++res;
            lcc.mark_cell<3>(ite, treated);
          }
        }
      }
    }
  }
  lcc.free_mark(treated);
  return res;
}
////////////////////////////////////////////////////////////////////////////////
void compute_marching_cubes(LCC3& lcc, size_type corner_mark,
                            const Tree& aabb_tree, bool no_signature)
{
  std::chrono::system_clock::time_point start=std::chrono::system_clock::now();
  std::chrono::system_clock::time_point start2=std::chrono::system_clock::now();
  Pattern_substituer<LCC3> ps2; // pattern substituer to create transitions
  
  ps2.load_vpatterns(CGAL::data_file_path("query_replace/marching-cubes/vpattern/"),
                     mark_vpattern_corners<LCC3>);

  ps2.load_fpatterns(CGAL::data_file_path("query_replace/marching-cubes/fpattern/"),
                     mark_fpattern_corners<LCC3>);
  std::chrono::duration<double> diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Load transition patterns: ", diff, "s.");

  start2 = std::chrono::system_clock::now();
  Side_of_mesh s(aabb_tree);
  for (auto it=lcc.attributes<0>().begin(); it!=lcc.attributes<0>().end(); ++it)
  {
    it->m_outside=(s(lcc.point_of_vertex_attribute(it))!=
        CGAL::ON_BOUNDED_SIDE);
  }
  diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Compute inside/outside: ", diff, "s.");

  // 1) Split edges
  start2 = std::chrono::system_clock::now();
  std::vector<Dart_handle> to_subdivide;
  to_subdivide.reserve(lcc.number_of_attributes<3>()); // a la louche ;)
  std::size_t nb_replaced_f=0;
  std::size_t nb_replaced_v=0;
  
  for(auto itd=lcc.one_dart_per_cell<1>().begin(),
      itdend=lcc.one_dart_per_cell<1>().end(); itd!=itdend; ++itd)
  {
    if(lcc.vertex_attribute(itd)->m_outside!=
       lcc.vertex_attribute(lcc.other_extremity(itd))->m_outside)
    { to_subdivide.push_back(itd); }
  }

  for(Dart_handle itd: to_subdivide)
  {
    FT d1=std::sqrt(aabb_tree.squared_distance(lcc.point(itd)));
    FT d2=std::sqrt(aabb_tree.squared_distance(lcc.point(lcc.other_extremity(itd))));
    Vector v(lcc.point(itd), lcc.point(lcc.other_extremity(itd)));
    FT lg=std::sqrt(v.squared_length());
    Point p1=lcc.point(itd);
    Point p2=lcc.point(lcc.other_extremity(itd));
    p1=p1.transform(Kernel::Aff_transformation_3(CGAL::TRANSLATION, (v*d1)/lg));
    p2=p2.transform(Kernel::Aff_transformation_3(CGAL::TRANSLATION, (-v*d2)/lg));
    auto newv=lcc.insert_point_in_cell<1>(itd, Point((p1.x()+p2.x())/2,
                                                     (p1.y()+p2.y())/2,
                                                     (p1.z()+p2.z())/2));
    lcc.vertex_attribute(newv)->m_outside=false; // new vertices belong to the surface => no outside
    lcc.vertex_attribute(newv)->m_new_vertex=true;
    // Complicated computation; but try to deal with the case where the surface
    // intersect the edge several times.
  }
  diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Split edges: ", diff, "s.");

  // 2) Replace faces
  to_subdivide.clear();
  start2=std::chrono::system_clock::now();
  for(auto itd=lcc.one_dart_per_cell<2>().begin(),
      itdend=lcc.one_dart_per_cell<2>().end(); itd!=itdend; ++itd)
  {
    if(lcc.beta<1,1,1,1>(itd)!=itd)
    { to_subdivide.push_back(itd); }
  }
  if(no_signature)
  {
    for(auto itf: to_subdivide)
    {
      ps2.query_replace_one_face_without_signature(lcc, itf, corner_mark);
      ++nb_replaced_f;
    }
  }
  else
  {
    for(auto itf: to_subdivide)
    {
      ps2.query_replace_one_face(lcc, itf, corner_mark);
      ++nb_replaced_f;
    }
  }
  diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Replace faces: ", diff, "s.");

  auto markfaces=lcc.get_new_mark();
  lcc.negate_mark(markfaces); // All darts of voxel faces are marked

  // 3) Replace volumes
  to_subdivide.clear();
  start2=std::chrono::system_clock::now();
  for(auto itv=lcc.attributes<3>().begin(),
      itvend=lcc.attributes<3>().end(); itv!=itvend; ++itv)
  {
    if(!lcc.is_volume_combinatorial_hexahedron(lcc.dart_of_attribute<3>(itv)))
    { to_subdivide.push_back(lcc.dart_of_attribute<3>(itv)); }
  }
  if(no_signature)
  {
    for(auto itv: to_subdivide)
    {
      ps2.query_replace_one_volume_without_signature(lcc, itv, corner_mark);
      ++nb_replaced_v;
    }
  }
  else
  {
    for(auto itv: to_subdivide)
    {
      ps2.query_replace_one_volume(lcc, itv, corner_mark);
      ++nb_replaced_v;
    }
  }
  diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Replace volumes: ", diff, "s.");

  // 4) Remove external faces
  start2=std::chrono::system_clock::now();
  lcc.set_automatic_attributes_management(false);
  for(auto it=lcc.darts().begin(); it!=lcc.darts().end(); ++it)
  {
    if(lcc.is_marked(it, markfaces))
    {
      bool newface=true;
      Dart_handle curdh=it;
      do
      {
        if(!lcc.vertex_attribute(curdh)->m_new_vertex)
        { newface=false; }
        else
        { curdh=lcc.beta<1>(curdh); }
      }
      while(newface && curdh!=it);
      if(!newface)
      { lcc.remove_cell<2>(it); }
      else
      { lcc.unmark_cell<2>(it, markfaces); }
    }
  }
  assert(lcc.is_whole_map_unmarked(markfaces));
  lcc.free_mark(markfaces);

  lcc.set_automatic_attributes_management(true);

  diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Remove outside faces: ", diff, "s.");

  // End
  diff=std::chrono::system_clock::now()-start;
  print_txt_with_endl("Compute marching cube TOTAL: ", diff, "s.");

  std::cout<<"#faces-replaced: "<<nb_replaced_f<<", #volumes-replaced: "<<nb_replaced_v<<std::endl;
}
////////////////////////////////////////////////////////////////////////////////
void smooth_mesh(LCC3& lcc, size_type corner_mark, const Tree& aabb_tree,
                 bool no_signature)
{
  std::chrono::system_clock::time_point start=std::chrono::system_clock::now();
  std::chrono::system_clock::time_point start2=std::chrono::system_clock::now();
  Pattern_substituer<LCC3> ps2; // pattern substituer to create transitions
  ps2.load_vpatterns(CGAL::data_file_path("query_replace/marching-cubes/vpattern/"),
                     mark_vpattern_corners<LCC3>);
  ps2.load_fpatterns(CGAL::data_file_path("query_replace/marching-cubes/fpattern/"),
                     mark_fpattern_corners<LCC3>);
  std::chrono::duration<double> diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Load transition patterns: ", diff, "s.");

  start2 = std::chrono::system_clock::now();
  Side_of_mesh s(aabb_tree);
  for (auto it=lcc.attributes<0>().begin(); it!=lcc.attributes<0>().end(); ++it)
  {
    it->m_outside=(s(lcc.point_of_vertex_attribute(it))!=
                   CGAL::ON_BOUNDED_SIDE);
  }
  diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Compute inside/outside: ", diff, "s.");

  bool mark_to_subdivide=lcc.get_new_mark();
  // 1) Split edges
  start2 = std::chrono::system_clock::now();
  size_type amark=lcc.get_new_mark();
  std::vector<Dart_handle> to_subdivide;
  to_subdivide.reserve(lcc.number_of_attributes<3>()); // a la louche ;)
  std::size_t nb_replaced_f=0;
  std::size_t nb_replaced_v=0;
  
  for(auto itd=lcc.darts().begin(), itdend=lcc.darts().end(); itd!=itdend; ++itd)
  {
    if(!lcc.is_marked(itd, amark) &&
       lcc.vertex_attribute(itd)->m_outside!=
       lcc.vertex_attribute(lcc.other_extremity(itd))->m_outside)
    {
      for(auto ite=lcc.darts_of_cell_basic<1>(itd, amark).begin(),
          iteend=lcc.darts_of_cell_basic<1>(itd, amark).end(); ite!=iteend; ++ite)
      {
        lcc.mark(ite, amark);
        if(!lcc.is_marked(ite, mark_to_subdivide) &&
           lcc.is_volume_combinatorial_hexahedron(ite))
        { lcc.mark_cell<3>(ite, mark_to_subdivide); }
      }
      to_subdivide.push_back(itd);
    }
  }

  for(Dart_handle itd: to_subdivide)
  {
    lcc.unmark_cell<1>(itd, amark);
    FT d1=std::sqrt(aabb_tree.squared_distance(lcc.point(itd)));
    FT d2=std::sqrt(aabb_tree.squared_distance(lcc.point(lcc.other_extremity(itd))));
    Vector v(lcc.point(itd), lcc.point(lcc.other_extremity(itd)));
    FT lg=std::sqrt(v.squared_length());
    Point p1=lcc.point(itd);
    Point p2=lcc.point(lcc.other_extremity(itd));
    p1=p1.transform(Kernel::Aff_transformation_3(CGAL::TRANSLATION, (v*d1)/lg));
    p2=p2.transform(Kernel::Aff_transformation_3(CGAL::TRANSLATION, (-v*d2)/lg));
    auto newv=lcc.insert_point_in_cell<1>(itd, Point((p1.x()+p2.x())/2,
                                                     (p1.y()+p2.y())/2,
                                                     (p1.z()+p2.z())/2));
    lcc.vertex_attribute(newv)->m_outside=false; // new vertices belong to the surface => no outside
    lcc.vertex_attribute(newv)->m_new_vertex=true;
    // Complicated computation; but try to deal with the case where the surface
    // intersect the edge several times.
  }
  diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Split edges: ", diff, "s.");
  assert(lcc.is_whole_map_unmarked(amark));

  // 2) Replace faces
  to_subdivide.clear();
  start2=std::chrono::system_clock::now();
  for(auto itd=lcc.darts().begin(), itdend=lcc.darts().end();
      itd!=itdend; ++itd)
  {
    if(!lcc.is_marked(itd, amark) && lcc.is_marked(itd, mark_to_subdivide) &&
       lcc.beta<1,1,1,1>(itd)!=itd)
    {
      lcc.mark_cell<2>(itd, amark);
      to_subdivide.push_back(itd);
    }
  }
  if(no_signature)
  {
    for(auto itf: to_subdivide)
    {
      lcc.unmark_cell<2>(itf, amark);
      ps2.query_replace_one_face_without_signature(lcc, itf, corner_mark);
      ++nb_replaced_f;
    }
  }
  else
  {
    for(auto itf: to_subdivide)
    {
      lcc.unmark_cell<2>(itf, amark);
      ps2.query_replace_one_face(lcc, itf, corner_mark);
      ++nb_replaced_f;
    }
  }
  diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Replace faces: ", diff, "s.");
  assert(lcc.is_whole_map_unmarked(amark));

  // 3) Replace volumes
  to_subdivide.clear();
  start2=std::chrono::system_clock::now();
  for(auto itd=lcc.darts().begin(), itdend=lcc.darts().end();
      itd!=itdend; ++itd)
  {
    if(lcc.is_marked(itd, mark_to_subdivide))
    {
      for(auto itv=lcc.darts_of_cell<3>(itd).begin(),
          itvend=lcc.darts_of_cell<3>(itd).end(); itv!=itvend; ++itv)
      { lcc.unmark(itv, mark_to_subdivide); }
      if(!lcc.is_volume_combinatorial_hexahedron(itd))
      { to_subdivide.push_back(itd); }
    }
  }
  if(no_signature)
  {
    for(auto itv: to_subdivide)
    {
      ps2.query_replace_one_volume_without_signature(lcc, itv, corner_mark);
      ++nb_replaced_v;
    }
  }
  else
  {
    for(auto itv: to_subdivide)
    {
      ps2.query_replace_one_volume(lcc, itv, corner_mark);
      ++nb_replaced_v;
    }
  }
  diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Replace volumes: ", diff, "s.");
  assert(lcc.is_whole_map_unmarked(amark));
  lcc.free_mark(amark);
  assert(lcc.is_whole_map_unmarked(mark_to_subdivide));
  lcc.free_mark(mark_to_subdivide);

  // 4) Remove external faces
  start2=std::chrono::system_clock::now();
  auto markfaces=lcc.get_new_mark();
  lcc.set_automatic_attributes_management(false);
  for(auto it=lcc.darts().begin(); it!=lcc.darts().end(); ++it)
  {
    if(!lcc.is_marked(it, markfaces))
    {
      bool outside=false;
      Dart_handle curdh=it;
      do
      {
        if(lcc.vertex_attribute(curdh)->m_outside &&
           !lcc.vertex_attribute(curdh)->m_new_vertex)
        { outside=true; }
        lcc.mark(curdh, markfaces);
        if(!lcc.is_free<3>(curdh))
        { lcc.mark(lcc.beta<3>(curdh), markfaces); }
        curdh=lcc.beta<1>(curdh);
      }
      while(curdh!=it);
      if(outside)
      { lcc.remove_cell<2>(it); }
    }
  }
  assert(lcc.is_whole_map_marked(markfaces));
  lcc.free_mark(markfaces);

  lcc.set_automatic_attributes_management(true);

  diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Remove outside faces: ", diff, "s.");

  // End
  diff=std::chrono::system_clock::now()-start;
  print_txt_with_endl("Compute marching cube TOTAL: ", diff, "s.");

  std::cout<<"#faces-replaced: "<<nb_replaced_f<<", #volumes-replaced: "<<nb_replaced_v<<std::endl;
}
////////////////////////////////////////////////////////////////////////////////
void create_transition_patterns(LCC3& lcc, size_type corner_mark,
                                bool no_signature)
{
  std::chrono::system_clock::time_point start=std::chrono::system_clock::now();
  std::chrono::system_clock::time_point start2=std::chrono::system_clock::now();
  Pattern_substituer<LCC3> ps2; // pattern substituer to create transitions
  ps2.load_vpatterns(CGAL::data_file_path("query_replace/hexa-325-patterns/"),
                     mark_vpattern_corners<LCC3>);
  ps2.load_fpatterns(CGAL::data_file_path("query_replace/square-5-patterns"),
                     mark_fpattern_corners<LCC3>);
  std::chrono::duration<double> diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Load transition patterns: ", diff, "s.");

  start2 = std::chrono::system_clock::now();
  std::vector<Dart_handle> faces;
  std::vector<Dart_handle> volumes;
  std::size_t nb_replaced_f=0;
  std::size_t nb_replaced_v=0;
  
  cell_topo t;
  for(auto it=lcc.attributes<3>().begin(); it!=lcc.attributes<3>().end(); ++it)
  {
    t=get_cell_topo(lcc, lcc.dart_of_attribute<3>(it));
    if(t==GENERIC_3D)
    { volumes.push_back(lcc.dart_of_attribute<3>(it)); }
  }
  for(auto curdh: volumes)
  {
    faces.clear();
    for(auto itf=lcc.one_dart_per_incident_cell<2,3>(curdh).begin(),
        itfend=lcc.one_dart_per_incident_cell<2,3>(curdh).end(); itf!=itfend; ++itf)
    { faces.push_back(itf); }
    if(no_signature)
    {
      for(auto itf: faces)
      {
        ps2.query_replace_one_face_without_signature(lcc, itf, corner_mark);
        ++nb_replaced_f;
      }
    }
    else
    {
      for(auto itf: faces)
      {
        ps2.query_replace_one_face(lcc, itf, corner_mark);
        ++nb_replaced_f;
      }
    }

    std::size_t res;
    if(no_signature)
    {
      res=ps2.query_replace_one_volume_without_signature(lcc, curdh, corner_mark);
      ++nb_replaced_v;
    }
    else
    {
      res=ps2.query_replace_one_volume(lcc, curdh, corner_mark);
      ++nb_replaced_v;
    }
    if(res==std::numeric_limits<std::size_t>::max())
    {
      std::cout<<"[ERROR] transition not found for volume "
                    <<lcc.attributes<3>().index(lcc.attribute<3>(curdh))
                   <<std::endl;
      static unsigned int nberror=0;
      LCC3 lcc_error;
      copy_cell<3>(lcc, curdh, lcc_error);
      CGAL::draw(lcc_error);
      save_object_3D(std::string("error"+std::to_string(nberror++)+
                                 ".mesh"), lcc_error);
    }
  }
  diff=std::chrono::system_clock::now()-start;
  print_txt_with_endl("Create transitions TOTAL: ", diff, "s.");

  std::cout<<"#faces-replaced: "<<nb_replaced_f<<", #volumes-replaced: "<<nb_replaced_v<<std::endl;
}
////////////////////////////////////////////////////////////////////////////////
void create_truncated_cubes(LCC3& lcc, size_type corner_mark,
                            bool no_signature)
{
  std::chrono::system_clock::time_point start=std::chrono::system_clock::now();
  std::chrono::system_clock::time_point start2=std::chrono::system_clock::now();
  Pattern_substituer<LCC3> ps2; // pattern substituer to create transitions
  ps2.load_vpatterns(CGAL::data_file_path("query_replace/hexa-to-truncated-cube/vpattern/"),
                     mark_vpattern_corners<LCC3>);
  ps2.load_fpatterns(CGAL::data_file_path("query_replace/hexa-to-truncated-cube/fpattern/"),
                     mark_fpattern_corners<LCC3>);
  std::chrono::duration<double> diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Load transition patterns: ", diff, "s.");

  // 1) Split edges
  start2 = std::chrono::system_clock::now();
  std::vector<Dart_handle> to_treat;
  to_treat.reserve(lcc.number_of_attributes<3>()); // a la louche ;)
  for(auto itd=lcc.one_dart_per_cell<1>().begin(),
      itdend=lcc.one_dart_per_cell<1>().end(); itd!=itdend; ++itd)
  { to_treat.push_back(itd); }

  for(Dart_handle itd: to_treat)
  {
    Dart_handle itd2=lcc.beta<2>(itd);
    Vector v(lcc.point(itd), lcc.point(itd2));
    Point p1=lcc.point(itd);
    Point p2=lcc.point(itd2);
    p1=p1.transform(Kernel::Aff_transformation_3(CGAL::TRANSLATION, v/3.6));
    p2=p2.transform(Kernel::Aff_transformation_3(CGAL::TRANSLATION, -v/3.6));
    lcc.insert_point_in_cell<1>(itd, p1);
    lcc.insert_point_in_cell<1>(itd2, p2);
  }

  to_treat.clear();
  for(auto itd=lcc.one_dart_per_cell<3>().begin(),
      itdend=lcc.one_dart_per_cell<3>().end(); itd!=itdend; ++itd)
  { to_treat.push_back(itd); }
  
  std::size_t nb_replaced_f=0;
  std::size_t nb_replaced_v=0;
  
  std::vector<Dart_handle> faces;
  for(auto curdh: to_treat)
  {
    faces.clear();
    for(auto itf=lcc.one_dart_per_incident_cell<2,3>(curdh).begin(),
        itfend=lcc.one_dart_per_incident_cell<2,3>(curdh).end(); itf!=itfend; ++itf)
    { faces.push_back(itf); }
    if(no_signature)
    {
      for(auto itf: faces)
      {
        ps2.query_replace_one_face_without_signature(lcc, itf, corner_mark);
        ++nb_replaced_f;
      }
    }
    else
    {
      for(auto itf: faces)
      {
        ps2.query_replace_one_face(lcc, itf, corner_mark);
        ++nb_replaced_f;
      }
    }

    std::size_t res;
    if(no_signature)
    {
      res=ps2.query_replace_one_volume_without_signature(lcc, curdh, corner_mark);
      ++nb_replaced_v;
    }
    else
    {
      res=ps2.query_replace_one_volume(lcc, curdh, corner_mark);
      ++nb_replaced_v;
    }
    if(res==std::numeric_limits<std::size_t>::max())
    {
      std::cout<<"[ERROR] transition not found for volume "
                    <<lcc.attributes<3>().index(lcc.attribute<3>(curdh))
                   <<std::endl;
      static unsigned int nberror=0;
      LCC3 lcc_error;
      copy_cell<3>(lcc, curdh, lcc_error);
      CGAL::draw(lcc_error);
      save_object_3D(std::string("error"+std::to_string(nberror++)+
                                 ".mesh"), lcc_error);
    }
  }
  to_treat.clear();
  for(auto itd=lcc.one_dart_per_cell<2>().begin(),
      itdend=lcc.one_dart_per_cell<2>().end(); itd!=itdend; ++itd)
  {
    if(lcc.beta<1,1,1>(itd)==itd && !lcc.is_free<3>(itd) &&
       itd<lcc.beta<3>(itd) &&
       lcc.is_volume_combinatorial_tetrahedron(itd) &&
       lcc.is_volume_combinatorial_tetrahedron(lcc.beta<3>(itd)))
    { to_treat.push_back(itd);  }
  }
  for(auto itd: to_treat)
  { lcc.remove_cell<2>(itd); }

  diff=std::chrono::system_clock::now()-start;
  print_txt_with_endl("Create truncated cubes TOTAL: ", diff, "s.");

  std::cout<<"#faces-replaced: "<<nb_replaced_f<<", #volumes-replaced: "<<nb_replaced_v<<std::endl;
}
////////////////////////////////////////////////////////////////////////////////
bool create_mesh_from_off(const std::string& filename,
                          bool create_transitions,
                          bool draw,
                          unsigned int initX,
                          unsigned int lmax,
                          bool marching_cubes,
                          bool no_remove_outside,
                          bool no_signature,
                          bool save,
                          bool smooth,
                          bool truncated_cubes
                          )
{
  LCC3 lcc;
  CGAL::Polyhedron_3<Kernel> surface;
  Tree aabb_tree;
  size_type corner_mark=lcc.get_new_mark();

  if (!create_initial_voxels(filename, lcc, initX, surface, aabb_tree))
  {
    lcc.clear();
    return false;
  }
  lcc.negate_mark(corner_mark); // At the beginning, all darts are corners

  std::chrono::system_clock::time_point start=std::chrono::system_clock::now();
  std::chrono::system_clock::time_point start2=std::chrono::system_clock::now();

  Pattern_substituer<LCC3> ps1; // pattern substituer to subdivide hexa
  ps1.load_fpatterns(CGAL::data_file_path("query_replace/hexa-8-subdivision/fpattern/"),
                     mark_fpattern_corners<LCC3>);
  ps1.load_vpatterns(CGAL::data_file_path("query_replace/hexa-8-subdivision/vpattern/"),
                     mark_vpattern_corners<LCC3>);

  bool cont=true;
  unsigned int curlevel=0; //, level_of_hexa;
  std::vector<std::vector<Dart_handle>> hexa_tosubdivide(lmax);
  std::size_t nbtosubdivide;
  std::size_t nb_replaced_f=0;
  std::size_t nb_replaced_v=0;
    
  while (curlevel<lmax && cont)
  {
    cont=false;
    nbtosubdivide=(truncated_cubes?
                     compute_hexa_to_subdivide_for_truncated_cubes(lcc, aabb_tree,
                                                                   no_remove_outside,
                                                                   hexa_tosubdivide):
                     compute_hexa_to_subdivide(lcc, aabb_tree, curlevel,
                                               no_remove_outside, marching_cubes,
                                               hexa_tosubdivide));
    if(nbtosubdivide>0)
    {
      for(std::size_t i=0; i<=curlevel; ++i)
      {
        subdivide_edges(lcc, hexa_tosubdivide[i], corner_mark);
        subdivide_faces(lcc, hexa_tosubdivide[i], corner_mark, ps1, nb_replaced_f);
        subdivide_volumes(lcc, hexa_tosubdivide[i], corner_mark, ps1, nb_replaced_v);
      }
      cont=true;
    }
    ++curlevel;
  }

  std::chrono::duration<double> diff=std::chrono::system_clock::now()-start2;
  print_txt_with_endl("    Subdivision: ", diff, "s.");

  if(!no_remove_outside)
  {
    start2=std::chrono::system_clock::now();
    for (auto it=lcc.attributes<3>().begin(); it!=lcc.attributes<3>().end(); ++it)
    {
      if(lcc.info_of_attribute<3>(it)==lmax)
      {
        if((marching_cubes &&
            !is_intersect(lcc, lcc.dart_of_attribute<3>(it), aabb_tree)) ||
           (!marching_cubes &&
            is_outside(lcc, lcc.dart_of_attribute<3>(it), aabb_tree)))
        { lcc.remove_cell<3>(lcc.dart_of_attribute<3>(it)); }
      }
    }
    diff=std::chrono::system_clock::now()-start2;
    print_txt_with_endl("    Remove outside: ", diff, "s.");
  }

  if(create_transitions)
  { create_transition_patterns(lcc, corner_mark, no_signature); }
  else if(marching_cubes)
  { compute_marching_cubes(lcc, corner_mark, aabb_tree, no_signature); }

  if(smooth)
  { smooth_mesh(lcc, corner_mark, aabb_tree, no_signature); }

  if(truncated_cubes)
  { create_truncated_cubes(lcc, corner_mark, no_signature); }

  diff=std::chrono::system_clock::now()-start;
  print_txt_with_endl("Create mesh TOTAL: ", diff, "s.");
  std::cout<<"Final map: "; display_stats(lcc);
  std::cout<<"#faces-replaced: "<<nb_replaced_f<<", #volumes-replaced: "<<nb_replaced_v<<std::endl;
  assert(lcc.is_valid());

  if(draw)
  { CGAL::draw(lcc); }

  if(save)
  {
    std::filesystem::path p(filename);
    save_object_3D(p.stem().string()+std::to_string(lmax)+".mesh", lcc);
  }

  lcc.free_mark(corner_mark);
  return true;
}
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  if (argc<2) // We need at least one filename
  { usage(argc, argv); }

  std::string filename;
  bool create_transitions;
  bool draw;
  unsigned int initX;
  unsigned int lmax;
  bool marching_cubes;
  bool no_remove_outside;
  bool no_signature;
  bool save;
  bool smooth;
  bool truncated_cubes;

  process_command_line(argc, argv,
                       filename,
                       create_transitions,
                       draw,
                       initX,
                       lmax,
                       marching_cubes,
                       no_remove_outside,
                       no_signature,
                       save,
                       smooth,
                       truncated_cubes);

  if (!create_mesh_from_off(filename,
                            create_transitions,
                            draw,
                            initX,
                            lmax,
                            marching_cubes,
                            no_remove_outside,
                            no_signature,
                            save,
                            smooth,
                            truncated_cubes))
  { return EXIT_FAILURE; }

  return EXIT_SUCCESS;
}
////////////////////////////////////////////////////////////////////////////////
