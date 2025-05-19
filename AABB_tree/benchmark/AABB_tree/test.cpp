#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/transform.h>

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/Rigid_triangle_mesh_collision_detection.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Surface_mesh;

typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> Primitive;
typedef CGAL::AABB_traits_3<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

namespace PMP = CGAL::Polygon_mesh_processing;

void naive_test(int k, const std::string& fname,
                int& nb_inter, int& nb_no_inter, int& nb_include)
{
  std::ifstream input(fname);
  Surface_mesh tm, tm2;
  input >> tm;
  copy_face_graph(tm, tm2);
  CGAL::Aff_transformation_3<K> init1(CGAL::SCALING, 6.0);
  PMP::transform(init1, tm);
  CGAL::Bbox_3 box = PMP::bbox(tm);

  Tree tmTree(tm.faces_begin(), tm.faces_end(), tm);
  Tree tmTree2(tm2.faces_begin(), tm2.faces_end(), tm2);
  CGAL::Aff_transformation_3<K> init2(CGAL::TRANSLATION, - K::Vector_3(
                                        (box.xmax()-box.xmin()),0,0));
  PMP::transform(init2, tm2);

  tmTree.build();
  K::Vector_3 unit_vec = (2.0/k * K::Vector_3((box.xmax()-box.xmin()),
                                              0,
                                              0));
  CGAL::Aff_transformation_3<K> T0(CGAL::IDENTITY);
  K::FT rot[9];
  rot[0] = 1.0;
  rot[1] = 0.0;
  rot[2] = 0.0;
  rot[3] = 0.0;
  rot[4] = std::cos(CGAL_PI/4.0);
  rot[5] = -std::sin(CGAL_PI/4.0);
  rot[6] = 0.0;
  rot[7] = std::sin(CGAL_PI/4.0);
  rot[8] = std::cos(CGAL_PI/4.0);
  CGAL::Aff_transformation_3<K> R(rot[0], rot[1], rot[2],
      rot[3], rot[4], rot[5],
      rot[6], rot[7], rot[8]);

  CGAL::Side_of_triangle_mesh<Surface_mesh, K> sotm1(tm);
  for(int i=1; i<k+1; ++i)
  {
    CGAL::Aff_transformation_3<K> T1 = CGAL::Aff_transformation_3<K>(CGAL::TRANSLATION, i*unit_vec);
    CGAL::Aff_transformation_3<K> transfo = T0*R*T1;
    PMP::transform(transfo, tm2);
    tmTree2.build();
    if(tmTree2.do_intersect(tmTree))
      ++nb_inter;
    else
    {
      if(sotm1(tm2.point(*tm2.vertices().begin())) != CGAL::ON_UNBOUNDED_SIDE)
      {
        ++nb_include;
      }
      else
      {
        CGAL::Side_of_triangle_mesh<Surface_mesh, K> sotm2(tm2);
        if(sotm2(tm.point(*tm.vertices().begin())) != CGAL::ON_UNBOUNDED_SIDE)
          ++nb_include;
        else
          ++nb_no_inter;
      }
    }
    T0 = CGAL::Aff_transformation_3<K>(CGAL::TRANSLATION, -i*unit_vec);
  }
}

void test_no_collision(int k, const std::string &fname,
                       int& nb_inter, int& nb_no_inter, int& nb_include)
{
  std::ifstream input(fname);
  Surface_mesh tm, tm2;
  input >> tm;
  copy_face_graph(tm, tm2);
  CGAL::Aff_transformation_3<K> init1(CGAL::SCALING, 6.0);
  PMP::transform(init1, tm);
  CGAL::Bbox_3 box = PMP::bbox(tm);
  Tree tmTree(tm.faces_begin(), tm.faces_end(), tm);

  Tree tmTree2(tm2.faces_begin(), tm2.faces_end(), tm2);
  CGAL::Aff_transformation_3<K> init2(CGAL::TRANSLATION, - K::Vector_3(
                                        (box.xmax()-box.xmin()),0,0));

  PMP::transform(init2, tm2);

  tmTree.build();
  tmTree2.build();
  typedef boost::property_map<Surface_mesh, CGAL::vertex_point_t>::type VPM;
  VPM vpm2 = get(CGAL::vertex_point, tm2);

  K::Vector_3 unit_vec = (2.0/k * K::Vector_3((box.xmax()-box.xmin()),
                                              0,
                                              0));

  CGAL::Side_of_triangle_mesh<Surface_mesh, K,
      VPM, Tree> sotm1(tmTree);

  CGAL::Rigid_triangle_mesh_collision_detection<Surface_mesh> collision_detection;

  collision_detection.add_mesh(tm);
  collision_detection.add_mesh(tm2);

  for(int i=1; i<k+1; ++i)
  {
    K::FT rot[9];
    rot[0] = 1.0;
    rot[1] = 0.0;
    rot[2] = 0.0;
    rot[3] = 0.0;
    rot[4] = std::cos(i*CGAL_PI/4.0);
    rot[5] = -std::sin(i*CGAL_PI/4.0);
    rot[6] = 0.0;
    rot[7] = std::sin(i*CGAL_PI/4.0);
    rot[8] = std::cos(i*CGAL_PI/4.0);
    CGAL::Aff_transformation_3<K> R(rot[0], rot[1], rot[2],
        rot[3], rot[4], rot[5],
        rot[6], rot[7], rot[8]);
    CGAL::Aff_transformation_3<K> T1 = CGAL::Aff_transformation_3<K>(CGAL::TRANSLATION, i*unit_vec);
    CGAL::Aff_transformation_3<K> transfo = R*T1;

    collision_detection.set_transformation(1, transfo);

    std::vector< std::pair<std::size_t, bool> > res = collision_detection.get_all_intersections_and_inclusions(0);

    if (res.empty())
      nb_no_inter++;
    else if(!res[0].second)
      ++nb_inter;
    else
      ++nb_include;
  }
}

int main(int argc, const char** argv)
{
  int k = (argc>1) ? atoi(argv[1]) : 20;
  std::string path = (argc>2)?argv[2]: CGAL::data_file_path("meshes/handle.off");

  std::cout<< k<<" steps in "<<path<<std::endl;
  int nb_inter(0), nb_no_inter(0), nb_include(0),
      naive_inter(0), naive_no_inter(0), naive_include(0);
  auto start = std::chrono::steady_clock::now();
  naive_test(k, path, naive_inter, naive_no_inter, naive_include);
  auto end = std::chrono::steady_clock::now();
  std::cout<<"Naive test           : "<<naive_inter<<" collisions, "<<naive_include<<" inclusions, "<<naive_no_inter<<" no collision, calculated in "
          <<std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms." << std::endl;
  start = std::chrono::steady_clock::now();
  test_no_collision(k, path,nb_inter, nb_no_inter, nb_include);
  end = std::chrono::steady_clock::now();
  std::cout<<"With transform_traits: "<<nb_inter<<" collisions, "<<nb_include<<" inclusions, "<<nb_no_inter<<" no collision, calculated in "
          <<std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms." << std::endl;
    return 0;
}
