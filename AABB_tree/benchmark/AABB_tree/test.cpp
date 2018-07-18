#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/transform.h>

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_do_intersect_transform_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>

#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>             Surface_mesh;

typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_do_intersect_transform_traits<AABB_triangle_traits, K> Traits;
typedef CGAL::AABB_tree<Traits> Tree;

namespace PMP = CGAL::Polygon_mesh_processing;

void test_no_collision(int k, const char* fname,
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
  
  std::ofstream out1("/home/gimeno/Data/tmp/m1.off"),
      out2("/home/gimeno/Data/tmp/m2.off");
  out1 << tm;
  out1.close();
  out2 << tm2;
  out2.close();
  tmTree.build();
  tmTree2.build();
  typedef boost::property_map<Surface_mesh, CGAL::vertex_point_t>::type VPM;
  //VPM vpm = get(CGAL::vertex_point, tm);
  
  K::Vector_3 unit_vec = (2.0/k * K::Vector_3((box.xmax()-box.xmin()),
                                              0,
                                              0));
  CGAL::Side_of_triangle_mesh<Surface_mesh, K,
      VPM, Tree> sotm1(tmTree);
  for(int i=1; i<k+1; ++i)
  {
    CGAL::Aff_transformation_3<K> trans22(CGAL::TRANSLATION, i*unit_vec);
    box = PMP::bbox(tm2);
    tmTree2.traits().set_transformation(trans22);
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
        CGAL::Side_of_triangle_mesh<Surface_mesh, K,
            VPM, Tree> sotm2(tmTree2);
        if(sotm2(tm2.point(*tm2.vertices().begin())) != CGAL::ON_UNBOUNDED_SIDE)
          ++nb_include; 
        else
          ++nb_no_inter;
      }
    }
  }
}

int main(int argc, const char** argv)
{
  int k = (argc>1) ? atoi(argv[1]) : 10;
  const char* path = (argc>2)?argv[2]:"data/handle"
                                      ".off";
  std::cout<< k<<" steps in "<<path<<std::endl;
  int nb_inter, nb_no_inter, nb_include;
  auto start = std::chrono::steady_clock::now();
  test_no_collision(k, path,nb_inter, nb_no_inter, nb_include);
  auto end = std::chrono::steady_clock::now();
  
  std::cout<<nb_inter<<" collisions, "<<nb_include<<" inclusions, "<<nb_no_inter<<" no collision, calculated in "
          <<std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "μs." << std::endl;
}
