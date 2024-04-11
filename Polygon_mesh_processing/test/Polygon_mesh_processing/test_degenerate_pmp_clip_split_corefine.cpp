#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// #define CGAL_DEBUG_PMP_CLIP

// TODO: test coref

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron;

template <class TriangleMesh>
void test()
{
  int i=0;
  auto run_a_test = [&i] (std::string f, double a, double b, double c, double d)
  {
    std::cout << "running test " << i << "\n";
    {
    std::cout << "   test clip with throw\n";
    TriangleMesh tm;
    std::ifstream(f) >> tm;
    try{
      PMP::clip(tm, K::Plane_3(a,b,c,d), params::throw_on_self_intersection(true));
    }
    catch(PMP::Corefinement::Self_intersection_exception&)
    {}
    }

    {
    std::cout << "   test split with SI allowed\n";
    TriangleMesh tm;
    std::ifstream(f) >> tm;
    PMP::split(tm, K::Plane_3(a,b,c,d), params::allow_self_intersections(true));
    }

    {
    std::cout << "   test corefine with SI allowed\n";
    TriangleMesh tm;
    std::ifstream(f) >> tm;
    TriangleMesh box_mesh;
    CGAL::Bbox_3 bbox = PMP::bbox(tm);

    //extend the bbox a bit to avoid border cases
    double xd=(std::max)(1.,(bbox.xmax()-bbox.xmin())/100);
    double yd=(std::max)(1.,(bbox.ymax()-bbox.ymin())/100);
    double zd=(std::max)(1.,(bbox.zmax()-bbox.zmin())/100);
    bbox=CGAL::Bbox_3(bbox.xmin()-xd, bbox.ymin()-yd, bbox.zmin()-zd,
                      bbox.xmax()+xd, bbox.ymax()+yd, bbox.zmax()+zd);
    PMP::internal::clip_to_bbox(K::Plane_3(a,b,c,d), bbox,
                                box_mesh, params::default_values());
    try{
      PMP::corefine(tm, box_mesh, params::throw_on_self_intersection(true));
    }
    catch(PMP::Corefinement::Self_intersection_exception&)
    {}

    PMP::corefine(tm, box_mesh, params::default_values(), params::do_not_modify(true));

    assert(CGAL::is_triangle_mesh(tm));
    assert(CGAL::is_triangle_mesh(box_mesh));
    }

    {
    std::cout << "   test clip with SI allowed\n";
    TriangleMesh tm;
    std::ifstream(f) >> tm;
#ifdef CGAL_DEBUG_PMP_CLIP
    std::ofstream("/tmp/input_"+std::to_string(i)+".off") << std::setprecision(17) << tm;
#endif
    PMP::clip(tm, K::Plane_3(a,b,c,d), params::allow_self_intersections(true));
#ifdef CGAL_DEBUG_PMP_CLIP
    std::ofstream("/tmp/output_"+std::to_string(i)+".off") << std::setprecision(17) << tm;
#endif
    }

    ++i;
  };

  run_a_test("data_degeneracies/deg_on_border.off",0.995127, 0.0129458, -0.0977449, -0.501428);
  run_a_test("data_degeneracies/deg_on_border.off",0.999671, 0.0211974, -0.0144363, -2.32495);
  run_a_test("data_degeneracies/degtri_2dt_1edge_split_twice.off",0.995127, 0.0129458, -0.0977449, -1.34958);
  run_a_test("data_degeneracies/degtri_edge.off",0.995127, 0.0129458, -0.0977449,-0.501428);
  run_a_test("data_degeneracies/degtri_four-2.off",1, 1.67292e-19, -2.22045e-16, -0.618478);
  run_a_test("data_degeneracies/degtri_four-2.off",1, 0, 0, -1);
  run_a_test("data_degeneracies/degtri_four.off",1, 1.67292e-19, -2.22045e-16, -0.618478);
  run_a_test("data_degeneracies/degtri_four.off",1, 0, 0, -1);
  run_a_test("data_degeneracies/degtri_on_border.off",1, 1.67292e-19, -2.22045e-16, -0.618478);
  run_a_test("data_degeneracies/degtri_on_border.off",1, 0, 0, -1);
  run_a_test("data_degeneracies/degtri_single.off",1, 1.67292e-19, -2.22045e-16, -0.618478);
  run_a_test("data_degeneracies/degtri_single.off",1, 0, 0, -1);
  run_a_test("data_degeneracies/degtri_three.off",1, 1.67292e-19, -2.22045e-16, -0.618478);
  run_a_test("data_degeneracies/degtri_three.off",1, 0, 0, -1);
  run_a_test("data_degeneracies/trihole.off",1, 1.67292e-19, -2.22045e-16, -0.618478);
  run_a_test("data_degeneracies/trihole.off",1, 0, 0, -1);
  run_a_test("data_degeneracies/degtri_nullface.off",1, 1.67292e-19, -2.22045e-16, -0.618478);
  run_a_test("data_degeneracies/degtri_nullface.off",1, 0, 0, -1);
  run_a_test("data_degeneracies/degtri_nullface_bis.off",1, 0, 0, -1);
  // TODO subdivide several times the central face
}

int main()
{
  test<Surface_mesh>();
//  test<Polyhedron>();
}
