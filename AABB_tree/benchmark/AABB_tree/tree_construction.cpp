#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>

int longest_axis(const CGAL::Bbox_3& bbox)
{
  const double dx = bbox.xmax() - bbox.xmin();
  const double dy = bbox.ymax() - bbox.ymin();
  const double dz = bbox.zmax() - bbox.zmin();
  return (dx>=dy) ? ((dx>=dz) ? 0 : 2) : ((dy>=dz) ? 1 : 2);
}

template <class RPM>
struct Split_primitives
{
  Split_primitives(RPM rpm)
    : rpm(rpm)
  {}

  template<typename PrimitiveIterator>
  void operator()(PrimitiveIterator first,
                  PrimitiveIterator beyond,
                  const CGAL::Bbox_3& bbox) const
    {
      PrimitiveIterator middle = first + (beyond - first)/2;
      typedef typename std::iterator_traits<PrimitiveIterator>::value_type Primitive;
      const int crd=longest_axis(bbox);
      const RPM& l_rpm=rpm;
      std::nth_element(first, middle, beyond,
                       [l_rpm, crd](const Primitive& p1, const Primitive& p2){ return get(l_rpm, p1.id())[crd] < get(l_rpm, p2.id())[crd];});
    }
   RPM rpm;
};

template <class BBM>
struct Compute_bbox {
  Compute_bbox(const BBM& bbm)
    : bbm(bbm)
  {}

  template<typename ConstPrimitiveIterator>
  CGAL::Bbox_3 operator()(ConstPrimitiveIterator first,
                          ConstPrimitiveIterator beyond) const
  {
    CGAL::Bbox_3 bbox = get(bbm, first->id());
    for(++first; first != beyond; ++first)
    {
      bbox += get(bbm, first->id());
    }
    return bbox;
  }
  BBM bbm;
};

template <class K>
void run(std::string input)
{
  typedef typename K::Point_3 Point_3;
  typedef CGAL::Surface_mesh<Point_3> Mesh;
  typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
  typedef CGAL::AABB_traits<K, Primitive> Traits;
  typedef CGAL::AABB_tree<Traits> Tree;

  Mesh tm;
  std::ifstream(input) >> tm;

  {
  Tree tree(faces(tm).begin(), faces(tm).end(), tm);
  CGAL::Timer time;
  time.start();
  tree.build();
  time.stop();
  std::cout << "  build() time: " << time.time() << "\n";
  }

  {
  Tree tree(faces(tm).begin(), faces(tm).end(), tm);
  CGAL::Timer time;
  time.start();

  typedef CGAL::Pointer_property_map<CGAL::Bbox_3>::type BBM;
  typedef CGAL::Pointer_property_map<CGAL::Epick::Point_3>::type RPM; // EPIC on purpose here

  std::vector<CGAL::Bbox_3> v_bb;
  std::vector<CGAL::Epick::Point_3> v_rp;

  std::size_t nbf = num_faces(tm);
  v_bb.resize(nbf);
  v_rp.resize(nbf);
  BBM bbm = CGAL::make_property_map(v_bb);
  RPM rpm = CGAL::make_property_map(v_rp);

  CGAL::Cartesian_converter<K, CGAL::Epick> to_input;
  for(typename Mesh::Face_index f : tm.faces())
  {
    v_bb[f]=CGAL::Polygon_mesh_processing::face_bbox(f, tm);
    v_rp[f]=to_input(tm.point(target(halfedge(f, tm), tm)));
  }

  Compute_bbox<BBM> compute_bbox(bbm);
  Split_primitives<RPM> split_primitives(rpm);
  tree.custom_build(compute_bbox, split_primitives);
  time.stop();
  std::cout << "  custom_build() time: " << time.time() << "\n";
  }
}

int main(int, char** argv)
{
  std::cout << "Build with Epick\n";
  run<CGAL::Epick>(argv[1]);
  std::cout << "Build with Epeck\n";
  run<CGAL::Epeck>(argv[1]);
  return EXIT_SUCCESS;
}
