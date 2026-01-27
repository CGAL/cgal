#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Octree.h>
#include <CGAL/Orthtree_traits_face_graph.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;

using OTraits = CGAL::Orthtree_traits_face_graph<Mesh, Mesh::Property_map<Mesh::Vertex_index, K::Point_3>>;
using Octree = CGAL::Orthtree<OTraits>;

void dump_as_polylines(const Octree& ot)
{
    std::ofstream out("octree.polylines.txt");
    for (Octree::Node_index node : ot.traverse(CGAL::Orthtrees::Leaves_traversal<Octree>(ot)))
    {
      if (!ot.is_leaf(node))
        continue;
      auto bb = ot.bbox(node);
      out << "2 " << bb.xmin() << " " << bb.ymin() << " " << bb.zmin()
          << "  " << bb.xmax() << " " << bb.ymin() << " " << bb.zmin() << "\n";
      out << "2 " << bb.xmin() << " " << bb.ymin() << " " << bb.zmin()
          << "  " << bb.xmin() << " " << bb.ymax() << " " << bb.zmin() << "\n";
      out << "2 " << bb.xmax() << " " << bb.ymin() << " " << bb.zmin()
          << "  " << bb.xmax() << " " << bb.ymax() << " " << bb.zmin() << "\n";
      out << "2 " << bb.xmin() << " " << bb.ymax() << " " << bb.zmin()
          << "  " << bb.xmax() << " " << bb.ymax() << " " << bb.zmin() << "\n";
//
      out << "2 " << bb.xmin() << " " << bb.ymin() << " " << bb.zmin()
          << "  " << bb.xmin() << " " << bb.ymin() << " " << bb.zmax() << "\n";
      out << "2 " << bb.xmax() << " " << bb.ymin() << " " << bb.zmin()
          << "  " << bb.xmax() << " " << bb.ymin() << " " << bb.zmax() << "\n";
      out << "2 " << bb.xmin() << " " << bb.ymax() << " " << bb.zmin()
          << "  " << bb.xmin() << " " << bb.ymax() << " " << bb.zmax() << "\n";
      out << "2 " << bb.xmax() << " " << bb.ymax() << " " << bb.zmin()
          << "  " << bb.xmax() << " " << bb.ymax() << " " << bb.zmax() << "\n";
//
      out << "2 " << bb.xmin() << " " << bb.ymin() << " " << bb.zmax()
          << "  " << bb.xmax() << " " << bb.ymin() << " " << bb.zmax() << "\n";
      out << "2 " << bb.xmin() << " " << bb.ymin() << " " << bb.zmax()
          << "  " << bb.xmin() << " " << bb.ymax() << " " << bb.zmax() << "\n";
      out << "2 " << bb.xmax() << " " << bb.ymin() << " " << bb.zmax()
          << "  " << bb.xmax() << " " << bb.ymax() << " " << bb.zmax() << "\n";
      out << "2 " << bb.xmin() << " " << bb.ymax() << " " << bb.zmax()
          << "  " << bb.xmax() << " " << bb.ymax() << " " << bb.zmax() << "\n";
    }
}

int main(int argc, char** argv)
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/elephant.off");

  Mesh mesh;
  if(!CGAL::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Error: cannot read file" << std::endl;
    return EXIT_FAILURE;
  }

  Octree tree(mesh, mesh.points());
  OTraits::Split_predicate_node_min_extent sp(0.01);
  tree.refine(sp);

  dump_as_polylines(tree);
}
