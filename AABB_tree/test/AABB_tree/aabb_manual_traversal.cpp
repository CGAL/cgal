#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/boost/graph/generators.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using Mesh = CGAL::Surface_mesh<Point>;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using Traits = CGAL::AABB_traits_3<K, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;

struct Node_property
{

};

int main()
{
  Point p(1.0, 0.0, 0.0);
  Point q(0.0, 1.0, 0.0);
  Point r(0.0, 0.0, 1.0);
  Point s(0.0, 0.0, 0.0);
  Mesh sm;
  CGAL::make_tetrahedron(p, q, r, s, sm);
  CGAL::make_tetrahedron(p, q, r, s, sm);

  Tree tree( faces(sm).first, faces(sm).second, sm);

  using Node = CGAL::AABB_node<Traits>;
  const Node* root = tree.root_node();

  auto node_id = [root](const Node* node) -> std::size_t
  {
    return std::size_t(node - root);
  };

  std::size_t nb_primitives = tree.size();
  std::vector<Node_property> node_properties(nb_primitives-1);



  std::vector< std::tuple<const Node*, Tree::Primitive_iterator, std::size_t> > traversal_queue;
  traversal_queue.emplace_back(root, tree.primitives_begin(), nb_primitives);

  while (!traversal_queue.empty())
  {
    auto [node, prim_it, nb_prim] = traversal_queue.back();
    traversal_queue.pop_back();

    std::cout << node_id(node) << "\n";
    CGAL::Bbox_3 bbox = node->bbox();
    CGAL_USE(bbox);
    Node_property& prop = node_properties[node_id(node)];
    CGAL_USE(prop);

    for (auto it = prim_it; it<prim_it+nb_prim; ++it)
    {
      Mesh::Face_index f = it->id();
      std::cout << "   " << f << "\n";
    }

    switch(nb_prim)
    {
    case 2:
    {
      // leaf node
      Mesh::Face_index fl = node->left_data().id();
      Mesh::Face_index fr = node->right_data().id();
      CGAL_USE(fl);
      CGAL_USE(fr);
    }
    break;
    case 3:
    {
      // partial leaf node
      // left contains data, right is a node with the two remaining primitives
      Mesh::Face_index fl = node->left_data().id();
      Mesh::Face_index frl = node->right_child().left_data().id();
      Mesh::Face_index frr = node->right_child().right_data().id();
      CGAL_USE(fl);
      CGAL_USE(frl);
      CGAL_USE(frr);
    }
    break;
    default:
    {
      std::size_t nb_left=nb_prim/2;
      std::size_t nb_right=nb_prim - nb_left;

      traversal_queue.emplace_back(std::addressof(node->left_child()), prim_it, nb_left);
      traversal_queue.emplace_back(std::addressof(node->right_child()), prim_it+nb_left, nb_right);
    }
    }
  }

  return EXIT_SUCCESS;
}
