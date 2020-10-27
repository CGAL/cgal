
#include <CGAL/Orthtree/Node.h>
#include <CGAL/Octree/IO.h>
#include <iostream>
#include <cassert>

typedef CGAL::Orthtree::Node<std::vector<int>::iterator> Node;

int main(void) {

  // Build a new node
  Node n = Node();

  // Check that its values are correct
  assert(n.is_root());
  assert(n.is_leaf());
  assert(!n.parent());
  assert(n.depth() == 0);
  assert(n.location()[0] == 0 && n.location()[1] == 0 && n.location()[2] == 0);

  // Split the node
  n.split();

  // Check that it's children's values are also correct
  for (std::size_t i = 0; i < 8; ++i) {

    assert(!n[i].is_root());
    assert(n[i].is_leaf());
    assert(*n[i].parent() == n);
    assert(n[i].depth() == 1);
  }

  // Check that the parent has updated
  assert(n.is_root());
  assert(!n.is_leaf());

  // Split one of the children
  n[1].split();

  // Check each of that child's children
  for (std::size_t i = 0; i < 8; ++i) {

    assert(!n[1][i].is_root());
    assert(n[1][i].is_leaf());
    assert(*n[1][i].parent() == n[1]);
    assert(*n[1][i].parent()->parent() == n);
    assert(n[1][i].depth() == 2);
  }

  // Check that the child's values have updated
  assert(!n[1].is_root());
  assert(!n[1].is_leaf());

  return 0;
}
