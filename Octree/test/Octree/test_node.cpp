
#include <CGAL/Octree/Node.h>
#include <CGAL/Octree/IO.h>
#include <iostream>

typedef CGAL::Octree::Node::Node<int> Node;

int main(void) {

  Node n{};

  n.split();

  std::cout << n;
  for (int i = 0; i < 8; ++i) {
    std::cout << n[i];
  }


  return 0;
}
