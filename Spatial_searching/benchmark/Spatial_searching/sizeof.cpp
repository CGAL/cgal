#include <iostream>
#include <boost/cstdint.hpp>


struct Node {
  bool b;
};


struct Plane_separator{ 
  unsigned int cutting_dim;
  double cutting_val;
  
};


  struct InternalA : public Node {
    Plane_separator ps;
    double d;
  };

  struct InternalB : public Node {
  unsigned int cutting_dim;
  double cutting_val;
    double d;
  };

  struct InternalC : public Node {
    //Plane_separator ps;
  double cutting_val;
  int cutting_dim;
    double d;
  };

  struct Leaf : public Node {
    unsigned int i;
    double d;
  };




int main()
{
  Node n;
  Plane_separator ps;
  InternalA ia;
  InternalB ib;
  InternalC ic;
  Leaf l;
  std::cerr << sizeof(int) << std::endl;
  std::cerr << sizeof(double) << std::endl;
  std::cerr << sizeof(n) << std::endl;
  std::cerr << sizeof(ps) << std::endl;
  std::cerr << sizeof(ia) << std::endl;
  std::cerr << sizeof(ib) << std::endl;
  std::cerr << sizeof(ic) << std::endl;
  std::cerr << sizeof(l) << std::endl;

  Leaf  leafs[10];
  std::cerr << sizeof(leafs) << std::endl;
  
  Node nodes[10];
  std::cerr << sizeof(nodes) << std::endl;
  return 0;
}
