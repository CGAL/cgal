#include <iostream>

namespace std {
template <typename T>
struct HashFct {
  void operator()(const T& t)
  {
    std::cerr << "generic HashFct" << std::endl;
  }
};
}

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

namespace std {
  template <>
struct HashFct<CGAL::SM_Halfedge_index > {
  void operator()(const CGAL::SM_Halfedge_index& t)
  {
    std::cerr << "SM_Halfegge_index HashFct" << std::endl;
  }
};
}

template <typename T>
struct Unordered
{
  Unordered()
  {
    std::HashFct<T> fct;
    T t;
    fct(t);
  }
};

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point;
typedef CGAL::Surface_mesh<Point>       Mesh;
typedef Mesh::Halfedge_index            Halfedge_index; 

int main()
{
  Unordered<Halfedge_index> U;
  return 0;
}
