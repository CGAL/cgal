#include <iostream>
#include <cassert>

#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>

#include "test_three_sites.h"

typedef CGAL::Exact_rational exact_nt;

typedef CGAL::Simple_cartesian<exact_nt> Rep;

typedef CGAL::Integral_domain_without_division_tag Method_tag;

typedef CGAL::Apollonius_graph_traits_2<Rep,Method_tag> Gt;

typedef CGAL::Apollonius_graph_2<Gt>  AG2;
typedef AG2::Point_2                  Point_2;
typedef AG2::Site_2                   Site_2;
typedef AG2::Vertex_handle            Vertex_handle;

int main()
{
  std::cout << "testing the Apollonius graph class for three sites..."
            << std::endl;

  Site_2 s1(Point_2(100,100),5), s2(Point_2(150,150),10),
    s3(Point_2(200,200),15);

  AG2 ag;

  {
    ag.clear();
    Vertex_handle v1 = ag.insert(s1);
    Vertex_handle v2 = ag.insert(s2);
    Vertex_handle v3 = ag.insert(s3);
    bool b = test_three_sites(ag, v1, v2, v3, MIDDLE_ON_CONVEX_HULL);
    std::cout << "Is AG okay? " << (b ? "YES" : "NO") << std::endl;
    std::cout << std::endl << "-------------------" << std::endl << std::endl;
    assert( b );
  }

  {
    ag.clear();
    Vertex_handle v1 = ag.insert(s1);
    Vertex_handle v3 = ag.insert(s3);
    Vertex_handle v2 = ag.insert(s2);
    bool b = test_three_sites(ag, v1, v2, v3, MIDDLE_ON_CONVEX_HULL);
    std::cout << "Is AG okay? " << (b ? "YES" : "NO") << std::endl;
    std::cout << std::endl << "-------------------" << std::endl << std::endl;
    assert( b );
  }

  Site_2 t1(Point_2(100,100),20), t2(Point_2(150,150),20),
    t3(Point_2(200,200),20);

  {
    ag.clear();
    Vertex_handle v1 = ag.insert(t1);
    Vertex_handle v2 = ag.insert(t2);
    Vertex_handle v3 = ag.insert(t3);
    bool b = test_three_sites(ag, v1, v2, v3, MIDDLE_ON_CONVEX_HULL);
    std::cout << "Is AG okay? " << (b ? "YES" : "NO") << std::endl;
    std::cout << std::endl << "-------------------" << std::endl << std::endl;
    assert( b );
  }

  {
    ag.clear();
    Vertex_handle v1 = ag.insert(t1);
    Vertex_handle v3 = ag.insert(t3);
    Vertex_handle v2 = ag.insert(t2);
    bool b = test_three_sites(ag, v1, v2, v3, MIDDLE_ON_CONVEX_HULL);
    std::cout << "Is AG okay? " << (b ? "YES" : "NO") << std::endl;
    std::cout << std::endl << "-------------------" << std::endl << std::endl;
    assert( b );
  }

  {
    ag.clear();
    Site_2 q(Point_2(150,150),25);
    Vertex_handle v1 = ag.insert(t1);
    Vertex_handle v2 = ag.insert(q);
    Vertex_handle v3 = ag.insert(t3);
    bool b = test_three_sites(ag, v1, v2, v3, MIDDLE_ON_CONVEX_HULL);
    std::cout << "Is AG okay? " << (b ? "YES" : "NO") << std::endl;
    std::cout << std::endl << "-------------------" << std::endl << std::endl;
    assert( b );
  }

  {
    ag.clear();
    Site_2 q(Point_2(150,150),25);
    Vertex_handle v1 = ag.insert(t1);
    Vertex_handle v3 = ag.insert(t3);
    Vertex_handle v2 = ag.insert(q);
    bool b = test_three_sites(ag, v1, v2, v3, MIDDLE_ON_CONVEX_HULL);
    std::cout << "Is AG okay? " << (b ? "YES" : "NO") << std::endl;
    std::cout << std::endl << "-------------------" << std::endl << std::endl;
    assert( b );
  }

  {
    ag.clear();
    Site_2 q(Point_2(150,150),15);
    Vertex_handle v1 = ag.insert(t1);
    Vertex_handle v2 = ag.insert(q);
    Vertex_handle v3 = ag.insert(t3);
    bool b = test_three_sites(ag, v1, v2, v3, MIDDLE_NOT_ON_CONVEX_HULL);
    std::cout << "Is AG okay? " << (b ? "YES" : "NO") << std::endl;
    std::cout << std::endl << "-------------------" << std::endl << std::endl;
    assert( b );
  }

  {
    ag.clear();
    Site_2 q(Point_2(150,150),15);
    Vertex_handle v1 = ag.insert(t1);
    Vertex_handle v3 = ag.insert(t3);
    Vertex_handle v2 = ag.insert(q);
    bool b = test_three_sites(ag, v1, v2, v3, MIDDLE_NOT_ON_CONVEX_HULL);
    std::cout << "Is AG okay? " << (b ? "YES" : "NO") << std::endl;
    std::cout << std::endl << "-------------------" << std::endl << std::endl;
    assert( b );
  }

  std::cout << "... testing is done!" << std::endl;
  return 0;
}
