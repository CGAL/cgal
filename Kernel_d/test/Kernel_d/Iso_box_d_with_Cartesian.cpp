#include <CGAL/Cartesian_d.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>

typedef CGAL::Cartesian_d<double> Kernel;
typedef Kernel::Point_d           Point;
typedef Kernel::Iso_box_d         Box;
typedef Kernel::Construct_iso_box_d Construct_iso_box_d;

int main() {

  const int dim=4;

  // define box b1==*b2==b3, b4
  // define degenerate box b5 of dimension dim+1

  double p1[dim];
  double q1[dim];
  double m1[dim];
  for (int i1=0; i1<dim; i1++) {
          p1[i1]=  0.1*(i1+1);
        q1[i1]=  0.2*(i1+1);
        m1[i1]= 0.15*(i1+1);
  }

  Point pp1(dim,p1,p1+dim);
  Point qq1(dim,q1,q1+dim);
  Point mm1(dim,m1,m1+dim);
  Box b1(pp1,qq1);

  double p2[dim];
  double q2[dim];
  for (int i2=0; i2<dim; i2++) {
          p2[i2]=  0.1*(i2+1);
        q2[i2]=  0.2*(i2+1);
  }
  Point pp2(dim,p2,p2+dim);
  Point qq2(dim,q2,q2+dim);
  Box *b2 = new Box(pp2,qq2);

  Box b3(*b2);

  delete b2;

  double p4[dim];
  double q4[dim];
  for (int i4=0; i4<dim; i4++) {
          p4[i4]=  1.0 + 0.1*(i4+1);
        q4[i4]=  1.0 + 0.2*(i4+1);
  }
  Point pp4(dim,p4,p4+dim);
  Point qq4(dim,q4,q4+dim);
  Box b4(pp4,qq4);

  double p5[dim+1];
  double q5[dim+1];
  for (int i5=0; i5<dim; i5++) {
          p5[i5]=  1.0 + 0.1*(i5+1);
        q5[i5]=  1.0 + 0.2*(i5+1);
  }
  p5[dim]=1.0;
  q5[dim]=1.0;
  Point pp5(dim+1,p5,p5+dim+1);
  Point qq5(dim+1,q5,q5+dim+1);
  Box b5(pp5,qq5);

  Construct_iso_box_d construct_iso_box_d = Kernel().construct_iso_box_d_object();
  Box b6 = construct_iso_box_d(pp5, qq5);
  assert(b5 == b6);
  assert(b1==b3);
  assert(b1!=b4);
  assert(b3!=b4);
  assert( !(b1.is_degenerate()));
  assert( b5.is_degenerate());
  assert(b1.dimension()==4);
  std::cout << "b1.min_coord(0)=" << b1.min_coord(0) << std::endl;
  std::cout << "b1.max_coord(0)=" << b1.max_coord(0) << std::endl;
  std::cout << "b1.min_coord(1)=" << b1.min_coord(1) << std::endl;
  std::cout << "b1.max_coord(1)=" << b1.max_coord(1) << std::endl;
  std::cout << "b1.min_coord(2)=" << b1.min_coord(2) << std::endl;
  std::cout << "b1.max_coord(2)=" << b1.max_coord(2) << std::endl;
  std::cout << "b1.min_coord(3)=" << b1.min_coord(3) << std::endl;
  std::cout << "b1.max_coord(3)=" << b1.max_coord(3) << std::endl;
  assert(b1.has_on_boundary(pp1));
  assert(b1.has_on_boundary(qq1));
  assert(b1.has_on_unbounded_side(qq4));
  assert(b1.has_on_bounded_side(mm1));
  assert(b1.bounded_side(pp1)==CGAL::ON_BOUNDARY);
  assert(b1.bounded_side(qq1)==CGAL::ON_BOUNDARY);
  assert(b1.bounded_side(qq4)==CGAL::ON_UNBOUNDED_SIDE);
  assert(b1.bounded_side(mm1)==CGAL::ON_BOUNDED_SIDE);
  // std::cout << "(b1.min)()=" << (b1.min)() << std::endl;
  // std::cout << "(b1.max)()=" << (b1.max)() << std::endl;
  std::cout << "volume of b1=" << b1.volume() << std::endl;
  std::cout << "volume of b4=" << b4.volume() << std::endl;
  std::cout << "volume of b5=" << b5.volume() << std::endl;
  return 0;
}

