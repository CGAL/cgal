#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Bbox_d.h>

#include <array>

int main()
{
  //Dimension d
  typedef CGAL::Bbox_d<CGAL::Dimension_tag<2>> BBox2;
  typedef CGAL::Bbox_d<CGAL::Dimension_tag<3>> BBox3;
  BBox3 bb3, bb3a(1.0);
  assert(bb3.dimension() == 3);
  assert(bb3 != bb3a);
  bb3 = bb3a;
  assert(bb3 == bb3a);

  std::array<std::pair<double,double>,3> coord = { std::make_pair(0.0, 0.0), std::make_pair(1.0, 1.1), std::make_pair(1.0, 20.0)};

  BBox3 bb3b(coord.begin(), coord.end());

  bb3b = bb3b + bb3b;
  bb3b += bb3;
  bb3b.dilate(15);
  assert(CGAL::do_overlap(bb3, bb3b));

  std::cout <<  bb3b << std::endl;

  BBox3::Cartesian_const_iterator beg = bb3b.cartesian_begin();
  BBox3::Cartesian_const_iterator end =  bb3b.cartesian_end();
  for(; beg != end; ++beg){
    std::cout << *beg << std::endl;
  }

  CGAL::Bbox_2 bb_2(0,0, 1, 1);
  BBox2 bb_d2(bb_2);

  CGAL::Bbox_3 bb_3(0,0, 0, 1, 1,1);
  BBox3 bb_d3(bb_3);

}
