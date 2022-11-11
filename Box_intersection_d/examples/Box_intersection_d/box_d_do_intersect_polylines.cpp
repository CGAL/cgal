#include <vector>
#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>


template <typename Poly>
class Polylines_do_intersect
{
  typedef typename Poly::const_iterator Iterator;
  typedef typename std::iterator_traits<Iterator>::value_type Point_2;
  typedef typename CGAL::Kernel_traits<Point_2>::Kernel::Segment_2 Segment_2;
  typedef CGAL::Bbox_2 Bbox_2;

  const Poly& pA, &pB;

  struct Box : public CGAL::Box_intersection_d::Box_with_handle_d<double,2,Iterator> {

    typedef CGAL::Box_intersection_d::Box_with_handle_d<double,2,Iterator> Base;

    Box( const Bbox_2& b, Iterator it)
      : Base(b,it)
    {}

  };


  class FirstIntersection
  {};


  class Overlap {

    const Poly& pA, &pB;
  public:
    Overlap(const Poly& pA, const Poly& pB)
      : pA(pA), pB(pB)
    {}


    void operator()(const Box& a, const Box& b) {
      Segment_2 sa(*(a.handle()), *(++(a.handle())));
      Segment_2 sb(*(b.handle()), *(++(b.handle())));
      if(CGAL::do_intersect(sa,sb)){
        throw FirstIntersection();
      }
    }
  };

public:

  Polylines_do_intersect(const Poly& pA, const Poly& pB)
    : pA(pA), pB(pB)
  {}

  bool operator()() const
  {
    std::vector<Box> bA, bB;
    bA.reserve(pA.size() - 1 );
    bB.reserve(pB.size() - 1 );

    Iterator begin = pA.begin();
    for(std::size_t j=0; j < pA.size()-1; j++){
      Bbox_2 bb = pA[j].bbox() + pA[j+1].bbox();
      bA.push_back(Box(bb, begin+j));
    }
    begin = pB.begin();
    for(std::size_t j=0; j < pB.size()-1; j++){
      Bbox_2 bb = pB[j].bbox() + pB[j+1].bbox();
      bB.push_back(Box(bb, begin+j));
    }

    Overlap overlap(pA,pB);

    try {
      CGAL::box_intersection_d( bA.begin(), bA.end(), bB.begin(), bB.end(), overlap);
    } catch(const FirstIntersection& ) {
      return true;
    }
    return false;
  }

};


template <typename Poly>
bool polylines_do_intersect(const Poly& pA,const Poly& pB)
{
  Polylines_do_intersect<Poly> pdi(pA,pB);
  return pdi();
}


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef std::vector<Point_2> CGAL_polyline;

int main()
{
  CGAL_polyline pA = { Point_2(0,0), Point_2(1,0) };
  CGAL_polyline pB = { Point_2(0,0), Point_2(0,1), Point_2(2,1) , Point_2(1,0)};

  if(polylines_do_intersect(pA,pB)){
    std::cout << "intersection" << std::endl;
  } else {
    std::cout << "no intersection" << std::endl;
  }
  return 0;
}
