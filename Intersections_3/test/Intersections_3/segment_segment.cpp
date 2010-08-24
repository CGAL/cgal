#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <cassert>

typedef CGAL::Quotient< CGAL::MP_Float >                    FT;
typedef CGAL::Cartesian<FT>                                 K;

// The construction test has to be working
// The equal test has to be working

template <class K>
void _test_intersection_construct(K)
{
  typedef typename K::FT                               FT;
  typedef typename K::Point_3                          Point_3;
  typedef typename K::Segment_3                        Segment_3;

  std::cout << "Testing intersection(Segment, Segment)..." << std::endl;

  Segment_3 s[5]={  Segment_3( Point_3( 0,1,0),Point_3(0,-1,0) ),
                    Segment_3( Point_3(-1,0,0),Point_3(1, 0,0) ),
                    Segment_3( Point_3(-2,0,0),Point_3(2, 0,0) ),
                    Segment_3( Point_3(-2,0,0),Point_3(0, 0,0) ),
                    Segment_3( Point_3( 2,0,0),Point_3(0, 0,0) )
  };

  for (int i=0;i<4;++i)
    for (int j=i+1;j<5;++j)
      CGAL_assertion( CGAL::do_intersect(s[i],s[j]) );
  
  CGAL::Object objs[6];
  objs[0] = CGAL::intersection(s[0],s[1]); //Point(0,0,0)
  objs[1] = CGAL::intersection(s[0],s[2]); //Point(0,0,0)
  objs[2] = CGAL::intersection(s[0],s[3]); //Point(0,0,0)
  objs[3] = CGAL::intersection(s[3],s[4]); //Point(0,0,0)
  objs[4] = CGAL::intersection(s[1],s[2]); //s[1]
  objs[5] = CGAL::intersection(s[1],s[3]); //Segment_3( Point(-1,0,0),Point(0,0,0) )

  for (int k=0;k<4;++k){
    const Point_3* p=CGAL::object_cast<Point_3>(&(objs[k]));
    CGAL_assertion(p!=NULL);
    CGAL_assertion(*p==Point_3(0,0,0));
  }
  
  const Segment_3* seg=CGAL::object_cast<Segment_3>(&(objs[4]));
  CGAL_assertion(seg!=NULL);
  CGAL_assertion(*seg==s[1]);

  seg=CGAL::object_cast<Segment_3>(&(objs[5]));
  CGAL_assertion(seg!=NULL);
  CGAL_assertion(*seg==Segment_3( Point_3(0,0,0),Point_3(-1,0,0) ) || *seg==Segment_3( Point_3(-1,0,0),Point_3(0,0,0) ));  

  
  CGAL_assertion( !CGAL::do_intersect(s[0],Segment_3(Point_3(-1,0,0),Point_3(-0.5,0,0))) );
  CGAL_assertion( !CGAL::do_intersect(s[4],Segment_3(Point_3(-1,0,0),Point_3(-0.5,0,0))) );
  
  std::cout << "OK!" << std::endl;
}

int main()
{
  K k;
  _test_intersection_construct(k);
  return 0;
}
