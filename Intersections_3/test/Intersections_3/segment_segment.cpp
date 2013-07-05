#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <cassert>

typedef CGAL::Quotient< CGAL::MP_Float >                    FT;
typedef CGAL::Cartesian<FT>                                 K;

// The construction test has to be working
// The equal test has to be working

#define TEST_CASE_SEGMENT_BASE(i,j,k,l,i1,i2) \
{\
  typename K::Segment_3 s1(pts[i],pts[j]); \
  typename K::Segment_3 s2(pts[k],pts[l]); \
  CGAL::Object obj= CGAL::intersection(s1,s2); \
  typename K::Segment_3 s; \
  typename K::Segment_3 res(pts[i1],pts[i2]);\
  if (!CGAL::assign(s,obj) || ( s!=res && s!=res.opposite() )  ) {\
    std::cerr << "ERROR test-case "<< i << j << k <<l <<std::endl; \
    exit (EXIT_FAILURE);\
  }\
}

#define TEST_CASE_SEGMENT(i,j,k,l,i1,i2) \
  TEST_CASE_SEGMENT_BASE(i,j,k,l,i1,i2) \
  TEST_CASE_SEGMENT_BASE(j,i,k,l,i1,i2) \
  TEST_CASE_SEGMENT_BASE(i,j,l,k,i1,i2) \
  TEST_CASE_SEGMENT_BASE(j,i,l,k,i1,i2) \
  TEST_CASE_SEGMENT_BASE(k,l,i,j,i1,i2) \
  TEST_CASE_SEGMENT_BASE(k,l,j,i,i1,i2) \
  TEST_CASE_SEGMENT_BASE(l,k,i,j,i1,i2) \
  TEST_CASE_SEGMENT_BASE(l,k,j,i,i1,i2) 
  
#define TEST_CASE_POINT_BASE(i,j,k,l,ind) \
{\
  typename K::Segment_3 s1(pts[i],pts[j]); \
  typename K::Segment_3 s2(pts[k],pts[l]); \
  CGAL::Object obj= CGAL::intersection(s1,s2); \
  typename K::Point_3 p; \
  if (!CGAL::assign(p,obj) || p!=pts[ind] ) {\
    std::cerr << "ERROR test-case "<< i << j << k <<l <<std::endl; \
    exit (EXIT_FAILURE);\
  }\
}

#define TEST_CASE_POINT(i,j,k,l,ind) \
  TEST_CASE_POINT_BASE(i,j,k,l,ind) \
  TEST_CASE_POINT_BASE(j,i,k,l,ind) \
  TEST_CASE_POINT_BASE(i,j,l,k,ind) \
  TEST_CASE_POINT_BASE(j,i,l,k,ind) \
  TEST_CASE_POINT_BASE(k,l,i,j,ind) \
  TEST_CASE_POINT_BASE(k,l,j,i,ind) \
  TEST_CASE_POINT_BASE(l,k,i,j,ind) \
  TEST_CASE_POINT_BASE(l,k,j,i,ind)
  

#define TEST_CASE_EMPTY_BASE(i,j,k,l) \
{\
  typename K::Segment_3 s1(pts[i],pts[j]); \
  typename K::Segment_3 s2(pts[k],pts[l]); \
  CGAL::Object obj= CGAL::intersection(s1,s2); \
  if ( !obj.empty() ) {\
    std::cerr << "ERROR test-case "<< i << j << k <<l <<std::endl; \
    exit (EXIT_FAILURE);\
  }\
}

#define TEST_CASE_EMPTY(i,j,k,l) \
  TEST_CASE_EMPTY_BASE(i,j,k,l)\
  TEST_CASE_EMPTY_BASE(i,j,l,k)\
  TEST_CASE_EMPTY_BASE(j,i,k,l)\
  TEST_CASE_EMPTY_BASE(j,i,l,k)\
  TEST_CASE_EMPTY_BASE(k,l,i,j)\
  TEST_CASE_EMPTY_BASE(l,k,i,j)\
  TEST_CASE_EMPTY_BASE(k,l,j,i)\
  TEST_CASE_EMPTY_BASE(l,k,j,i)

template <class K>
void all_cases_collinear(typename K::Point_3 pts[4]){
  std::cout << "4 points cases\n";
  TEST_CASE_EMPTY(0,1,2,3)
  TEST_CASE_SEGMENT(0,2,1,3,1,2)
  TEST_CASE_SEGMENT(0,3,2,1,1,2)
  std::cout << "3 points cases\n";
  TEST_CASE_SEGMENT(0,1,0,2,0,1)
  TEST_CASE_SEGMENT(0,2,1,2,1,2)
  TEST_CASE_POINT(1,2,0,1,1)
  std::cout << "2 points cases\n";
  TEST_CASE_SEGMENT(1,2,1,2,1,2)
}

template <class K>
void _test_intersection_construct(K)
{
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
    for (int j=i+1;j<5;++j){
      assert( CGAL::do_intersect(s[i],s[j]) );
      assert( CGAL::do_intersect(s[i].supporting_line(),s[j]) );
      assert( CGAL::do_intersect(s[i],s[j].supporting_line()) );
    }
  
  CGAL::Object objs[8];
  objs[0] = CGAL::intersection(s[0],s[1]); //Point(0,0,0)
  objs[1] = CGAL::intersection(s[0],s[2]); //Point(0,0,0)
  objs[2] = CGAL::intersection(s[0],s[3]); //Point(0,0,0)
  objs[3] = CGAL::intersection(s[3],s[4]); //Point(0,0,0)
  objs[4] = CGAL::intersection(s[0].supporting_line(),s[1]); //Point(0,0,0)    
  objs[5] = CGAL::intersection(s[1],s[2]); //s[1]
  objs[6] = CGAL::intersection(s[1],s[3]); //Segment_3( Point(-1,0,0),Point(0,0,0) )
  objs[7] = CGAL::intersection(s[2],s[1].supporting_line()); //s[2]    
  for (int k=0;k<5;++k){
    const Point_3* p=CGAL::object_cast<Point_3>(&(objs[k]));
    assert(p!=NULL);
    assert(*p==Point_3(0,0,0));
  }
  
  const Segment_3* seg=CGAL::object_cast<Segment_3>(&(objs[5]));
  assert(seg!=NULL);
  assert(*seg==s[1]);
  seg=CGAL::object_cast<Segment_3>(&(objs[7]));
  assert(*seg==s[2]);  
  
  seg=CGAL::object_cast<Segment_3>(&(objs[6]));
  assert(seg!=NULL);
  assert(*seg==Segment_3( Point_3(0,0,0),Point_3(-1,0,0) ) || *seg==Segment_3( Point_3(-1,0,0),Point_3(0,0,0) ));  

  
  assert( !CGAL::do_intersect(s[0],Segment_3(Point_3(-1,0,0),Point_3(-0.5,0,0))) );
  assert( !CGAL::do_intersect(s[4],Segment_3(Point_3(-1,0,0),Point_3(-0.5,0,0))) );
  
  std::cout << "OK!" << std::endl;
}

int main()
{
  K k;
  _test_intersection_construct(k);
  
   
  K::Point_3 pts[4] = {K::Point_3(0,0,0),K::Point_3(1,0,0),K::Point_3(2,0,0),K::Point_3(3,0,0)};
  all_cases_collinear<K>(pts);
  return 0;
}
