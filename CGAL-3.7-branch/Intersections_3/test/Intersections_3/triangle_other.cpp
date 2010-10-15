#include <cassert>


#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Triangle_2_Triangle_2_do_intersect.h>
#include <CGAL/Triangle_3_Line_3_do_intersect.h>
#include <CGAL/Triangle_3_Plane_3_do_intersect.h>
#include <CGAL/Triangle_3_Segment_3_do_intersect.h>
#include <CGAL/Triangle_3_Ray_3_do_intersect.h>
#include <CGAL/Triangle_3_Triangle_3_do_intersect.h>
#include <CGAL/Triangle_3_Tetrahedron_3_do_intersect.h>

#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

//typedef double Coord_type;
typedef CGAL::Lazy_exact_nt<CGAL::Quotient<CGAL::MP_Float> > Coord_type;
typedef CGAL::Cartesian<Coord_type>  Rep;

typedef Rep::Point_3         Point;
typedef Rep::Segment_3       Segment;
typedef Rep::Triangle_3      Triangle;
typedef Rep::Tetrahedron_3   Tetrahedron;
typedef Rep::Line_3          Line;
typedef Rep::Plane_3         Plane;
typedef Rep::Ray_3           Ray;
typedef Rep::Vector_3        Vector;

typedef Rep::Triangle_2     Triangle_2;
typedef Rep::Point_2        Point_2;


bool tl_correct_answer(int k, int i) {
    switch ( k ) {
    case 0:
      return (( i == 1) || ( i == 8 ) || (i > 10 ));
    case 1:
       return (( i == 2 ) || (i > 10 ));
    case 2:
      return true;
    case 3:
      return (( i != 1) && ( i != 6) && ( i != 7) && ( i != 8));
    case 4:
      return true;
    default:
      return false;
    }
  }

bool ts_correct_answer(int k, int i) {
    switch ( k ) {
    case 0:
      return (( i == 8) || (i > 10 ));
    case 1:
       return (i > 10 );
    case 2:
      return true;
    case 3:
      return (( i == 4) || ( i == 5) || ( i > 8));
    case 4:
      return true;
    default:
      return false;
    }
  }

bool tr_correct_answer(int k, int i) {
    switch ( k ) {
    case 0:
      return ( i == 1);
    case 1:
       return (i == 2);
    case 2:
      return true;
    case 3:
      return (( i == 0) || (i == 2) || ( i == 3) );

    case 4:
      return true;
    default:
      return false;
    }
  }


bool _tr_correct_answer(int k, int i) {
    switch ( k ) {
    case 0:
      return (( i == 1) || ( i == 8) || (i > 10 ));
    case 1:
       return ((i == 2) || (i > 10 ));
    case 2:
      return true;
    case 3:
      return (( i == 0) || (i == 2) || ( i == 3)
	      || ( i == 4) || ( i == 5 ) || ( i > 8));
    case 4:
      return true;
    default:
      return false;
    }
  }


Point_2 project(const Point & p) {
  return Point_2(p.x(),p.y());
}

Triangle_2 project (const Triangle & t) {
  return Triangle_2(project(t.vertex(0)),project(t.vertex(1)),project(t.vertex(2)));
}


int
main()
{

  /////////////////////////////////////////////////////////////////////
  //
  // Triangle Line Predicate
  //
  /////////////////////////////////////////////////////////////////////




  Point a(0,0,2), b(4,4,2), c(4,-4,2);
  Triangle t(a,b,c);

  Point  p(0,0,0);

  Triangle tris[6] = {
    Triangle(a,b,c), Triangle(a,c,b), Triangle(b,a,c),
    Triangle(b,c,a), Triangle(c,a,b), Triangle(c,b,a)
  };

  Point pts[19] = {
    Point(-1,1,1),
    Point(-1,-1,1),
    Point(-1,0,1),
    Point(0,1,1),
    Point(0,-1,1),
    Point(2,-3,1),
    Point(2,3,1),
    Point(4,6,1),
    Point(4,4,1),
    Point(4,-6,1),
    Point(4,-4,1),
    Point(4,0,1),
    Point(0,0,1),
    Point(2,-2,1),
    Point(2,2,1),
    Point(1,0,1),
    Point(1,-1,1),
    Point(2,0,1),
    Point(1,1,1)
  };

   Point copl[5] = {
    Point(-1,-1,2),
    Point(-1,0,2),
    Point(0,0,2),
    Point(0,1,2),
    Point(1,1,2)
  };

   std::cout <<"Testing Triangle-Line Overlap Test ...."<<std::endl ;

   std::cout <<"Three dimensional case " <<std::endl;
     for (int k = -3 ; k < 4 ; k = k + 2 ) {
       for(int i = 0; i<19 ;++i) {
	
	 Line l(p,p+k*(pts[i]-p));
	
	 for ( int j = 0 ; j < 6 ; j++) {
	   assert(CGAL::do_intersect(l,tris[j]) == (i<12 ? false : true));
	   assert(CGAL::do_intersect(l.opposite(),tris[j])
		  == (i<12 ? false : true));
	 }
       }
     }

  // coplanar cases

     std::cout <<"Coplanar case"<<std::endl;

     for (int k = 0 ; k < 5 ; k ++ ) {
    for(int i = 0; i<19 ;++i) {
        Line l(copl[k],  CGAL::ORIGIN + 2*(pts[i]- CGAL::ORIGIN) );

	bool answer = tl_correct_answer(k,i);


      for ( int j = 0 ; j < 6 ; j++) {
	if (! l.is_degenerate() ) {
	  assert(CGAL::do_intersect(l,tris[j]) ==  answer);
	  assert(CGAL::do_intersect(l.opposite(),tris[j]) == answer);
	}
      }
    }
  }


  /////////////////////////////////////////////////////////////////////
  //
  // Triangle Plane Predicate
  //
  /////////////////////////////////////////////////////////////////////



  Plane h(0,1,0,5);
  std::cout <<"Testing Triangle-Plane Overlap Test ...."<<std::endl ;

  std::cout <<"Three dimensional case " <<std::endl;
  for (int j = 0 ; j < 6 ; j++ ) {
    assert(CGAL::do_intersect(h,tris[j]) == false);
    assert(CGAL::do_intersect(h.opposite(),tris[j]) == false);
  }

  h = Plane(0,1,0,4);
  for (int j = 0 ; j < 6 ; j++ ) {
    assert(CGAL::do_intersect(h,tris[j]) == true);
    assert(CGAL::do_intersect(h.opposite(),tris[j]) == true);
  }

  h = Plane(0,1,0,2);
  for (int j = 0 ; j < 6 ; j++ ) {
    assert(CGAL::do_intersect(h,tris[j]) == true);
    assert(CGAL::do_intersect(h.opposite(),tris[j]) == true);
  }

  h = Plane(0,1,0,0);
  for (int j = 0 ; j < 6 ; j++ ) {
    assert(CGAL::do_intersect(h,tris[j]) == true);
    assert(CGAL::do_intersect(h.opposite(),tris[j]) == true);
  }

  h = Plane(1,0,0,-4);
  for (int j = 0 ; j < 6 ; j++ ) {
    assert(CGAL::do_intersect(h,tris[j]) == true);
    assert(CGAL::do_intersect(h.opposite(),tris[j]) == true);
  }


  // coplanar cases
  std::cout <<"Coplanar case"<<std::endl ;

  h = Plane(0,0,1,-2);
  for (int j = 0 ; j < 6 ; j++ ) {
    assert(CGAL::do_intersect(h,tris[j]) == true);
    assert(CGAL::do_intersect(h.opposite(),tris[j]) == true);
  }

   /////////////////////////////////////////////////////////////////////
  //
  // Triangle Segment Predicate
  //
  /////////////////////////////////////////////////////////////////////

  std::cout <<"Testing Triangle-Segment Overlap Test ...."<<std::endl ;

  std::cout <<"Three dimensional case " <<std::endl;

  for (int k = -3 ; k < 4 ; k = k + 2 ) {
    for(int i = 0; i<19 ;++i) {

      Segment s(p,p+k*(pts[i]-p));
      	
      for ( int j = 0 ; j < 6 ; j++) {
	assert(CGAL::do_intersect(s,tris[j]) ==
	       (( i < 12 || k < 2 ) ? false : true));
	assert(CGAL::do_intersect(s.opposite(),tris[j])
	       == ((( i < 12 ) || ( k < 2 ) ) ? false : true));
      }
    }
  }

  // coplanar cases
  std::cout <<"Coplanar case"<<std::endl ;

  for (int k = 0 ; k < 5 ; k ++ ) {
    for(int i = 0; i<19 ;++i) {
        Segment s(copl[k],  CGAL::ORIGIN + 2*(pts[i]- CGAL::ORIGIN) );

	bool answer = ts_correct_answer(k,i);


      for ( int j = 0 ; j < 6 ; j++) {
	if (! s.is_degenerate() ) {
	  assert(CGAL::do_intersect(s,tris[j]) ==  answer);
	  assert(CGAL::do_intersect(s.opposite(),tris[j]) == answer);
	}
      }
    }
  }

  /////////////////////////////////////////////////////////////////////
  //
  // Triangle Ray Predicate
  //
  /////////////////////////////////////////////////////////////////////

  std::cout <<"Testing Triangle-Ray Overlap Test ...."<<std::endl ;

  std::cout <<"Three dimensional case " <<std::endl;

  for (int k = -3 ; k < 4 ; k = k + 2 ) {
    for(int i = 0; i<19 ;++i) {
      Ray r(p,p+k*(pts[i]-p));
      	
      for ( int j = 0 ; j < 6 ; j++) {
	assert(CGAL::do_intersect(r,tris[j]) ==
	       (( i < 12 || k < 0 ) ? false : true));
	assert(CGAL::do_intersect(r.opposite(),tris[j])
	       == ((( i < 12 ) || ( k > 0)) ? false : true));
      }
    }
  }

  // coplanar cases
  std::cout <<"Coplanar case"<<std::endl ;

  for (int k = 0 ; k < 5 ; k ++ ) {
    for(int i = 0; i<19 ;++i) {
        Ray r(copl[k],  CGAL::ORIGIN + 2*(pts[i]- CGAL::ORIGIN) );

	bool answer = ts_correct_answer(k,i);
	bool answer_opposite = tr_correct_answer(k,i);

      for ( int j = 0 ; j < 6 ; j++) {
	if (! r.is_degenerate() ) {
	  assert(CGAL::do_intersect(r,tris[j]) ==  answer);
	  assert(CGAL::do_intersect(r.opposite(),tris[j])
		 == answer_opposite);
	}
      }
    }
  }

  // parallel case
  Ray r(Point(-1,0,1),Point(4,0,1));
  assert(CGAL::do_intersect(r,tris[0]) ==  false);

  /////////////////////////////////////////////////////////////////////
  //
  // Triangle Tetrahedron Predicate
  //
  /////////////////////////////////////////////////////////////////////

  std::cout <<"Testing Triangle-Tetrahedron Overlap Test ...."<<std::endl ;

  {
    Tetrahedron tet(Point(1.0962000000000000632383034826489165425300598144531,
			  0.97450000000000003286260152890463359653949737548828,
			  -1.6816999999999999726441046732361428439617156982422),
	 	    Point(0.70588000000000006295408638834487646818161010742188,
			  0.83156000000000007688072400924284011125564575195312,
			  -1.3078650000000000552802248421357944607734680175781),
		    Point(0.33200527798677881285982493864139541983604431152344,
			  1.6713309295164104906206148370984010398387908935547,
			  -1.3771306515070147469259609351865947246551513671875),
		    Point(0.29692468967934992907231617209617979824542999267578,
			  1.654795228488893599205766804516315460205078125,
			  -1.4200807988528463265964774109306745231151580810547));
    Triangle tr(Point(0.4993996236189303661312521853687940165400505065918,
		      1.2953416022081787328801283365464769303798675537109,
		      -1.3461184494882023621187272510724142193794250488281),
		Point(0.48104952476832690821950677673157770186662673950195,
		      1.2841482485875792551865970381186343729496002197266,
		      -1.3695576366961845771186290221521630883216857910156),
		Point(0.53749717115046513615794765428290702402591705322266,
		      1.4500343543104132759680169328930787742137908935547,
		      -1.4988251058616692823477478668792173266410827636719));
    assert(CGAL::do_intersect(tet,tr) == true);

    tet = Tetrahedron(Point(0,0,0), Point(1,0,0), Point(0,1,0), Point(0,0,1));

    assert(CGAL::do_intersect(tet,
		 Triangle(Point(-1,-1,.5), Point(-1,2,.5), Point(2,-1,.5))) == true);
    assert(CGAL::do_intersect(tet,
		 Triangle(Point(0,0,.5), Point(.5,0,.5), Point(0,.5,.5))) == true);
    assert(CGAL::do_intersect(tet,
		 Triangle(Point(1,1,.5), Point(.5,0,.5), Point(0,.5,.5))) == true);
    assert(CGAL::do_intersect(tet,
		 Triangle(Point(.2,.3,.5), Point(1,1,1), Point(2,1,1))) == true);
    assert(CGAL::do_intersect(tet,
		 Triangle(Point(-1,0,1), Point(-1,1,0), Point(1,0,-1))) == true);
  }

  /////////////////////////////////////////////////////////////////////
  //
  // Triangle Triangle Predicate
  //
  /////////////////////////////////////////////////////////////////////


  std::cout <<"Testing three-dimensional Triangle-Triangle Overlap Test ...."<<std::endl ;

  // all tris1 lie on the plane z = 0;

  Triangle tris1[7] = {
    // included in the halfspace y < 0  does not intersect the x axis
    Triangle(Point(0,-3,0),Point(4,-5,0),Point(4,-1,0)),
    // one vertex  y = 0, others y < 0 intersects the x axis in 0
    Triangle(Point(0,0,0),Point(0,-4,0),Point(4,-2,0)),
    // one vertex y > 0 others y < 0 intersect the x axis in [ 0 1 ]
    Triangle(Point(0,1,0),Point(0,-3,0),Point(2,-1,0)),
    Triangle(Point(-1,-1,0),Point(1,-3,0),Point(1,1,0)),
    // one vertex y < 0 one y = 0 one y > 0 intersect the x axis in [ 0 2 ]
    Triangle(Point(0,2,0),Point(0,-2,0),Point(2,0,0)),
    Triangle(Point(0,0,0),Point(2,2,0),Point(2,-2,0)),
    // two vertices y = 0, one y < 0 intersects the x axis in [ 0 3 ]
    Triangle(Point(0,0,0),Point(3,0,0),Point(1.5,-1.5,0))

  };

  // all tris2 lie on the plane y = 0;
  Triangle tris2[7] = {
     // included in the halfspace z > 0  does not intersect the x axis
    Triangle(Point(0,0,3),Point(4,0,5),Point(4,0,1)),
    // one vertex  z = 0, others z > 0 intersects the x axis in 0
    Triangle(Point(0,0,0),Point(0,0,4),Point(4,0,2)),
    // one vertex z < 0 others z > 0 intersect the x axis in [ 0 1 ]
    Triangle(Point(0,0,-1),Point(0,0,3),Point(2,0,1)),
    Triangle(Point(-1,0,1),Point(1,0,3),Point(1,0,-1)),
    // one vertex z < 0 one z = 0 one z > 0 intersect the x axis in [ 0 2 ]
    Triangle(Point(0,0,-2),Point(0,0,2),Point(2,0,0)),
    Triangle(Point(0,0,0),Point(2,0,-2),Point(2,0,2)),
    // two vertices z = 0, one z > 0 intersects the x axis in [ 0 3 ]
    Triangle(Point(0,0,0),Point(3,0,0),Point(1.5,0,1.5))

  };

  Coord_type x_translation[9] = { -4, -3, -2 , -1, 0 , 1, 2, 3, 4 };

  // answer = i != 0 and j != 0
  // Translation of input1 by vector (-4,0,0)  no intersection
  // Translation  of input1 by vector (-3,0,0), answer && i > 6
  // Translation  of input1 by vector (-2,0,0), answer && i > 3
  // Translation  of input1 by vector (-1,0,0), answer && i > 1
  // No translation , answer
  // Translation  of input1 by vector (1,0,0), answer && j > 1
  // Translation  of input1 by vector (2,0,0), answer && j > 3
  // Translation  of input1 by vector (3,0,0), answer && j > 5
  // Translation  of input1 by vector (4,0,0), no intersection

  int i_condition[9] = { 8 , 5 , 3 , 1 , -1 , -1 , -1, -1, -1 };
  int j_condition[9] = { -1 , -1 , -1 , -1 , -1 , 1 , 3, 5, 8 };


  for ( int tr = 0 ; tr < 9 ; tr++ ) {
    for (int i = 0 ; i < 7 ; i++ ) {
      for (int j = 0 ; j < 7 ; j++) {

	Triangle input1 = tris1[i];
	Triangle input2 = tris2[j];
	
	Point v1a = input1.vertex(0) + Vector(x_translation[tr],0,0);
	Point v1b = input1.vertex(1)+  Vector(x_translation[tr],0,0);
	Point v1c = input1.vertex(2)+  Vector(x_translation[tr],0,0);
	Point v2a = input2.vertex(0);
	Point v2b = input2.vertex(1);
	Point v2c = input2.vertex(2);
	
	bool answer = ( ( i != 0 ) && ( j != 0 ) )
	  && ( ( i > i_condition[tr]) && ( j > j_condition[tr] ) );
	
	
	input2 = Triangle(v2a,v2b,v2c);

	assert(CGAL::do_intersect(Triangle(v1a,v1b,v1c),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1a,v1b,v1c))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1a,v1c,v1b),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1a,v1c,v1b))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1b,v1a,v1c),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1b,v1a,v1c))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1b,v1c,v1a),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1b,v1c,v1a))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1c,v1a,v1b),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1c,v1a,v1b))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1c,v1b,v1a),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1c,v1b,v1a))
	       == answer);
	
	input2 = Triangle(v2a,v2c,v2b);

	assert(CGAL::do_intersect(Triangle(v1a,v1b,v1c),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1a,v1b,v1c))
	       == answer);		
	assert(CGAL::do_intersect(Triangle(v1a,v1c,v1b),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1a,v1c,v1b))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1b,v1a,v1c),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1b,v1a,v1c))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1b,v1c,v1a),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1b,v1c,v1a))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1c,v1a,v1b),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1c,v1a,v1b))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1c,v1b,v1a),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1c,v1b,v1a))
	       == answer);
	
	input2 = Triangle(v2b,v2a,v2c);
	
	assert(CGAL::do_intersect(Triangle(v1a,v1b,v1c),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1a,v1b,v1c))
	       == answer);		
	assert(CGAL::do_intersect(Triangle(v1a,v1c,v1b),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1a,v1c,v1b))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1b,v1a,v1c),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1b,v1a,v1c))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1b,v1c,v1a),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1b,v1c,v1a))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1c,v1a,v1b),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1c,v1a,v1b))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1c,v1b,v1a),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1c,v1b,v1a))
	       == answer);
	
	input2 = Triangle(v2b,v2c,v2a);
	
	assert(CGAL::do_intersect(Triangle(v1a,v1b,v1c),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1a,v1b,v1c))
	       == answer);	
	assert(CGAL::do_intersect(Triangle(v1a,v1c,v1b),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1a,v1c,v1b))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1b,v1a,v1c),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1b,v1a,v1c))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1b,v1c,v1a),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1b,v1c,v1a))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1c,v1a,v1b),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1c,v1a,v1b))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1c,v1b,v1a),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1c,v1b,v1a))
	       == answer);	
	
	input2 = Triangle(v2c,v2a,v2b);
	
	assert(CGAL::do_intersect(Triangle(v1a,v1b,v1c),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1a,v1b,v1c))
	       == answer);		
	assert(CGAL::do_intersect(Triangle(v1a,v1c,v1b),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1a,v1c,v1b))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1b,v1a,v1c),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1b,v1a,v1c))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1b,v1c,v1a),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1b,v1c,v1a))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1c,v1a,v1b),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1c,v1a,v1b))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1c,v1b,v1a),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1c,v1b,v1a))
	       == answer);
	
	input2 = Triangle(v2c,v2b,v2a);
	
	assert(CGAL::do_intersect(Triangle(v1a,v1b,v1c),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1a,v1b,v1c))
	       == answer);		
	assert(CGAL::do_intersect(Triangle(v1a,v1c,v1b),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1a,v1c,v1b))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1b,v1a,v1c),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1b,v1a,v1c))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1b,v1c,v1a),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1b,v1c,v1a))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1c,v1a,v1b),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1c,v1a,v1b))
	       == answer);
	assert(CGAL::do_intersect(Triangle(v1c,v1b,v1a),input2)
	       == answer);
	assert(CGAL::do_intersect(input2,Triangle(v1c,v1b,v1a))
	       == answer);		
      }
    }
  }
  std::cout <<"Testing two-dimensional Triangle-Triangle Overlap Test"<<std::endl;
  std::cout <<"and three-dimensional coplanar pairs ..."<<std::endl ;

   Triangle tris1tt[6] = {
     Triangle(Point(0,0,1), Point(4,-4,1),Point(4,4,1)),
     Triangle(Point(0,0,1), Point(4,4,1),Point(4,-4,1)),
     Triangle(Point(4,-4,1), Point(0,0,1),Point(4,4,1)),
     Triangle(Point(4,-4,1), Point(4,4,1),Point(0,0,1)),
     Triangle(Point(4,4,1), Point(4,-4,1),Point(0,0,1)),
     Triangle(Point(4,4,1), Point(0,0,1),Point(4,-4,1))
   };



  Point ptstt2[50] = {
    Point(7,5,1), Point(6,5,1), Point(6,6,1), Point(6,7,1),
    Point(5,5,1), Point(5,7,1),
    Point(4,6,1), Point(4,7,1),
    Point(2,5,1), Point(2,6,1),Point(2,7,1),
    Point(1,4,1), Point(1,5,1),
    Point(0,3,1), Point(0,4,1), Point(0,4.5,1), Point(0,6,1),Point(0,7,1),
    Point(-1,3,1),
    Point(-1,5,1),Point(-1,7,1),Point(-1,8,1),
    Point(-2,7,1),
    Point(-5,5,1),Point(-5,4,1),


    Point(0,-3,1),Point(0,-5,1),
    Point(2,-5,1), Point(4,-7,1),Point(4,-5,1),Point(5,-7,1),

    Point(6,-6,1),  Point(5,-5,1),

    Point(0,0,1),
    Point(1,1,1),Point(1,-1,1),
    Point(2,2,1),Point(2,1,1),Point(2,0,1),
    Point(3,3,1),Point(3,0,1),Point(3,-3,1),
    Point(4,4,1), Point(4,1,1),
    Point(4,-1,1),Point(4,-4,1),
    Point(6,3,1),Point(6,1,1),Point(6,-1,1),
    Point(8,2,1)
  };

  Point p_overlap[3] = {
    Point(0,0,1), Point(2,1,1), Point(2,2,1)
  };


  for (int i = 0 ; i < 3; i++) {
    for(int j = 0 ; j < 50 ; j++) {
      for(int k = j + 1 ; k < 50 ; k++) {
	Triangle t1(p_overlap[i], ptstt2[j],ptstt2[k]);
	Triangle t2(p_overlap[i], ptstt2[k],ptstt2[j]);
	Triangle t3(ptstt2[k],p_overlap[i],ptstt2[j]);
	Triangle t4(ptstt2[k],ptstt2[j],p_overlap[i]);
	Triangle t5(ptstt2[j],p_overlap[i],ptstt2[k]);
	Triangle t6(ptstt2[j],ptstt2[k],p_overlap[i]);


	
	if (!t1.is_degenerate()) {
	  for (int l = 0 ; l < 6 ; l++) {
	    assert(CGAL::do_intersect(tris1tt[l],t1) == true);
	    assert(CGAL::do_intersect(project(tris1tt[1]),project(t1)) == true);
	    assert(CGAL::do_intersect(t1,tris1tt[l]) == true);
	    assert(CGAL::do_intersect(project(t1),project(tris1tt[l])) == true);
	    assert(CGAL::do_intersect(tris1tt[l],t2) == true);
	    assert(CGAL::do_intersect(project(tris1tt[l]),project(t2)) == true);
	    assert(CGAL::do_intersect(t2,tris1tt[l]) == true);
	    assert(CGAL::do_intersect(project(t2),project(tris1tt[l])) == true);
	    assert(CGAL::do_intersect(tris1tt[l],t3) == true);
	    assert(CGAL::do_intersect(project(tris1tt[l]),project(t3)) == true);
	    assert(CGAL::do_intersect(t3,tris1tt[l]) == true);
	    assert(CGAL::do_intersect(project(t3),project(tris1tt[l])) == true);
	    assert(CGAL::do_intersect(tris1tt[l],t4) == true);
	    assert(CGAL::do_intersect(project(tris1tt[l]),project(t4)) == true);
	    assert(CGAL::do_intersect(t4,tris1tt[l]) == true);
	    assert(CGAL::do_intersect(project(t4),project(tris1tt[l])) == true);
	    assert(CGAL::do_intersect(tris1tt[l],t5) == true);
	    assert(CGAL::do_intersect(project(tris1tt[l]),project(t5)) == true);
	    assert(CGAL::do_intersect(t5,tris1tt[l]) == true);
	    assert(CGAL::do_intersect(project(t5),project(tris1tt[l])) == true);
	    assert(CGAL::do_intersect(tris1tt[l],t6) == true);
	    assert(CGAL::do_intersect(project(tris1tt[l]),project(t6)) == true);
	    assert(CGAL::do_intersect(t6,tris1tt[l]) == true);
	    assert(CGAL::do_intersect(project(t6),project(tris1tt[l])) == true);
	  }
	}
      }
    }
  }



  // Tests with point p = Point(0,6,1)

  p = Point(0,6,1);

  for(int j = 0 ; j < 50 ; j++) {
    for(int k = j + 1 ; k < 50 ; k++) {
      Triangle t1(p, ptstt2[j],ptstt2[k]);
      Triangle t2(p, ptstt2[k],ptstt2[j]);
      Triangle t3(ptstt2[k],p,ptstt2[j]);
      Triangle t4(ptstt2[k],ptstt2[j],p);
      Triangle t5(ptstt2[j],p,ptstt2[k]);
      Triangle t6(ptstt2[j],ptstt2[k],p);

      if (!t1.is_degenerate()) {
	for (int l = 0 ; l < 6 ; l++) {

	  assert(CGAL::do_intersect(tris1tt[l],t1) ==  (k >= 25));
	  assert(CGAL::do_intersect(project(tris1tt[1]),project(t1)) ==  (k >= 25));
	  assert(CGAL::do_intersect(t1,tris1tt[l]) ==  (k >= 25));
	  assert(CGAL::do_intersect(project(t1),project(tris1tt[l])) ==  (k >= 25));
	  assert(CGAL::do_intersect(tris1tt[l],t2) ==  (k >= 25));
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t2)) ==  (k >= 25));
	  assert(CGAL::do_intersect(t2,tris1tt[l]) ==  (k >= 25));
	  assert(CGAL::do_intersect(project(t2),project(tris1tt[l])) ==  (k >= 25));
	  assert(CGAL::do_intersect(tris1tt[l],t3) ==  (k >= 25));
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t3)) ==  (k >= 25));
	  assert(CGAL::do_intersect(t3,tris1tt[l]) ==  (k >= 25));
	  assert(CGAL::do_intersect(project(t3),project(tris1tt[l])) ==  (k >= 25));
	  assert(CGAL::do_intersect(tris1tt[l],t4) ==  (k >= 25));
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t4)) ==  (k >= 25));
	  assert(CGAL::do_intersect(t4,tris1tt[l]) ==  (k >= 25));
	  assert(CGAL::do_intersect(project(t4),project(tris1tt[l])) ==  (k >= 25));
	  assert(CGAL::do_intersect(tris1tt[l],t5) ==  (k >= 25));
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t5)) ==  (k >= 25));
	  assert(CGAL::do_intersect(t5,tris1tt[l]) ==  (k >= 25));
	  assert(CGAL::do_intersect(project(t5),project(tris1tt[l])) ==  (k >= 25));
	  assert(CGAL::do_intersect(tris1tt[l],t6) ==  (k >= 25));
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t6)) ==  (k >= 25));
	  assert(CGAL::do_intersect(t6,tris1tt[l]) ==  (k >= 25));
	  assert(CGAL::do_intersect(project(t6),project(tris1tt[l])) ==  (k >= 25));
	
	}
      }
    }
  }

  // Tests with point p = Point(-5,5,1)

  p =  Point(-5.5,5.5,1);

  for(int j = 0 ; j < 50 ; j++) {
    for(int k = j + 1 ; k < 50 ; k++) {
      Triangle t1(p, ptstt2[j],ptstt2[k]);
      Triangle t2(p, ptstt2[k],ptstt2[j]);
      Triangle t3(ptstt2[k],p,ptstt2[j]);
      Triangle t4(ptstt2[k],ptstt2[j],p);
      Triangle t5(ptstt2[j],p,ptstt2[k]);
      Triangle t6(ptstt2[j],ptstt2[k],p);


      bool answer = !(( k < 25) ||
		     (( k <= 26) && (j >=18) && (j <= 22)) ||
		     (( k <= 30) && (j >=23)));

      if (!t1.is_degenerate()) {
	for (int l = 0 ; l < 6 ; l++) {
	  assert(CGAL::do_intersect(tris1tt[l],t1) ==  answer);
	  assert(CGAL::do_intersect(project(tris1tt[1]),project(t1))== answer);
	  assert(CGAL::do_intersect(t1,tris1tt[l]) ==  answer);
	  assert(CGAL::do_intersect(project(t1),project(tris1tt[l])) ==  answer);
	  assert(CGAL::do_intersect(tris1tt[l],t2) ==  answer);
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t2)) ==  answer);
	  assert(CGAL::do_intersect(t2,tris1tt[l]) ==  answer);
	  assert(CGAL::do_intersect(project(t2),project(tris1tt[l])) ==  answer);
	  assert(CGAL::do_intersect(tris1tt[l],t3) ==  answer);
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t3)) ==  answer);
	  assert(CGAL::do_intersect(t3,tris1tt[l]) ==  answer);
	  assert(CGAL::do_intersect(project(t3),project(tris1tt[l])) ==  answer);
	  assert(CGAL::do_intersect(tris1tt[l],t4) ==  answer);
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t4)) ==  answer);
	  assert(CGAL::do_intersect(t4,tris1tt[l]) ==  answer);
	  assert(CGAL::do_intersect(project(t4),project(tris1tt[l])) ==  answer);
	  assert(CGAL::do_intersect(tris1tt[l],t5) ==  answer);
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t5)) ==  answer);
	  assert(CGAL::do_intersect(t5,tris1tt[l]) ==  answer);
	  assert(CGAL::do_intersect(project(t5),project(tris1tt[l])) ==  answer);
	  assert(CGAL::do_intersect(tris1tt[l],t6) ==  answer);
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t6)) ==  answer);
	  assert(CGAL::do_intersect(t6,tris1tt[l]) ==  answer);
	  assert(CGAL::do_intersect(project(t6),project(tris1tt[l])) ==  answer);

	}
      }
    }
  }

  // Tests with point p = Point(-5,4,1)

  p =  Point(-5,4,1);

  for(int j = 0 ; j < 50 ; j++) {
    for(int k = j + 1 ; k < 50 ; k++) {
      Triangle t1(p, ptstt2[j],ptstt2[k]);
      Triangle t2(p, ptstt2[k],ptstt2[j]);
      Triangle t3(ptstt2[k],p,ptstt2[j]);
      Triangle t4(ptstt2[k],ptstt2[j],p);
      Triangle t5(ptstt2[j],p,ptstt2[k]);
      Triangle t6(ptstt2[j],ptstt2[k],p);


      bool answer = !(( k < 25) ||
		      (( k <= 26) && (j >=18) && (j <= 22)) ||
		      (( k <= 30) && (j >=23)) ||
		      (( k <= 32) && (j >=24)));

       if (!t1.is_degenerate()) {
	for (int l = 0 ; l < 6 ; l++) {
	  assert(CGAL::do_intersect(tris1tt[l],t1) ==  answer);
	  assert(CGAL::do_intersect(project(tris1tt[1]),project(t1))== answer);
	  assert(CGAL::do_intersect(t1,tris1tt[l]) ==  answer);
	  assert(CGAL::do_intersect(project(t1),project(tris1tt[l])) ==  answer);
	  assert(CGAL::do_intersect(tris1tt[l],t2) ==  answer);
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t2)) ==  answer);
	  assert(CGAL::do_intersect(t2,tris1tt[l]) ==  answer);
	  assert(CGAL::do_intersect(project(t2),project(tris1tt[l])) ==  answer);
	  assert(CGAL::do_intersect(tris1tt[l],t3) ==  answer);
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t3)) ==  answer);
	  assert(CGAL::do_intersect(t3,tris1tt[l]) ==  answer);
	  assert(CGAL::do_intersect(project(t3),project(tris1tt[l])) ==  answer);
	  assert(CGAL::do_intersect(tris1tt[l],t4) ==  answer);
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t4)) ==  answer);
	  assert(CGAL::do_intersect(t4,tris1tt[l]) ==  answer);
	  assert(CGAL::do_intersect(project(t4),project(tris1tt[l])) ==  answer);
	  assert(CGAL::do_intersect(tris1tt[l],t5) ==  answer);
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t5)) ==  answer);
	  assert(CGAL::do_intersect(t5,tris1tt[l]) ==  answer);
	  assert(CGAL::do_intersect(project(t5),project(tris1tt[l])) ==  answer);
	  assert(CGAL::do_intersect(tris1tt[l],t6) ==  answer);
	  assert(CGAL::do_intersect(project(tris1tt[l]),project(t6)) ==  answer);
	  assert(CGAL::do_intersect(t6,tris1tt[l]) ==  answer);
	  assert(CGAL::do_intersect(project(t6),project(tris1tt[l])) ==  answer);

	}
      }
    }
  }




  return 0;
}


