#include <predicates.h>
int main(int, char *[]) {

  initialize_rational_points();
 
  





  {
    Summary sum;
    std::cout << "\n\nTesting center_bounded_side_of_sphere_c" << std::endl;
    std::vector<S> ss;
    ss.push_back(S(P(0,0,0), 1));
    std::cout << ss[0] << std::endl;
    for (int i=0; i< 50; ++i) {
      ss.push_back(random_sphere(ss.back().bbox()));
      V2 v= random_vector_on_unit_circle();
      if (i%2==0) {
	ss.push_back(S(P(v[0], 10, v[1]), random_coordinate()));
      } else {
	ss.push_back(S(P(v[0], v[1], 10), random_coordinate()));
      }
    }
    Traits tr(ss.begin(), ss.end());
    CGAL::Bbox_3 bb= ss.back().bbox();
    for (unsigned int i=1; i< ss.size(); ++i) {
      {
	FT x= random_coordinate(bb.xmin(), bb.xmax());
	C2 c2;
	if (slice(ss[0], x, c2)){
	  P pt(x, ss[i].center().y(), ss[i].center().z());
	  SP sp= random_sp(x);
	  //SP spx= random_sp(ss[i].center().x());
	  
	  CGAL::Bounded_side mbs0= tr.center_bounded_side_of_circle_c(K(i),sp, K(0), plane_coordinate(0));
	  CGAL::Bounded_side mbs1= tr.center_bounded_side_of_circle_c(K(i),sp, K(0), plane_coordinate(1));
	  CGAL::Bounded_side cbs0= check_bounded_side(c2, pt, plane_coordinate(0));
	  CGAL::Bounded_side cbs1= check_bounded_side(c2, pt, plane_coordinate(1));
	  sum(mbs0);
	  sum(mbs1);
	  
	  CGAL_assertion(mbs0==cbs0);
	  CGAL_assertion(mbs1==cbs1);
	}
      }
      if (i%2 == 0) {
	C2 c2;
	if (slice(ss[0], ss[i].center().x(), c2)){
	  P pt= ss[i].center();
	  SP sp= random_sp(ss[i].center().x());
	  //SP spx= random_sp(ss[i].center().x());
	  
	  CGAL::Bounded_side mbs0= tr.center_bounded_side_of_circle_c(K(i),sp, K(0), 
								      plane_coordinate(0));
	  CGAL::Bounded_side mbs1= tr.center_bounded_side_of_circle_c(K(i),sp, K(0), 
								      plane_coordinate(1));
	  CGAL::Bounded_side cbs0= check_bounded_side(c2, pt, plane_coordinate(0));
	  CGAL::Bounded_side cbs1= check_bounded_side(c2, pt, plane_coordinate(1));
	  sum(mbs0);
	  sum(mbs1);
	  
	  CGAL_assertion(mbs0==cbs0);
	  CGAL_assertion(mbs1==cbs1);
	  CGAL_assertion(mbs0== CGAL::ON_BOUNDARY || mbs1 == CGAL::ON_BOUNDARY);
	  //CGAL_assertion(mbs0== CGAL::ON_BOUNDARY);
	  //CGAL_assertion(mbs1== CGAL::ON_BOUNDARY);
	} else {
	  CGAL_assertion(0);
	}
	
      }
    }

  }

  {
    Summary sum;
    std::cout << "\n\nTesting point_bounded_side_of_sphere_c" << std::endl;
    std::vector<S> ss;
    ss.push_back(S(P(0,0,0), 1));
    std::cout << ss[0] << std::endl;
    Traits tr(ss.begin(), ss.end());
    CGAL::Bbox_3 bb= ss.back().bbox();
    for (int i=0; i< 100; ++i) {
      {
	FT x= random_coordinate(bb.xmin(), bb.xmax());
	C2 c2;
	if (slice(ss[0], x, c2)){
	  P pt(x, random_coordinate(bb.ymin(), bb.ymax()),
	       random_coordinate(bb.zmin(), bb.zmax()));
	  SP sp=random_sp(pt);
	  CGAL::Bounded_side mbs0= tr.point_bounded_side_of_sphere_c(sp, K(0), plane_coordinate(0));
	  CGAL::Bounded_side mbs1= tr.point_bounded_side_of_sphere_c(sp, K(0), plane_coordinate(1));
	  CGAL::Bounded_side cbs0= check_bounded_side(c2, pt, plane_coordinate(0));
	  CGAL::Bounded_side cbs1= check_bounded_side(c2, pt, plane_coordinate(1));
	  sum(mbs0); sum(mbs1);
	  CGAL_assertion(mbs0==cbs0);
	  CGAL_assertion(mbs1==cbs1);
	}
      }
      {
	V2 v= random_vector_on_unit_circle();
	P pt;
	if (i%2==0) {
	  pt=P(v[0], 10, v[1]);
	} else {
	  pt=P(v[0], v[1], 10);
	}
	
	C2 c2;
	if (slice(ss[0], pt.x(), c2)){
	  SP sp=random_sp(pt);
	  CGAL::Bounded_side mbs0= tr.point_bounded_side_of_sphere_c(sp, K(0), plane_coordinate(0));
	  CGAL::Bounded_side mbs1= tr.point_bounded_side_of_sphere_c(sp, K(0), plane_coordinate(1));
	  CGAL::Bounded_side cbs0= check_bounded_side(c2, pt, plane_coordinate(0));
	  CGAL::Bounded_side cbs1= check_bounded_side(c2, pt, plane_coordinate(1));
	  sum(mbs0); sum(mbs1);
	  CGAL_assertion(mbs0==cbs0);
	  CGAL_assertion(mbs1==cbs1);
	}

      }
    }
  }

  {
    Summary sum;
    std::cout << "\n\nTesting point_bounded_side_of_sphere" << std::endl;
 
    std::vector<S> ss;
    ss.push_back(S(P(0,2,-2), 1));
    std::cout << ss[0] << std::endl;
    Traits tr(ss.begin(), ss.end());
    CGAL::Bbox_3 bb= ss.back().bbox();
    for (unsigned int i=0; i< 100; ++i) {
      P p= random_point(bb);
      SP sp=random_sp(p);
      std::cout << "SP is " << sp << std::endl;
      CGAL::Bounded_side cbs= ss.back().bounded_side(p);
      CGAL::Bounded_side mbs= tr.point_bounded_side_of_sphere(sp, K(0));
      sum(mbs);
      CGAL_assertion(cbs== mbs);
    }
    for (unsigned int i=0; i< 20; ++i) {
      SP sp= random_sp(ss.back());
      CGAL_assertion(sp.sphere() == ss.front());
      std::cout << "SP is " << sp << std::endl;
      CGAL::Bounded_side mbs= tr.point_bounded_side_of_sphere(sp, K(0));
      sum(mbs);
      CGAL_assertion(mbs== CGAL::ON_BOUNDARY);
    }
  }

  {
    Summary sum, sumr;
    std::cout << "\n\nTesting rules/center_bounded_side_of_sphere" << std::endl;
    {
      
      std::vector<S> ss;
      ss.push_back(S(P(1,3,1), 1));
      for (int i=0; i< 100; ++i) {
	ss.push_back(random_sphere(ss.back().bbox()));
      }
      Traits tr(ss.begin(), ss.end());
      
      for (int i=0; i< 100; ++i) {
	FT x= random_coordinate(ss.front().bbox().xmin(), ss.front().bbox().xmax());
	P cpt(x, ss[i+1].center().y(),ss[i+1].center().z());
	SP sp = random_sp(x);
	CGAL::Bounded_side cbs= ss.front().bounded_side(cpt);
	CGAL::Bounded_side mbs= tr.center_bounded_side_of_sphere(sp, K(i+1), K(0));
	CGAL::Bounded_side rbs= tr.rules_bounded_side_of_sphere(sp, K(i+1), K(i+1), K(0));
	sum(mbs);
	sumr(rbs);
	CGAL_assertion(cbs==mbs);
	CGAL_assertion(rbs==mbs);
      }
    }
    {
      std::vector<S> ss;
      ss.push_back(S(P(-2,1.1,3), 1));
      for (int i=0; i< 20; ++i) {
	ss.push_back(S(random_point_on_sphere(ss.front()), random_radius()));
      }
      Traits tr(ss.begin(), ss.end());
      for (int i=0; i< 20; ++i ) {
	SP sp= random_sp(ss[i+1].center().x());
	CGAL::Bounded_side cbs= ss.front().bounded_side(ss[i+1].center());
	CGAL::Bounded_side mbs= tr.center_bounded_side_of_sphere(sp, K(i+1), K(0));
	CGAL::Bounded_side rbs= tr.rules_bounded_side_of_sphere(sp, K(i+1), K(i+1), K(0));
	sum(mbs);
	sumr(rbs);
	CGAL_assertion(cbs== CGAL::ON_BOUNDARY);
	CGAL_assertion(mbs== CGAL::ON_BOUNDARY);
	CGAL_assertion(rbs==mbs);
      }
    }
  }
  {
    Summary sum;
    std::cout << "\n\nTesting compare_point_to_equipower_line_c" << std::endl;
    for (int j=0; j< 20; ++j) {
      Traits tr;
      std::vector<S> ss;
      do {
	ss.clear();
	ss.push_back(random_sphere(0));
	ss.push_back(random_sphere(0));
	tr= Traits(ss.begin(), ss.end());
      } while (!tr.intersects(K(0), K(1)));
      CGAL::Bbox_3 bb=ss[0].bbox()+ss[1].bbox();
      
      std::cout << "Spheres are " << ss[0] << ": " << ss[1] << std::endl;

      tr= Traits(ss.begin(), ss.end());
      for (unsigned int i=0; i< 20; ++i) {
	FT x;
	C2 ca, cb;

	do {
	  x= random_coordinate(bb.xmin(), bb.xmax());
	} while (!slice(ss[0], x, ca) || !slice(ss[1], x, cb));

	P2 eqp= equipower_point(ca, cb);
	{
	  P cp(x, random_coordinate(bb.ymin(), bb.ymax()), 
	       random_coordinate(bb.zmin(), bb.zmax()));
	  SP pt=random_sp(cp);
	  CGAL::Comparison_result c0= tr.compare_point_to_equipower_line_c(pt,
									   K(0),
									   K(1),
									   plane_x);
	  CGAL::Comparison_result c1= tr.compare_point_to_equipower_line_c(pt,
									   K(0),
									 K(1),
									   plane_y);
	  sum(c0); sum(c1);
	  CGAL::Comparison_result cc0= CGAL::compare(cp[1], eqp.x());
	  CGAL::Comparison_result cc1= CGAL::compare(cp[2], eqp.y());
	  if (cc0 != c0 || cc1 != c1) {
	    std::cout << "Error : " << x << std::endl;
	    std::cout << c0 << " " << c1 << " " << cc0 << " " << cc1 << std::endl;
	    std::cout << ss[0] << std::endl; //0
	    std::cout << ss[1] << std::endl; //1
	    std::cout << x << " " << eqp << " 1/10000 1" << std::endl; //2
	    std::cout << cp << " 1/10000 1" << std::endl;
	    CGAL_assertion(0);
	  }
	}
	{
	  P cp0(x, eqp.x(), 
	       random_coordinate(bb.zmin(), bb.zmax()));
	  SP pt0=random_sp(cp0);
	  CGAL::Comparison_result c0= tr.compare_point_to_equipower_line_c(pt0,
									   K(0),
									   K(1),
									   plane_x);
	  P cp1(x, random_coordinate(bb.ymin(), bb.ymax()), eqp.y());
	  SP pt1=random_sp(cp1);
	  CGAL::Comparison_result c1= tr.compare_point_to_equipower_line_c(pt1,
									   K(0),
									   K(1),
									   plane_y);
	  sum(c0); sum(c1);
	  CGAL::Comparison_result cc0= CGAL::compare(cp0[1], eqp.x());
	  CGAL::Comparison_result cc1= CGAL::compare(cp1[2], eqp.y());
	  if (cc0 != c0 || cc1 != c1) {
	    std::cout << "Error : " << x << std::endl;
	    std::cout << c0 << " " << c1 << " " << cc0 << " " << cc1 << std::endl;
	    std::cout << ss[0] << std::endl; //0
	    std::cout << ss[1] << std::endl; //1
	    std::cout << x << " " << eqp << " 1/10000 1" << std::endl; //2
	    std::cout << cp0 << " 1/10000 1" << std::endl;
	    std::cout << cp1 << " 1/10000 1" << std::endl;
	    CGAL_assertion(0);
	  }

	}
      }
    }
  }
  {
    Summary sum;
    std::cout << "\n\nTesting compare_center_to_equipower_line_c" << std::endl;
    for (int j=0; j< 20; ++j) {
      Traits tr;
      std::vector<S> ss;
      do {
	ss.clear();
	ss.push_back(random_sphere(0));
	ss.push_back(random_sphere(0));
	tr= Traits(ss.begin(), ss.end());
      } while (!tr.intersects(K(0), K(1)));
      CGAL::Bbox_3 bb=ss[0].bbox()+ss[1].bbox();
      
      std::cout << "Spheres are " << ss[0] << ": " << ss[1] << std::endl;

      for (unsigned int i=0; i< 20; ++i) {
	S qs=random_sphere(0, bb);// random_point(), random_radius());
	ss.push_back(qs);
	FT x;
	C2 ca, cb;
	do {
	  x= random_coordinate(bb.xmin(), bb.xmax());
	} while (!slice(ss[0], x, ca) || !slice(ss[1], x, cb));

	P2 eqp= equipower_point(ca, cb);
	if (i%2==0) {
	  ss.push_back(S(P(x, eqp.x(), random_coordinate()), random_coordinate()));
	} else {
	  ss.push_back(S(P(x, random_coordinate(), eqp.y()), random_coordinate()));
	}
      }
      
      tr= Traits(ss.begin(), ss.end());
      for (unsigned int i=2; i< ss.size(); ++i) {
	{
	  FT x;
	  C2 ca, cb;
	  
	  do {
	    x= random_coordinate(bb.xmin(), bb.xmax());
	  } while (!slice(ss[0], x, ca) || !slice(ss[1], x, cb));

	  P2 eqp= equipower_point(ca, cb);
	  SP t=random_sp(x);
	  CGAL::Comparison_result c0= tr.compare_center_to_equipower_line_c(K(i),
									    t,
									    K(0),
									    K(1),
									    plane_x);
	  CGAL::Comparison_result c1= tr.compare_center_to_equipower_line_c(K(i),
									    t,
									    K(0),
									    K(1),
									    plane_y);
	  sum(c0); sum(c1);
	  CGAL::Comparison_result cc0= CGAL::compare(ss[i].center()[1], eqp.x());
	  CGAL::Comparison_result cc1= CGAL::compare(ss[i].center()[2], eqp.y());
	  if (cc0 != c0 || cc1 != c1) {
	    std::cout << "Error : " << x << std::endl;
	    std::cout << c0 << " " << c1 << " " << cc0 << " " << cc1 << std::endl;
	    std::cout << ss[0] << std::endl; //0
	    std::cout << ss[1] << std::endl; //1
	    std::cout << x << " " << eqp << " 1/10000 1" << std::endl; //2
	    std::cout << ss[i] << std::endl; //3
	    std::cout << x << " " << ss[i].center().y() << " "
		      << ss[i].center().z() << " 1/10000 1" << std::endl; //4
	    CGAL_assertion(0);
	  }
	}
	if (i%2==1) {
	  FT x=ss[i].center().x();
	  C2 ca,cb;
	  bool ba=slice(ss[0], x, ca);
	  bool bb=slice(ss[1], x, cb);
	  CGAL_assertion(ba && bb);
	  P2 eqp= equipower_point(ca, cb);
	  SP t=random_sp(x);
	  CGAL::Comparison_result c0= tr.compare_center_to_equipower_line_c(K(i),
									    t,
									    K(0),
									    K(1),
									    plane_x);
	  CGAL::Comparison_result c1= tr.compare_center_to_equipower_line_c(K(i),
									    t,
									    K(0),
									    K(1),
									    plane_y);
	  sum(c0); sum(c1);
	  CGAL::Comparison_result cc0= CGAL::compare(ss[i].center()[1], eqp.x());
	  CGAL::Comparison_result cc1= CGAL::compare(ss[i].center()[2], eqp.y());
	  if (cc0 != c0 || cc1 != c1) {
	    std::cout << "Error : " << x << std::endl;
	    std::cout << c0 << " " << c1 << " " << cc0 << " " << cc1 << std::endl;
	    std::cout << ss[0] << std::endl; //0
	    std::cout << ss[1] << std::endl; //1
	    std::cout << x << " " << eqp << " 1/10000 1" << std::endl; //2
	    std::cout << ss[i] << std::endl; //3
	    std::cout << x << " " << ss[i].center().y() << " "
		      << ss[i].center().z() << " 1/10000 1" << std::endl; //4
	    CGAL_assertion(0);
	  }
	}
      }
    }
  }
 
  {
    Summary sum;
    Traits tr;
    std::vector<std::pair<P, SP> > pts;
    for (unsigned int i=0; i< 50; ++i) {
      P p= random_point();
      pts.push_back(std::make_pair(p, random_sp(p)));
      pts.push_back(std::make_pair(p, random_sp(p)));
    }

    std::cout << "\n\nTesting compare_points_c" << std::endl;
    for (unsigned int i=0; i< pts.size(); ++i) {
      for (unsigned int j=0; j<i; ++j) {
	for (unsigned int k=0; k<3; ++k) {
	  CGAL_AOS3_INTERNAL_NS::Coordinate_index ci(k);
	  CGAL::Comparison_result cr=tr.compare_points_c(pts[i].second, pts[j].second, ci);
	  CGAL_assertion(cr
			 == CGAL::compare(pts[i].first[ci.index()],
					  pts[j].first[ci.index()]));
	  sum(cr);
	}
      }
    }
  }
 
  {
    Summary sum;
    std::vector<S> spheres;
    P cp= random_point();
    spheres.push_back(S(cp, 2));
    Traits tr(spheres.begin(), spheres.end());
    {
      std::cout << "\n\nTesting compare_point_to_rule_c" << std::endl;
      for (unsigned int i=0; i< 100; ++i) {
	P p=random_point(tr.bbox_3());
	SP sp= random_sp(p);
	for (unsigned int k=0; k<3; ++k) {
	  CI ci(k);
	  CGAL::Comparison_result cr=tr.compare_point_to_rule_c(sp, K(0), ci);
	  sum(cr);
	  CGAL_assertion(cr
			 == CGAL::compare(p[ci.index()], cp[ci.index()]));
	}
      }
      for (unsigned int i=0; i< 100; ++i) {
	FT pv[3]={random_coordinate(), random_coordinate(), random_coordinate()};
	pv[i%3]= spheres[0].center()[i%3];
	P p(pv[0], pv[1], pv[2]);
	SP sp= random_sp(p);
	for (unsigned int k=0; k<3; ++k) {
	  CI ci(k);
	  CGAL::Comparison_result cr= tr.compare_point_to_rule_c(sp, K(0), ci);
	  CGAL_assertion(cr
			 == CGAL::compare(p[ci.index()], cp[ci.index()]));
	  sum(cr);
	}
      }
      //SP m0= random_sp(P(cp[1],random_number.get_double(-10,10),random_number.get_double(-10,10)));
      SP m1= random_sp(P(random_coordinate(),cp[1],random_coordinate()));
      SP m2= random_sp(P(random_coordinate(),random_coordinate(),cp[2]));
      K k(0);
      CGAL_assertion(tr.compare_point_to_rule_c(m1, k, plane_x)== CGAL::EQUAL);
      CGAL_assertion(tr.compare_point_to_rule_c(m2, k, plane_y)== CGAL::EQUAL);    
      sum(CGAL::EQUAL); sum(CGAL::EQUAL);
    }
  }
  





  return 0;
}
