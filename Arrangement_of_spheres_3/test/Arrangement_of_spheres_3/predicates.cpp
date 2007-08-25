#include <CGAL/Arrangement_of_spheres_traits_3.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>


#ifdef CGAL_AOS3_USE_TEMPLATES
typedef CGAL::Cartesian<CGAL::Gmpq> K;
typedef CGAL::Arrangement_of_spheres_traits_3<K> Traits;
#else 
typedef CGAL::Arrangement_of_spheres_traits_3 Traits;
#endif
typedef Traits::Sphere_point_3 SP;
typedef Traits::Event_point_3 EP;
typedef Traits::Point_3 P;
typedef Traits::Line_3 L;
typedef Traits::Point_2 P2;
typedef Traits::Vector_2 V2;
typedef Traits::Circle_2 C2;
typedef Traits::Line_2 L2;
typedef Traits::Vector_3 V;
typedef Traits::Sphere_3 S;
typedef Traits::FT FT;
typedef Traits::Sphere_3_key K;
typedef CGAL_AOS3_INTERNAL_NS::Coordinate_index CI;

using CGAL_AOS3_INTERNAL_NS::plane_coordinate;
using CGAL_AOS3_INTERNAL_NS::sweep_coordinate;
using CGAL_AOS3_INTERNAL_NS::other_plane_coordinate;
CI plane_x= plane_coordinate(0);
CI plane_y= plane_coordinate(1);
CI sweep= sweep_coordinate();


std::vector<V> rational_points;
CGAL::Random_points_in_cube_3<V> random_point_gen(1);
CGAL::Random random_number_gen;


P random_point(CGAL::Bbox_3 bb=CGAL::Bbox_3(0,0,0,1,1,1)) {
  V rp=.5*(*(++random_point_gen) + V(1,1,1));
  return P(rp[0]*FT(bb.xmax()-bb.xmin())+FT(bb.xmin()),
	   rp[1]*FT(bb.ymax()-bb.ymin())+FT(bb.ymin()),
	   rp[2]*FT(bb.zmax()-bb.zmin())+FT(bb.zmin()));
}

V random_vector() {
  return *(++random_point_gen);
}

FT random_radius() {
  return random_number_gen.get_double(.1, .5);
}

FT random_coordinate(double min=0, double max=1) {
  return random_number_gen.get_double(min, max);
}

S random_sphere(FT depth, CGAL::Bbox_3 bb=CGAL::Bbox_3(-1,-1,-1,1,1,1)) {
  FT rad= random_radius();
  FT x= random_number_gen.get_double(CGAL::to_interval(depth-rad).second, 
				     CGAL::to_interval(depth+rad).first);
  return S(P(x, random_coordinate(bb.ymin(), bb.ymax()),
	     random_coordinate(bb.zmin(), bb.zmax())), rad*rad);
}

void initialize_rational_points() {
  std::ifstream rpf("data/rational_points");
 
  while (true) {
    char buf[1000];
    rpf.getline(buf,1000);
    if (!rpf) break;
    std::istringstream iss(buf);
    FT c[4];
    for (unsigned int i=0; i< 4; ++i) {
      int t;
      iss >> t;
      c[i]= CGAL::Gmpq(t,1);
    }
    c[0]/=c[3];
    c[1]/=c[3];
    c[2]/=c[3];
    CGAL_assertion(CGAL::square(c[0])+CGAL::square(c[1])+CGAL::square(c[2])==CGAL::Gmpq(1,1));
    rational_points.push_back(V(c[0], c[1], c[2]));
  }
  std::cout << "Loaded " << rational_points.size() << " points" << std::endl;
}


S random_sphere_through_point(const P&pt, FT radius=random_radius()) {
  int index= random_number_gen.get_int(0, rational_points.size());
  V cv=rational_points[index];
  return S(P(pt+cv*radius), radius*radius);
}

SP random_sp(const P &pt) {
  V lv;
  do {
    lv=random_vector();
  } while (lv == V(0,0,0));
  S sphere= random_sphere_through_point(pt);
  if ((sphere.center()-pt)*lv < 0) lv =-lv;
  L line(pt, lv);
  SP ret(sphere, line);
  CGAL_assertion(ret.compare(pt, plane_x)==CGAL::EQUAL);
  CGAL_assertion(ret.compare(pt, plane_y)==CGAL::EQUAL);
  CGAL_assertion(ret.compare(pt, sweep)==CGAL::EQUAL);
  return ret;
}


SP random_sp(FT x) {
  return random_sp(P(x, random_coordinate(), random_coordinate()));
}

std::pair<S,S> random_sr_point(const P&pt, CI rule_c) {
  FT sc[3]={random_coordinate(), 
	    random_coordinate(),
	    random_coordinate()};
  sc[rule_c.index()]=pt[rule_c.index()];
  S sr(P(sc[0], sc[1], sc[2]), random_radius());
  S ss;
  do {
    ss= random_sphere_through_point(pt);
  } while (CGAL::compare(sr.center()[other_plane_coordinate(rule_c).index()],
			 pt[other_plane_coordinate(rule_c).index()])
	   == CGAL::compare(ss.center()[other_plane_coordinate(rule_c).index()],
			    pt[other_plane_coordinate(rule_c).index()]));
  std::cout << "For SR point " << pt << " spheres are " << sr << " and " << ss << std::endl;
  return std::make_pair(sr, ss);
}


std::pair<S,S> random_ss_point(const P&pt) {
  S sa, sb;
  bool ok;
  do {
    sa= random_sphere_through_point(pt);
    sb= random_sphere_through_point(pt);
    P2 sac2(sa.center()[plane_x.index()],
	    sa.center()[plane_y.index()]);
    P2 sbc2(sb.center()[plane_x.index()],
	    sb.center()[plane_y.index()]);
    P2 pt2(pt[plane_x.index()],
	   pt[plane_y.index()]);
    L2 l(sac2, sbc2);
    ok=l.has_on_negative_side(pt2);
  } while(!ok);
  return std::make_pair(sa,sb);
}

/*SP exact_random_sp(const P &center) {
  CGAL::Random rnd;
  double radius= rnd.get_double(.1,2);
  CGAL::Random_points_on_sphere_3<P> rp(radius);
  V cv(*rp-CGAL::ORIGIN); ++rp;
  V lv=cv;
  S sphere(P(center+cv), cv*cv);
  //if (cv*lv < 0) lv =-lv;
  L line(center, lv);
  SP ret(sphere, line);
  CGAL_postcondition(ret.compare(center, p0) == CGAL::EQUAL);
  CGAL_postcondition(ret.compare(center, p1) == CGAL::EQUAL);
  return ret;
  }*/

bool slice(S s, FT z, C2 &c) {
  CGAL_AOS3_TYPENAME Traits::FT r2= s.squared_radius() -  CGAL::square(s.center()[sweep_coordinate().index()]-z);
  if (r2 >=0) {
    c=  C2(CGAL_AOS3_TYPENAME Traits::Point_2(s.center()[plane_coordinate(0).index()], 
					      s.center()[plane_coordinate(1).index()]), r2);
    return true;
  } else {
    return false;
  }
}

FT discriminant(C2 c) {
  V2 v= c.center()-CGAL::ORIGIN;
  return v*v - c.squared_radius();
}


P2 equipower_point(C2 a, C2 b) {
  V2 n=2*(a.center()-b.center());
  FT d= discriminant(b) - discriminant(a);
  L2 eqp(n[0], n[1], d);

  L2 l(a.center(), (a.center()-b.center()));
  CGAL::Object o= Traits::Geom_traits().intersect_2_object()(eqp, l);
  P2 pt;
  if (CGAL::assign(pt, o)) {
    return pt;
  } else {
    CGAL_assertion(0);
    return P2();
  }
}

int main(int, char *[]) {

  initialize_rational_points();
 
   {
     std::cout << "\n\nTesting compare_point_to_equipower_line_c" << std::endl;
 
     std::vector<S> ss;
     for (unsigned int i=0; i< 100; ++i) {
       //ss.push_back(random
     }

   }

  {
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
    }
  }
 {
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

	for (unsigned int i=0; i< 100; ++i) {
	  S qs=random_sphere(0, bb);// random_point(), random_radius());
	  ss.push_back(qs);
	}
      
	tr= Traits(ss.begin(), ss.end());
	for (unsigned int i=0; i< 20; ++i) {
	  FT x;
	  C2 ca, cb;

	  do {
	    x= random_coordinate(bb.xmin(), bb.xmax());
	  } while (!slice(ss[0], x, ca) || !slice(ss[1], x, cb));

	  P2 eqp= equipower_point(ca, cb);
	  SP t=random_sp(x);
	  CGAL::Comparison_result c0= tr.compare_center_to_equipower_line_c(K(i+2),
									    t,
									    K(0),
									    K(1),
									    plane_x);
	  CGAL::Comparison_result c1= tr.compare_center_to_equipower_line_c(K(i+2),
									    t,
									    K(0),
									    K(1),
									    plane_y);
	  CGAL::Comparison_result cc0= CGAL::compare(ss[i+2].center()[1], eqp.x());
	  CGAL::Comparison_result cc1= CGAL::compare(ss[i+2].center()[2], eqp.y());
	  if (cc0 != c0 || cc1 != c1) {
	    std::cout << "Error : " << x << std::endl;
	    std::cout << c0 << " " << c1 << " " << cc0 << " " << cc1 << std::endl;
	    std::cout << ss[0] << std::endl; //0
	    std::cout << ss[1] << std::endl; //1
	    std::cout << x << " " << eqp << " 1/10000 1" << std::endl; //2
	    std::cout << ss[i+2] << std::endl; //3
	    std::cout << x << " " << ss[i+2].center().y() << " "
		      << ss[i+2].center().z() << " 1/10000 1" << std::endl; //4
	    CGAL_assertion(0);
	  }
	}
      }
    }
  {
    std::cout << "\n\nTesting compare_center_to_circle_circle_c" << std::endl;
    for (int j=0; j< 100; ++j) {
      P pt= random_point();
      SP t= random_sp(pt);
      std::cout << "Point is " << pt << " 1/10000 1" << std::endl;
      std::pair<S,S> rp=random_ss_point(pt);
      std::cout << "Spheres are " << rp.first << ": " << rp.second << std::endl;
      std::vector<S> spheres;
      spheres.push_back(rp.first);
      spheres.push_back(rp.second);
      for (unsigned int i=0; i< 100; ++i) {
	S qs=random_sphere(pt[sweep.index()], rp.first.bbox()+rp.second.bbox());// random_point(), random_radius());
	spheres.push_back(qs);
      }
      Traits tr(spheres.begin(), spheres.end());
      for (unsigned int i=0; i< 100; ++i) {
	CGAL::Comparison_result cr0= tr.compare_center_to_circle_circle_c(K(i+2),
									  t,
									  K(0), K(1),
									  plane_coordinate(0));
	CGAL::Comparison_result cr1= tr.compare_center_to_circle_circle_c(K(i+2),
									  t,
									  K(0), K(1),
									  plane_coordinate(1));
	CGAL::Comparison_result ccr0=CGAL::compare(spheres[i+2].center()[plane_coordinate(0).index()], 
						   pt[plane_coordinate(0).index()]);
	CGAL::Comparison_result ccr1=CGAL::compare(spheres[i+2].center()[plane_coordinate(1).index()], 
						   pt[plane_coordinate(1).index()]);
	if (cr0 != ccr0 || cr1 != ccr1) {
	  std::cout << "Error: " << pt.x() << std::endl;
	  std::cout << cr0 << " " << cr1 << " " << ccr0 << " " << ccr1 << std::endl;
	  std::cout <<  tr.debug_separating_plane(K(0), K(1)).oriented_side(pt) 
		    << " " << tr.debug_separating_plane(K(1), K(0)).oriented_side(pt)  << std::endl;
	  std::cout << spheres[0] << std::endl << std::endl;
	  std::cout << spheres[1] << std::endl << std::endl;
	  std::cout << spheres[i+2] << std::endl << std::endl;
	  std::cout << pt.x() << " " << spheres[i+2].center().y() << " "
		    << spheres[i+2].center().z() << " 1/10000 1"<< std::endl << std::endl;
	  std::cout << "Point is " << pt << " 1/10000 1" << std::endl << std::endl;
	  std::cout << "p " << tr.debug_separating_plane(K(0), K(1)) << std::endl << std::endl;
	  std::cout << "p " << tr.debug_equipower_plane(K(0), K(1)) << std::endl << std::endl;
	  tr.compare_center_to_circle_circle_c(K(i+2),
					       t,
					       K(0), K(1),
					       plane_coordinate(0));
	  tr.compare_center_to_circle_circle_c(K(i+2),
					       t,
					       K(0), K(1),
					       plane_coordinate(1));
	}
	CGAL_assertion(cr0==ccr0);
	CGAL_assertion(cr1==ccr1);
      }
    }
    
   


    std::cout << "\n\nTesting compare_point_to_circle_circle_c" << std::endl;
    for (int j=0; j< 100; ++j) {
      P pt= random_point();
      std::cout << "Point is " << pt << " 1/10000 1" << std::endl;
      SP t= random_sp(pt);
      std::pair<S,S> rp=random_ss_point(pt);
      std::cout << "Spheres are " << rp.first << ": " << rp.second << std::endl;
      std::vector<S> spheres;
      spheres.push_back(rp.first);
      spheres.push_back(rp.second);
      CGAL::Bbox_3 bb= rp.first.bbox()+rp.second.bbox();
      Traits tr(spheres.begin(), spheres.end());
      for (unsigned int i=0; i< 100; ++i) {
	P qpt(pt.x(), random_coordinate(bb.ymin(), bb.ymax()), 
	      random_coordinate(bb.zmin(), bb.zmax()));
	SP qsp= random_sp(qpt);
	CGAL::Comparison_result cr0= tr.compare_point_to_circle_circle_c(qsp,
									 K(0), K(1),
									 
									 plane_coordinate(0));
	CGAL::Comparison_result cr1= tr.compare_point_to_circle_circle_c(qsp,
									 K(0), K(1),
									 
									 plane_coordinate(1));
	CGAL::Comparison_result ccr0=CGAL::compare(qpt[plane_coordinate(0).index()], 
						   pt[plane_coordinate(0).index()]);
	CGAL::Comparison_result ccr1=CGAL::compare(qpt[plane_coordinate(1).index()], 
						   pt[plane_coordinate(1).index()]);
	if (cr0 != ccr0 || cr1 != ccr1) {
	  std::cout << "Error: ";
	  std::cout << qpt << " 1/10000 1" << std::endl;
	  std::cout << cr0 << " " << cr1 << " " << ccr0 << " " << ccr1 << std::endl;
	  tr.compare_point_to_circle_circle_c(qsp,
					      K(0), K(1),
					      plane_coordinate(0));
	  tr.compare_point_to_circle_circle_c(qsp,
					      K(0), K(1),
					      plane_coordinate(1));
	}
	CGAL_assertion(cr0==ccr0);
	CGAL_assertion(cr1==ccr1);
      }
    }
    
  }
  {
    Traits tr;
    std::vector<std::pair<P, SP> > pts;
    for (unsigned int i=0; i< 100; ++i) {
      P p= random_point();
      pts.push_back(std::make_pair(p, random_sp(p)));
    }

    std::cout << "\n\nTesting compare_points_c" << std::endl;
    for (unsigned int i=0; i< pts.size(); ++i) {
      for (unsigned int j=0; j<i; ++j) {
	for (unsigned int k=0; k<3; ++k) {
	  CGAL_AOS3_INTERNAL_NS::Coordinate_index ci(k);
	  CGAL_assertion(tr.compare_points_c(pts[i].second, pts[j].second, ci)
			 == CGAL::compare(pts[i].first[ci.index()],
					  pts[j].first[ci.index()]));
	}
      }
    }
  }
 
  {
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
	  CGAL_assertion(tr.compare_point_to_rule_c(sp, K(0), ci)
			 == CGAL::compare(p[ci.index()], cp[ci.index()]));
	}
      }
      //SP m0= random_sp(P(cp[1],random_number.get_double(-10,10),random_number.get_double(-10,10)));
      SP m1= random_sp(P(random_coordinate(),cp[1],random_coordinate()));
      SP m2= random_sp(P(random_coordinate(),random_coordinate(),cp[2]));
      K k(0);
      CGAL_assertion(tr.compare_point_to_rule_c(m1, k, plane_x)== CGAL::EQUAL);
      CGAL_assertion(tr.compare_point_to_rule_c(m2, k, plane_y)== CGAL::EQUAL);      
    }
  }
  
  {
    std::cout << "\n\nTesting compare_center_to_circle_rule_c" << std::endl;
    for (int j=0; j< 100; ++j) {
      P pt= random_point();
      SP t= random_sp(pt);
      CI rule_coordinate(plane_coordinate(j%2));
      std::pair<S,S> rp=random_sr_point(pt, rule_coordinate);
      std::vector<S> spheres;
      spheres.push_back(rp.first);
      spheres.push_back(rp.second);
      for (unsigned int i=0; i< 100; ++i) {
	S qs=random_sphere(pt.x(), rp.first.bbox()+rp.second.bbox());
	spheres.push_back(qs);
      }
      Traits tr(spheres.begin(), spheres.end());
      for (unsigned int i=0; i< 100; ++i) {
	CGAL::Comparison_result cr0= tr.compare_center_to_circle_rule_c(K(i+2),
									t,
									K(1), K(0),
									rule_coordinate,
									plane_coordinate(0));
	CGAL::Comparison_result cr1= tr.compare_center_to_circle_rule_c(K(i+2),
									t,
									K(1), K(0),
									rule_coordinate,
									plane_coordinate(1));
	CGAL::Comparison_result ccr0=CGAL::compare(spheres[i+2].center()[plane_coordinate(0).index()], 
						   pt[plane_coordinate(0).index()]);
	CGAL::Comparison_result ccr1=CGAL::compare(spheres[i+2].center()[plane_coordinate(1).index()], 
						   pt[plane_coordinate(1).index()]);
	if (cr0 != ccr0 || cr1 != ccr1) {
	  std::cout << "Error: " << j%2 << " ";
	  std::cout << spheres[i+2] << std::endl;
	  tr.compare_center_to_circle_rule_c(K(i+2),
					     t,
					     K(1), K(0),
					     plane_coordinate(i%2),
					     plane_coordinate(0));
	  tr.compare_center_to_circle_rule_c(K(i+2),
					     t,
					     K(1), K(0),
					     plane_coordinate(i%2),
					     plane_coordinate(1));
	}
	CGAL_assertion(cr0==ccr0);
	CGAL_assertion(cr1==ccr1);
      }
    }
    
    std::cout << "\n\nTesting compare_point_to_circle_rule_c" << std::endl;
    for (int j=0; j< 100; ++j) {
      P pt= random_point();
      std::cout << "Point is " << pt << " 1/10000 1" << std::endl;
      SP t= random_sp(pt);
      CI rule_coordinate(plane_coordinate(j%2));
      std::pair<S,S> rp=random_sr_point(pt, rule_coordinate);
      std::vector<S> spheres;
      spheres.push_back(rp.first);
      spheres.push_back(rp.second);
      CGAL::Bbox_3 bb=rp.first.bbox() + rp.second.bbox();
      Traits tr(spheres.begin(), spheres.end());
      for (unsigned int i=0; i< 100; ++i) {
	P qpt(pt.x(), random_coordinate(bb.ymin(), bb.ymax()), random_coordinate(bb.zmin(), bb.zmax()));
	SP qsp= random_sp(qpt);
	CGAL::Comparison_result cr0= tr.compare_point_to_circle_rule_c(qsp,
								       K(1), K(0),
								       rule_coordinate,
								       plane_coordinate(0));
	CGAL::Comparison_result cr1= tr.compare_point_to_circle_rule_c(qsp,
								       K(1), K(0),
								       rule_coordinate,
								       plane_coordinate(1));
	CGAL::Comparison_result ccr0=CGAL::compare(qpt[plane_coordinate(0).index()], 
						   pt[plane_coordinate(0).index()]);
	CGAL::Comparison_result ccr1=CGAL::compare(qpt[plane_coordinate(1).index()], 
						   pt[plane_coordinate(1).index()]);
	if (cr0 != ccr0 || cr1 != ccr1) {
	  std::cout << "Error: " << j%2 << " ";
	  std::cout << qpt << " 1/10000 1" << std::endl;
	  std::cout << cr0 << " " << cr1 << " " << ccr0 << " " << ccr1 << std::endl;
	  tr.compare_point_to_circle_rule_c(qsp,
					    K(1), K(0),
					    plane_coordinate(i%2),
					    plane_coordinate(0));
	  tr.compare_point_to_circle_rule_c(qsp,
					    K(1), K(0),
					    plane_coordinate(i%2),
					    plane_coordinate(1));
	}
	CGAL_assertion(cr0==ccr0);
	CGAL_assertion(cr1==ccr1);
      }
    }
    
  }









  return 0;
}
