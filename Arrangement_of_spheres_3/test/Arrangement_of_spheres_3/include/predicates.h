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
typedef Traits::Center_point_3 CP;
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
using CGAL_AOS3_INTERNAL_NS::Rule_direction;
using CGAL_AOS3_INTERNAL_NS::sweep_coordinate;
using CGAL_AOS3_INTERNAL_NS::other_plane_coordinate;
CI plane_x= plane_coordinate(0);
CI plane_y= plane_coordinate(1);
CI sweep= sweep_coordinate();


std::vector<V> rational_points;
std::vector<V2> rational_points_2;
CGAL::Random_points_in_cube_3<V> random_point_gen(1);
CGAL::Random_points_in_sphere_3<V> random_point_in_sphere_gen(1);
CGAL::Random random_number_gen;


P random_point(CGAL::Bbox_3 bb=CGAL::Bbox_3(0,0,0,1,1,1)) {
  V rp=.5*(*(++random_point_gen) + V(1,1,1));
  return P(rp[0]*FT(bb.xmax()-bb.xmin())+FT(bb.xmin()),
	   rp[1]*FT(bb.ymax()-bb.ymin())+FT(bb.ymin()),
	   rp[2]*FT(bb.zmax()-bb.zmin())+FT(bb.zmin()));
}

P random_point(S s) {
  V rp=*random_point_in_sphere_gen; ++random_point_in_sphere_gen;
  FT rub= std::sqrt(CGAL::to_interval(s.squared_radius()).first);
  return P(s.center().x()+rp.x()*rub,
	   s.center().y()+rp.y()*rub,
	   s.center().z()+rp.z()*rub);
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

P2 random_point_on_circle(C2 c) {
  CGAL_assertion(c.squared_radius() ==1);
    int ri= random_number_gen.get_int(0, rational_points_2.size());
    return c.center()+ rational_points_2[ri];

}

V2 random_vector_on_unit_circle() {
  int ri= random_number_gen.get_int(0, rational_points_2.size());
  return rational_points_2[ri];

}

S random_sphere(CGAL::Bbox_3 bb=CGAL::Bbox_3(-1,-1,-1,1,1,1)) {
  FT rad= random_radius();
  return S(P(random_coordinate(bb.xmin(), bb.xmax()),
	     random_coordinate(bb.ymin(), bb.ymax()),
	     random_coordinate(bb.zmin(), bb.zmax())), rad*rad);
}

P random_point_on_sphere(S s) {
  CGAL_assertion(s.squared_radius()==1);
  int ri= random_number_gen.get_int(0, rational_points.size());
  return s.center()+ rational_points[ri];
}


void initialize_rational_points() {
  {
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
  {
 std::ifstream rpf("data/rational_points_2");
 
  while (true) {
    char buf[1000];
    rpf.getline(buf,1000);
    if (!rpf) break;
    std::istringstream iss(buf);
    FT c[3];
    for (unsigned int i=0; i< 3; ++i) {
      int t;
      iss >> t;
      c[i]= CGAL::Gmpq(t,1);
    }
    c[0]/=c[2];
    c[1]/=c[2];
    CGAL_assertion(CGAL::square(c[0])+CGAL::square(c[1])==CGAL::Gmpq(1,1));
    rational_points_2.push_back(V2(c[0], c[1]));
  }
  std::cout << "Loaded " << rational_points.size() << " point 2s" << std::endl;
  }
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
  CGAL_assertion(ret.compare_c(pt, plane_x)==CGAL::EQUAL);
  CGAL_assertion(ret.compare_c(pt, plane_y)==CGAL::EQUAL);
  CGAL_assertion(ret.compare_c(pt, sweep)==CGAL::EQUAL);
  return ret;
}


SP random_sp(FT x) {
  return random_sp(P(x, random_coordinate(), random_coordinate()));
}

SP random_sp(S s) {
  V lv=random_vector();
  P pt=random_point(s);
  L line(pt, lv);
  SP ret(s, line);				
  CGAL_assertion(ret.is_valid());
  return ret;
}


std::pair<S,S> random_sr_point(const P&pt, CI rule_c, bool smaller) {
  S sr;
  do {
    FT r=random_radius();
    FT sc[3]={random_coordinate(CGAL::to_interval(pt.x()-r).first,
				CGAL::to_interval(pt.x()+r).second), 
	      random_coordinate(),
	      random_coordinate()};
    sc[rule_c.index()]=pt[rule_c.index()];
    sr=S(P(sc[0], sc[1], sc[2]), r);
  } while (sr.bounded_side(pt) != CGAL::ON_UNBOUNDED_SIDE);

  S ss;
  bool this_smaller;
  do {
    ss= random_sphere_through_point(pt);
    this_smaller= (pt[other_plane_coordinate(rule_c).index()] 
		   < ss.center()[other_plane_coordinate(rule_c).index()]);
  } while (this_smaller != smaller);
  std::cout << "For SR point " << pt << " spheres are " << sr << " and " << ss << std::endl;
  return std::make_pair(sr, ss);
}


std::pair<S,S> random_ss_point(const P&pt) {
  S sa, sb;
  //bool ok;
  Traits tr;
  do {
    sa= random_sphere_through_point(pt);
    sb= random_sphere_through_point(pt);
    /*P2 sac2(sa.center()[plane_x.index()],
	    sa.center()[plane_y.index()]);
    P2 sbc2(sb.center()[plane_x.index()],
	    sb.center()[plane_y.index()]);
    P2 pt2(pt[plane_x.index()],
	   pt[plane_y.index()]);
    L2 l(sac2, sbc2);
    ok=l.has_on_negative_side(pt2);*/
    std::vector<S> ss; ss.push_back(sa); ss.push_back(sb);
    tr= Traits(ss.begin(), ss.end());
  } while(tr.debug_separating_plane(K(0), K(1)).oriented_side(pt) == CGAL::NEGATIVE);
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
    std::cout << a << std::endl << b << std::endl;
    CGAL_assertion(0);
    return P2();
  }
}
FT squared_depth(P p, FT x){
  return CGAL::square(p[0]-x);
}
/*bool slice(S s, FT x, C2 &c) {
  FT r2= s.squared_radius() - squared_depth(s.center(), x);
  if (r2 < 0) return false;
  CGAL_assertion(r2>=0);
  c= C2(P2(s.center()[1],
  s.center()[2]), r2);
  return true;
  }*/

CGAL::Bounded_side check_bounded_side(C2 c, P p, CI ci) {
  FT d2= CGAL::square(c.center()[CGAL_AOS3_INTERNAL_NS::project(ci)]- p[ci.index()]);
  if (d2 < c.squared_radius()) return CGAL::ON_BOUNDED_SIDE;
  else if (d2 == c.squared_radius()) return CGAL::ON_BOUNDARY;
  else return CGAL::ON_UNBOUNDED_SIDE;
}




struct Summary {
  Summary():m1_(0), z_(0), p1_(0), o_(0){}
  void operator()(int i) {
    switch (i) {
    case -1: 
      ++m1_; break;
    case 0:
      ++z_; break;
    case +1: 
      ++p1_; break;
    default:
      ++o_;
    }
  }
  ~Summary() {
    std::cout << "Summary " << m1_ << " " << z_ << " " << p1_ << " " << o_ << std::endl;
  }

  int m1_, z_, p1_, o_;
};

