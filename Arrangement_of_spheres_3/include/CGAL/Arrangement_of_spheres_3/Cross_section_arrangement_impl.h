/*#ifdef CGAL_CHECK_EXPENSIVE
#undef CGAL_CHECK_EXPENSIVE
#endif
#ifdef CGAL_CHECK_EXPENSIVE
#undef CGAL_CHECK_EXACTNESS
#endif*/

#include <CGAL/Arrangement_of_spheres_3_basic.h>

#include <CGAL/Arrangement_of_spheres_3/Cross_section_arrangement.h>

#include <CGAL/Arrangement_of_spheres_3/Filtered_sphere_line_intersection.h>
#include <CGAL/IO/Qt_examiner_viewer_2.h>
#include <CGAL/Arrangement_of_spheres_3/Rule_direction.h>



CGAL_AOS3_BEGIN_INTERNAL_NAMESPACE

template <class K>
struct Center: public boost::static_visitor<> {
  typedef typename K::Point_2 result_type;
  result_type operator()(const typename K::Circular_arc_2 &k) const {
    return k.supporting_circle().center();
  }
  result_type operator()(const typename K::Line_arc_2 &k) const {
    CGAL_assertion(0);
    return k.supporting_line().point();
  }
};

CGAL_AOS3_TEMPLATE
template <class C>
struct Cross_section_arrangement CGAL_AOS3_TARG::Rule{
  typedef Filtered_sphere_line_intersection<Circular_k,C> SLI;
  static Sphere_3 unproject(Circle_2 c){
    return Sphere_3(unproject(c.center()), c.squared_radius());
  }
  static Line_3 unproject(Line_2 c){
    return Line_3(unproject(c.point()), unproject(c.to_vector()));
  }
  static Point_3 unproject(Point_2 c){
    NT pt[3]={666,666,666};
    pt[plane_coordinate(0).index()]= c.x();
    pt[plane_coordinate(1).index()]= c.y();
    return Point_3(pt[0], pt[1], pt[2]);
  }
  static Vector_3 unproject(Vector_2 c){
    NT pt[3]={0,0,0};
    pt[plane_coordinate(0).index()]= c.x();
    pt[plane_coordinate(1).index()]= c.y();
    return Vector_3(pt[0], pt[1], pt[2]);
  }
  static Point_2 project(Point_3 c) {
    CGAL_assertion(c[Sweep_coordinate::index()]== 666);
    return Point_2(c[Plane_coordinate_0::index()], c[Plane_coordinate_1::index()]);
  }
  static Vector_2 project(Vector_3 c){
    CGAL_assertion(c[Sweep_coordinate::index()]== 0);
    return Vector_2(c[Plane_coordinate_0::index()], c[Plane_coordinate_1::index()]);
  }
  Rule(Circle_2 s, Curve f, NT inf): f_(f) {
    CGAL_precondition(f.is_rule());
    Vector_3 v= C::template basis_3<Vector_3>();
    if (!f.is_negative()) v=-v;
   
    /*if (C== 0) {
      vc[plane_coordinate(0).index()]= v;
    } else {
      vc[plane_coordinate(1).index()]= v;
      }*/
    //NT pt[3]=unproject(s.center()); //{0,0,0};
    //pt[plane_coordinate(0).index()]=s.center().x();
    //pt[plane_coordinate(1).index()]=s.center().y();
    Point_3 ppt= unproject(s.center()); //(pt[0], pt[1], pt[2]);
    //std::cout << f_ << std::endl;
    //std::cout << "ppt is " << ppt << std::endl;

    Line_3 l= Line_3(ppt, v);

    //std::cout << "l is " << l << std::endl;
   

    Sphere_3 s3(ppt, 
		s.squared_radius());
    start_= SLI(s3,l);
    //std::cout << "start_ is " << start_ << std::endl;

    NT pc[3]={666,666,666};
    pc[C::Other_plane_coordinate::index()] =s.center()[C::Other_plane_coordinate::plane_index()];
    pc[C::index()]= (f.is_negative())? -inf: inf; // pc[C::index()]= -pc[C::index()];

    Point_3 ppc(pc[0], pc[1], pc[2]);
    //std::cout << "ppc is " << ppc << std::endl;
   
    end_= SLI(ppc, l);
    //std::cout << "Type: " << f_ << " from " << start_ << " to " << end_ << std::endl;
  }

  bool is_on(const Point_2 &p) const {
    CGAL_precondition(p[C::Other_plane_coordinate::plane_index()] == constant_coordinate());// return false;
    SLI pc(unproject(p), start_.line());
    CGAL_assertion(pc != SLI());
    return is_on(pc);
  }

  bool is_on(const SLI &pc) const {
    if (f_.is_negative()) {
      return  end_.compare(pc, C::object())== CGAL::SMALLER && start_.compare(pc, C::object()) == CGAL::LARGER;
    } else {
      return end_.compare(pc, C::object())== CGAL::LARGER && start_.compare(pc, C::object()) == CGAL::SMALLER;
    }
  }

  NT constant_coordinate() const {
    return line().point()[C::Other_plane_coordinate::plane_index()];
  }

  typename SLI::Quadratic_NT source_coordinate() const {
    return start_.exact_coordinate(C::object());
  }

  typename SLI::Quadratic_NT target_coordinate() const {
    return end_.exact_coordinate(C::object());
  }

  void clip(SLI o) {
    std::cout << "Clipping " << f_ << " " << start_ << " to " << end_ 
      << " against " << o << std::endl;
    if (!o.is_valid()) return;
    if (is_on(o)) end_=o;

  }
  Line_2 line() const {
    //std::cout << start_.line() << std::endl;
    //std::cout << end_.line() << std::endl;
    CGAL_assertion(start_.line() == end_.line()  || start_.line().opposite() == end_.line());
    return Line_2(project(start_.line().point()),
		  project(start_.line().to_vector()));
  }
  Curve curve() const {
    return f_;
  }

  void clip(Circle_2 c) {
    //std::cout << "Clipping " << f_ << " against " << c << std::endl;
    SLI i0(unproject(c), start_.line());
    SLI i1(unproject(c), start_.line().opposite());
    clip(i0);
    clip(i1);
  }
  
  void clip(Rule<typename C::Other_plane_coordinate> ra) {
    /*std::cout << "Clipping " << f_ << " from " << start_ << " to " 
	      << end_ << " against " << ra.f_ << " from " << ra.start_ 
	      << " to " << ra.end_ << std::endl;*/
    //CGAL_precondition(CA!=CB);
    NT pc[3]={666,666, 666};
    pc[C::index()]= ra.constant_coordinate();
    pc[C::Other_plane_coordinate::index()]= constant_coordinate();
    Point_2 pi(pc[plane_coordinate(0).index()], pc[plane_coordinate(1).index()]);
    //if (CGAL::assign(pi, oi)){
    SLI ib(unproject(pi), start_.line()); 
    clip(ib);
  }

  SLI start_, end_;
  Curve f_;
};


CGAL_AOS3_TEMPLATE
void Cross_section_arrangement CGAL_AOS3_TARG::build_arrangement(const std::vector<Circular_k::Circle_2> &circles, 
					  const std::vector<int> &names){

  Circular_k::Point_2 lb(-inf_, -inf_);
  Circular_k::Point_2 lt(-inf_, inf_);
  Circular_k::Point_2 rt(inf_, inf_);
  Circular_k::Point_2 rb(inf_, -inf_);
  Circular_k::Line_2 l(lb, Circular_k::Vector_2(0,1));
  Circular_k::Line_2 r(rt, Circular_k::Vector_2(0,1));
  Circular_k::Line_2 b(lb, Circular_k::Vector_2(1,0));
  Circular_k::Line_2 t(rt, Circular_k::Vector_2(1,0));
  audit();
  insert(Circular_k::Line_arc_2(l, lb, lt), Curve(Curve::Key::bl_key(), Curve::T_RULE));
  insert(Circular_k::Line_arc_2(r, rb, rt), Curve(Curve::Key::tr_key(), Curve::B_RULE));
  insert(Circular_k::Line_arc_2(b, lb, rb), Curve(Curve::Key::bl_key(), Curve::R_RULE));
  insert(Circular_k::Line_arc_2(t, lt, rt), Curve(Curve::Key::tr_key(), Curve::L_RULE));



  std::vector<Rule<Plane_coordinate_0> > horizontal_rules;
  std::vector<Rule<Plane_coordinate_1> > vertical_rules;

  //std::vector<std::vector<Line_3> > lines;
  //std::cout << "Constructing circles..." << std::flush;
  for (unsigned int i=0; i< circles.size(); ++i){
    horizontal_rules.push_back(Rule<Plane_coordinate_0>(circles[i], 
				       Curve(names[i], Curve::L_RULE), inf_));
    horizontal_rules.push_back(Rule<Plane_coordinate_0>(circles[i],
				       Curve(names[i], Curve::R_RULE), inf_));
    vertical_rules.push_back(Rule<Plane_coordinate_1>(circles[i],
				     Curve(names[i], Curve::B_RULE), inf_));
    vertical_rules.push_back(Rule<Plane_coordinate_1>(circles[i], 
				     Curve(names[i], Curve::T_RULE), inf_));
  }
  std::cout << "done." << std::endl;
  std::cout << "Intersecting segments with circles..." << std::flush;
  for (unsigned int j=0; j < circles.size(); ++j){
    for (unsigned int i=0; i< horizontal_rules.size(); ++i){
      horizontal_rules[i].clip(circles[j]);
    }
    for (unsigned int i=0; i< vertical_rules.size(); ++i){
      vertical_rules[i].clip(circles[j]);
    }
  }
  std::cout << "done." << std::endl;
    
  std::cout << "Intersecting segments with segments..." << std::flush;
  for (unsigned int i=0; i< horizontal_rules.size(); ++i){
    for (unsigned int j=0; j< i; ++j){
      vertical_rules[i].clip(horizontal_rules[j]);
      horizontal_rules[i].clip(vertical_rules[j]);
    }
  }
  std::cout << "done." << std::endl;
  
  std::cout << "Building arrangement (circles)." << std::flush;
   
  typedef std::vector< Arc> ArcContainer;
    
    
  for (unsigned int i=0; i< circles.size(); ++i){
    Circular_k::Root_of_2 r= CGAL::make_root_of_2(NT(1),NT(0),
						  -circles[i].squared_radius(),
						  false);
    //std::cout << CGAL::to_double(r) << " " << circles_[i].squared_radius() << std::endl;
    CGAL_assertion(CGAL::square(r) == Circular_k::Root_of_2(circles[i].squared_radius()));
    CGAL_assertion(r>0);
    NT x(circles[i].center().x());
    NT y(circles[i].center().y());
    Circular_k::Root_of_2 xr(circles[i].center().x());
    Circular_k::Root_of_2 yr(circles[i].center().y());
    typedef Circular_k::Root_for_circles_2_2 RFC;
    RFC tr(x, y+r);
    Circular_k::Circular_arc_point_2 tp(tr);
    Circular_k::Circular_arc_point_2 bp(RFC(x, circles[i].center().y()-r));
    Circular_k::Circular_arc_point_2 lp(RFC(circles[i].center().x()-r, y));
    Circular_k::Circular_arc_point_2 rp(RFC(circles[i].center().x()+r, y));
    /*std::cout << circles_[i] << std::endl;
      std::cout << tp << std::endl;
      std::cout << bp << std::endl;
      std::cout << lp << std::endl;
      std::cout << rp << std::endl;*/
    Circular_k::Circle_2 c(Point_2(circles[i].center().x(),
				   circles[i].center().y()),
			   circles[i].squared_radius());
    /*Arr::Curve_handle chul=*/ 

    int nm= names[i];
    insert(Circular_k::Circular_arc_2(c, lp, bp), 
	   Curve(nm, Curve::LB_ARC));
    //circle_map[chul]=std::pair<int, Part>(i, UL);
    /*Arr::Curve_handle chur=*/ 
    insert(Circular_k::Circular_arc_2(c, bp, rp), 
	   Curve(nm, Curve::RB_ARC));
    //circle_map[chur]=std::pair<int, Part>(i, UR);
    /*Arr::Curve_handle chll=*/ 
    insert(Circular_k::Circular_arc_2(c, rp, tp),
	   Curve(nm, Curve::RT_ARC));
    //circle_map[chll]=std::pair<int, Part>(i, LL);
    /*Arr::Curve_handle chlr=*/ 
    insert(Circular_k::Circular_arc_2(c, tp, lp),
	   Curve(nm, Curve::LT_ARC));
    //circle_map[chlr]=std::pair<int, Part>(i, LR);
    if (i%10 ==0) std::cout << "." << std::flush;
  }
  std::cout << "done. (segments)..." << std::flush;
  typedef Circular_k::Circular_arc_point_2 CAP;
  typedef Circular_k::Root_for_circles_2_2 RFC;
  
  std::vector<std::pair<Line_arc_2, Curve > > arcs;
  for (unsigned int i=0; i< horizontal_rules.size(); ++i){
    CAP p0(RFC(horizontal_rules[i].source_coordinate(),
	       horizontal_rules[i].constant_coordinate()));
    CAP p1(RFC(horizontal_rules[i].target_coordinate(),
	       horizontal_rules[i].constant_coordinate()));
    Circular_k::Line_2 l(Circular_k::Point_2(0, horizontal_rules[i].constant_coordinate()),
			 Circular_k::Vector_2(1, 0));
    /*if (!l.has_on(p0) || !l.has_on(p1)){
      std::cout << "Failed for i=" << i<< std::endl;
      std::cout << p0 << std::endl;
      std::cout << p1 << std::endl;
      std::cout << l << std::endl;
      }*/
    arcs.push_back(std::pair<Line_arc_2, Curve>(Circular_k::Line_arc_2(l,
								       p0,p1),
						Curve(horizontal_rules[i].curve())));
  }
  for (unsigned int i=0; i< vertical_rules.size(); ++i){
    CAP p0(RFC(vertical_rules[i].constant_coordinate(),
	       vertical_rules[i].source_coordinate()));
    CAP p1(RFC(vertical_rules[i].constant_coordinate(),
	       vertical_rules[i].target_coordinate()));
    /*Circular_k::Line_arc_2 l(0, 1, -vertical_rules_[i].constant_coordinate());*/
    Circular_k::Line_2 l(Circular_k::Point_2(vertical_rules[i].constant_coordinate(), 0),
			 Circular_k::Vector_2(0,1));
    /*if (!l.has_on(p0) || !l.has_on(p1)){
      std::cout << "Failed for i=" << i<< std::endl;
      std::cout << p0 << std::endl;
      std::cout << p1 << std::endl;
      std::cout << l << std::endl;
      }*/
    arcs.push_back(std::pair<Line_arc_2, Curve >(Circular_k::Line_arc_2(l,
									p0,p1),
						 vertical_rules[i].curve()));

    //arcs.push_back(Circular_k::Line_arc_2(l, p0,p1));
  }
  
  std::cout << "done.\n Inserting segments." << std::flush;
  //std::map<Arr::Curve_handle, int, Handle_less> segments_map;
  
  for (unsigned int i=0; i< arcs.size(); ++i){
    //segments_map[arr.insert(arcs[i])]=i;
    insert(arcs[i].first, arcs[i].second);
    if (i%10 ==0) std::cout << "." << std::flush;
  }
  
  std::cout << "done." << std::endl;
  
  std::cout << "The arrangement has " << number_of_faces() << " faces." << std::endl;
  
  /*for (Arr::Vertex_iterator vit= arr.vertices_begin(); vit != arr.vertices_end(); ++vit){
    CAP ip= vit->point();
    Point_2 dp(CGAL::to_double(ip.x()), CGAL::to_double(ip.y()));
    if (CGAL::abs(dp.x()) <10000 && CGAL::abs(dp.y()) <10000) {
    points_.push_back(dp);
    //std::cout << dp << std::endl;
    }
    }*/

  //Slice_data_structure sds;
  //sds.reserve(arr.number_of_vertices(), arr.number_of_halfedges(), arr.number_of_faces()+1);
  //sds.set_building(true);
  /*std::cout << "Traversing the arrangement..." << std::flush;
  int fi=0;
  for (Face_iterator fit = faces_begin(); fit != faces_end(); ++fit){
    std::cout << fi << ": ";
    ++fi;
    for (unsigned int i=0; i< fit->size(); ++i){
      std::cout << fit->at(i).first << " -- " << fit->at(i).second << "--";
    }
    std::cout << std::endl;
    //sds.new_face(fit->begin(), fit->end());
  }
  
  std::cout << "done." << std::endl;*/
 
  //sds.set_building(false);
  
  // sds.audit();
    

  /*std::cout << "There were " << num_tests_ << " comparisons of which " 
    << num_exact_tests_ << " had to be evaluated exactly." << std::endl;*/
  //if (overwrite) sds_=sds;
}


CGAL_AOS3_TEMPLATE
bool Cross_section_arrangement CGAL_AOS3_TARG::Line_arc_less::operator()(Circular_k::Line_arc_2 a, 
									 Circular_k::Line_arc_2 b) const {
  if (a.source().x() < b.source().x()) return true;
  else if (a.source().x() > b.source().x()) return false;
  else if (a.source().y() < b.source().y()) return true;
  else if (a.source().y() > b.source().y()) return false;
  else if (a.supporting_line().point().x() 
	   < b.supporting_line().point().x()) return true;
  else if (a.supporting_line().point().x() 
	   > b.supporting_line().point().x()) return false;
  else if (a.supporting_line().point().y() 
	   < b.supporting_line().point().y()) return true;
  else if (a.supporting_line().point().y() 
	   > b.supporting_line().point().y()) return false;
  else if (a.supporting_line().to_vector().x() 
	   < b.supporting_line().to_vector().x()) return true;
  else if (a.supporting_line().to_vector().x() 
	   > b.supporting_line().to_vector().x()) return false;
  else return (a.supporting_line().to_vector().y() 
	       < b.supporting_line().to_vector().y());
}



CGAL_AOS3_TEMPLATE
bool Cross_section_arrangement CGAL_AOS3_TARG::Circular_arc_less::operator()(Circular_k::Circular_arc_2 a, 
						      Circular_k::Circular_arc_2 b) const {
  if (a.source().x() < b.source().x()) return true;
  else if (a.source().x() > b.source().x()) return false;
  else if (a.source().y() < b.source().y()) return true;
  else if (a.source().y() > b.source().y()) return false;
  else if (a.supporting_circle().center().x() 
	   < b.supporting_circle().center().x()) return true;
  else if (a.supporting_circle().center().x() 
	   > b.supporting_circle().center().x()) return false;
  else if (a.supporting_circle().center().y() 
	   < b.supporting_circle().center().y()) return true;
  else if (a.supporting_circle().center().y() 
	   > b.supporting_circle().center().y()) return false;
  else return (a.supporting_circle().squared_radius() 
	       < b.supporting_circle().squared_radius());
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Cross_section_arrangement CGAL_AOS3_TARG::Vertex_iterator
Cross_section_arrangement CGAL_AOS3_TARG::vertices_begin() {
  return arr_.vertices_begin();
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Cross_section_arrangement CGAL_AOS3_TARG::Vertex_iterator 
Cross_section_arrangement CGAL_AOS3_TARG::vertices_end()  {
  return arr_.vertices_end();
}

CGAL_AOS3_TEMPLATE
void Cross_section_arrangement CGAL_AOS3_TARG::audit() const {
  /*for (CArr::Halfedge_const_iterator hit = arr_.halfedges_begin(); hit != arr_.halfedges_end(); ++hit){
    std::cout << "\"" << feature(hit) << "\" ";
    }
    std::cout << std::endl;*/
}


/*Cross_section_arrangement::Curve_handle*/ 
CGAL_AOS3_TEMPLATE
void Cross_section_arrangement CGAL_AOS3_TARG::insert(const Circular_k::Line_arc_2 &la, Curve f){
  //std::cout << "Inserting: ";
  //write(std::cout, la);
  //std::cout << std::endl;
  CArr::Curve_handle h= CGAL::insert_curve(arr_, Arc(la), pl_);
  map_[h]=f;
  audit();
}
/*Arr::Curve_handle Arr::insert(const Circular_k::Circle_2 &la) {
  return CGAL::insert_curve(arr_, Arc(la), pl_);
  }*/
//Arr::Curve_handle
CGAL_AOS3_TEMPLATE
void Cross_section_arrangement CGAL_AOS3_TARG::insert(const Circular_k::Circular_arc_2 &la,
						      Curve f) {
  //std::cout << "Inserting: ";
  //write(std::cout, la);
  //std::cout << std::endl;
  CArr::Curve_handle h= CGAL::insert_curve(arr_, Arc(la), pl_);
  map_[h]=f;
  audit();
  //map_[ch]= std::pair<int, Part>(curve, pt);
  //return ch;
}

CGAL_AOS3_TEMPLATE
std::size_t Cross_section_arrangement CGAL_AOS3_TARG::number_of_faces() const {
  return arr_.number_of_faces()-1;
}

CGAL_AOS3_TEMPLATE
std::size_t Cross_section_arrangement CGAL_AOS3_TARG::number_of_halfedges() const {
  return arr_.number_of_halfedges();
}
CGAL_AOS3_TEMPLATE
std::size_t Cross_section_arrangement CGAL_AOS3_TARG::number_of_vertices() const {
  return arr_.number_of_vertices();
}

CGAL_AOS3_TEMPLATE
void Cross_section_arrangement CGAL_AOS3_TARG::write(std::ostream& out, const Circular_k::Circular_arc_2 &k) const{
  out << "c((" << CGAL::to_double(k.supporting_circle().center().x()) 
      << ", " << CGAL::to_double(k.supporting_circle().center().y()) 
      << "), " << CGAL::to_double(k.supporting_circle().squared_radius())
      << ") i((" << CGAL::to_double(k.source().x())
      << ", " << CGAL::to_double(k.source().y())
      << ")...(" << CGAL::to_double(k.target().x())
      << ", " << CGAL::to_double(k.target().y()) << ")";
}
CGAL_AOS3_TEMPLATE
void Cross_section_arrangement CGAL_AOS3_TARG::write(std::ostream& out, const Circular_k::Line_arc_2 &k) const{
  out << "l((" << CGAL::to_double(k.supporting_line().point().x()) 
      << ", " << CGAL::to_double(k.supporting_line().point().y()) 
      << "), v(" << CGAL::to_double(k.supporting_line().to_vector().x())
      << ", " << CGAL::to_double(k.supporting_line().to_vector().y())
      << ")) i((" << CGAL::to_double(k.source().x())
      << ", " << CGAL::to_double(k.source().y())
      << ")...(" << CGAL::to_double(k.target().x())
      << ", " << CGAL::to_double(k.target().y()) << ")";
}

/**/

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Cross_section_arrangement CGAL_AOS3_TARG::Point
Cross_section_arrangement CGAL_AOS3_TARG::point(CArr::Vertex_const_handle h) const{
  //std::cout << "Looking for vertex (" << vmap_.size() << ", " << vset_.size() << ", " << map_.size() << ", " << this << ")..." << std::flush;
  if (vmap_.find(h) != vmap_.end()){
    //std::cout << "found." << std::endl;
    CGAL_assertion(vmap_.find(h)->second.is_valid());
    return vmap_.find(h)->second;
  } else {
    CArr::Halfedge_around_vertex_const_circulator c= h->incident_halfedges();
    Curve fa= curve(c);
    CArr::Halfedge_const_handle ha=c, hb=c;
    CGAL_precondition(fa.is_valid());
    Curve fb=fa;
    
    CGAL_assertion(h->degree() > 1);
    if (h->degree()==2){
      CGAL_assertion(!fa.is_finite());
    }
    ++c;
    do {
      Curve fc= curve(c);
      if ((fc.key() != fa.key() || fa.is_rule() || fc.is_rule()) && fa != fc) {
	fb=fc;
	hb=c;
	break;
      }
      ++c;
    } while (c != h->incident_halfedges());
    CGAL_assertion(fa!=fb);
    CGAL_precondition(fb.is_valid());
    
    Point pv;
    if (fa.is_arc() && fb.is_arc()) { 
      //std::cout << "\nSorting out " << fa << ", " << fb;
      Line_2 l(boost::apply_visitor(Center<Circular_k>(), ha->curve()),
	       boost::apply_visitor(Center<Circular_k>(), hb->curve()));
      Line_2 lo(boost::apply_visitor(Center<Circular_k>(), hb->curve()),
		boost::apply_visitor(Center<Circular_k>(), ha->curve()));
      CGAL_assertion(l.opposite()==lo);
      /*std::cout << " with coordinates " << CGAL::to_double(h->point().x()) << " (" << h->point().x() << "), " 
	<< CGAL::to_double(h->point().y()) << " (" << h->point().y() << ")" << std::endl;*/
      Exact_NT va= l.a()*h->point().x();
      Exact_NT vb= l.b()*h->point().y()+l.c();
      /*std::cout << l << std::endl;
      std::cout << "va is " << CGAL::to_double(va) << "(" << va << ")" << std::endl;
      std::cout << "vb is " << CGAL::to_double(vb) << "(" << vb << ")" << std::endl;*/
      if ((va <0 && vb <0) || (va<0 && CGAL::abs(va) > vb) || (vb<0 && CGAL::abs(vb) > va) || va==vb && fa.key() > fb.key()) {
	std::swap(fa, fb);
      }
      //std::cout << "Final is " << fa << ", " << fb << std::endl;
    } 
    if (fa.key() == fb.key() && fa.key().is_input() && fb.key().is_input()) {
      Rule_direction ri;
      if (fa.is_rule()) ri= fa.rule_direction();
      else ri= fb.rule_direction();
      pv= Point::make_extremum(fa.key(), ri);
    } else {
      pv=Point(fa, fb);
    }
    CGAL_assertion(pv.is_valid());

    //std::cout << "creating " << pv << std::endl;
    vmap_[h]=pv;
    //vset_.insert(pv);
    return pv;
  }

  //return Vertex(std::min(fa,fb), std::max(fa,fb));
}


CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Cross_section_arrangement CGAL_AOS3_TARG::Curve
Cross_section_arrangement CGAL_AOS3_TARG::curve(CArr::Halfedge_const_handle h) const {
  CGAL_assertion(map_.find(arr_.originating_curves_begin(h)) != map_.end());
  CGAL_AOS3_TYPENAME Cross_section_arrangement::Curve c= map_.find(arr_.originating_curves_begin(h))->second;
  if (std::distance(arr_.originating_curves_begin(h),
		    arr_.originating_curves_end(h))!=1) {
    CGAL_assertion(c.is_rule());
  }
  return c;
  //return boost::apply_visitor(lookup_curve_data_, h->curve());
  
  /*if (Circular_k::Circular_arc_2& a=boost::get<Circular_k::Circular_arc_2>(h->curve())) {
    return circle_map_[*a];
    } else if (Circular_k::Linear_arc_2& a=boost::get<Circular_k::Circular_arc_2>(h->curve())) {){
    return linear_map_[*a];
    } else {
    CGAL_precodintion(0);
    }*/ 
}

CGAL_AOS3_TEMPLATE
void Cross_section_arrangement CGAL_AOS3_TARG::Face_iterator::increment() {
  ++c_;
  cache_.clear();
  if (c_ != e_ && c_->is_unbounded()) ++c_;
  if (c_ != e_) {
    fill_cache();
  }
}


CGAL_AOS3_TEMPLATE
Cross_section_arrangement CGAL_AOS3_TARG::Face_iterator::Face_iterator(CArr::Face_const_iterator c, 
						CArr::Face_const_iterator e,
						const Cross_section_arrangement *a): c_(c),
									     e_(e),
									     a_(a){
  if (c_ != e_) {
    if (c_->is_unbounded()) increment();
    else fill_cache();
  }
}


CGAL_AOS3_TEMPLATE
void Cross_section_arrangement CGAL_AOS3_TARG::Face_iterator::fill_cache() {
  CGAL_precondition(cache_.empty());
  CArr::Ccb_halfedge_const_circulator cc= c_->outer_ccb();
  /*if (!cc) {
    std::cout << "Empty face?" << std::endl;
    return;
    }*/
  do {
    CArr::Vertex_const_handle t= cc->target();
    Point v= a_->point(t);
    Curve f= a_->curve(cc);
    //std::cout << "Vertex: " << v.first << " " << v.second << std::endl;
    //std::cout << "Curve: " << f << std::endl;
    CGAL::Comparison_result d= cc->direction();
    /*if (d != CGAL::LARGER && d != CGAL::SMALLER) {
      std::cout << "They lied." << std::endl;
      }*/
    if (f.is_arc()) {
      if (d== CGAL::LARGER && f.is_top()
	  || d== CGAL::SMALLER && !f.is_top()){
	f.set_is_inside(true);
      }
    } else if (f.is_rule()){
      if (f.is_finite()) {
	if (f.is_vertical()) {
	  f.set_is_inside( d== CGAL::SMALLER);
	} else {
	  f.set_is_inside( d== CGAL::LARGER);
	}
      } else {
	if (f.is_vertical()){
	  //CGAL_assertion((f.part()& Curve::L_BIT) || (f.part()& Curve::R_BIT));
	  f.set_is_inside( d== CGAL::SMALLER);
	} else {
	  //CGAL_assertion((f.part()& Curve::T_BIT) || (f.part()& Curve::B_BIT));
	  f.set_is_inside( d== CGAL::LARGER);
	}
      }
      /*if (d== CGAL::LARGER && !f.is_vertical()
	|| d== CGAL::SMALLER && f.is_vertical()) {
	f.set_is_inside(true);
	}*/
			       
    }
    //std::cout << "io: " << f.is_inside() << std::endl;
    cache_.push_back(FP(f, v));
    ++cc;
  } while (cc != c_->outer_ccb());
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Cross_section_arrangement CGAL_AOS3_TARG::Face_iterator 
Cross_section_arrangement CGAL_AOS3_TARG::faces_begin() const {
  return Face_iterator(arr_.faces_begin(), arr_.faces_end(), this);
}

CGAL_AOS3_TEMPLATE
CGAL_AOS3_TYPENAME Cross_section_arrangement CGAL_AOS3_TARG::Face_iterator 
Cross_section_arrangement CGAL_AOS3_TARG::faces_end() const {
  return Face_iterator(arr_.faces_end(), arr_.faces_end(), this);
}

CGAL_AOS3_TEMPLATE
 struct Cross_section_arrangement CGAL_AOS3_TARG::Get_circle {
    typedef Circle_2 result_type;
   
    result_type operator()(Line_arc_2) const {
      CGAL_assertion(0);
      return result_type();
    }
    result_type operator()(Circular_arc_2 c) const {
      return c.supporting_circle();
    }
  };
 

CGAL_AOS3_TEMPLATE
struct Cross_section_arrangement CGAL_AOS3_TARG::Get_segment {
  typedef Line_arc_2 result_type;
  
  result_type operator()(Circular_arc_2) const {
    CGAL_assertion(0);
    return result_type();
  }
  result_type operator()(Line_arc_2 c) const {
    /*typename CK::Point_2 s(CGAL::to_double(c.source().x()), CGAL::to_double(c.source().y()));
    typename CK::Point_2 t(CGAL::to_double(c.target().x()), CGAL::to_double(c.target().y()));
    return result_type(s, t);*/
    return c;
  }
};


CGAL_AOS3_TEMPLATE
void Cross_section_arrangement CGAL_AOS3_TARG::draw(Qt_examiner_viewer_2 *qtv) const {
  typedef Qt_examiner_viewer_2 Q;
   
  std::vector<bool> drawn(nums_, false);
    
  for (CArr::Halfedge_const_iterator cit= arr_.halfedges_begin(); cit != arr_.halfedges_end(); ++cit){
    if (curve(cit).is_arc()){
      if (!drawn[curve(cit).key().input_index()]){
	drawn[curve(cit).key().input_index()]=true;
	qtv->set_updating_box(true);
	*qtv << CGAL::BLACK;
	*qtv <<boost::apply_visitor(Get_circle(), cit->curve());
      }
    } else if (&*cit < &*cit->twin()){
      qtv->set_updating_box(false);
      *qtv << CGAL::GRAY;
      *qtv <<  boost::apply_visitor(Get_segment(), cit->curve());
    }
  }

   
   
    
  for (CArr::Vertex_const_iterator vit= arr_.vertices_begin(); vit != arr_.vertices_end(); ++vit){
    //GK::Point_2 dp(CGAL::to_double(vit->point().x()), CGAL::to_double(vit->point().y()));
    if (CGAL::abs(vit->point().x()) == inf_ || CGAL::abs(vit->point().y())==inf_) {
      qtv->set_updating_box(false);
    } else {
      qtv->set_updating_box(true);
    }
    *qtv << CGAL::GREEN;
    *qtv << vit->point();
    Point pt= point(vit);
    std::ostringstream out;
    if (pt.is_sphere_extremum()) {
      out << pt.sphere_key();
    } else if (pt.is_sphere_rule()) {
      out << pt.sphere_key() << ":" << pt.rule_key();
    } else if (pt.is_rule_rule()) {
      out << pt.rule_key(plane_coordinate(0)) << ":" << pt.rule_key(plane_coordinate(1));
    } else {
      out << pt.sphere_key(0) << ":" << pt.sphere_key(1);
    }
    /*Curve::Key sa= pt.first().key();
    Curve::Key sb= pt.second().key();
   
    if (sa != sb) {
      out << sa << ":" << sb;
    } else {
      out << sa;
      }*/
    *qtv << CGAL::GRAY;
    *qtv << out.str().c_str();
  }
}
CGAL_AOS3_END_INTERNAL_NAMESPACE
