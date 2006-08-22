#ifndef UPDATABLE_DELAUNAY_TRIANGULATION_2_H
#define UPDATABLE_DELAUNAY_TRIANGULATION_2_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Simulation_traits.h>
#include <CGAL/Kinetic/Free_event_base.h>
#include <CGAL/Kinetic/Delaunay_triangulation_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_vertex_base_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_recent_edges_visitor_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_visitor_base_2.h>
#include <CGAL/Kinetic/Active_objects_update_vector.h>
#include <CGAL/Indirect_point_2_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polynomial/Interval_root.h>

#include <vector>
#include <boost/dynamic_bitset.hpp>

CGAL_BEGIN_NAMESPACE


/*
  Where to go next:

  - compute interval bound on root and see if that is good enough 
    (use compare_concurrent if they overlap).

  - maintain a list of certificates which are valid at the next event
     and just check at each successive event to make sure that they
     can be advanced. This requires being able to create certificates
     from the Traits which is an easy change to Delaunay but still a
     change.
*/
template <class IndirectKernel>
struct Updatable_Delaunay_triangulation_2 {

  typedef typename IndirectKernel::Point_2 Point_key;
  
  typedef typename Kinetic::Suggested_exact_simulation_traits_base SimTraits_base;





  struct Simulation_traits: public SimTraits_base {
    typedef SimTraits_base P;
    typedef typename Kinetic::Active_objects_update_vector<P::Kinetic_kernel::Point_2, Point_key> Active_points_2_table;
    Active_points_2_table* active_points_2_table_handle() {
      return ap_.get();
    }
    const Active_points_2_table* active_points_2_table_handle() const {
      return ap_.get();
    }
    
    typedef typename Kinetic::Cartesian_instantaneous_kernel<Active_points_2_table, P::Static_kernel> Instantaneous_kernel;
    Instantaneous_kernel instantaneous_kernel_object() const {
      return Instantaneous_kernel(ap_, P::static_kernel_object());
    }
    
    Simulation_traits(const P::Time &lb,
		      const P::Time &ub): P(lb,ub), 
					  ap_(new Active_points_2_table()){}
    
    
  protected:
    typename Active_points_2_table::Handle ap_;
  };

  //  typedef typename Simulation_traits::Kinetic_kernel::Point_2 Kinetic_point_2;
  typedef CGAL::Delaunay_triangulation_2<IndirectKernel,
					 CGAL::Triangulation_data_structure_2<
    CGAL::Kinetic::Delaunay_triangulation_vertex_base_2<IndirectKernel>,
    CGAL::Kinetic::Delaunay_triangulation_face_base_2<Simulation_traits> > > Triangulation;
  typedef CGAL::Kinetic::internal::Triangulation_data_structure_helper_2<typename Triangulation::Triangulation_data_structure> TDS_helper;
  typedef IndirectKernel Indirect_kernel;

  typedef typename Triangulation::Edge Edge;

  typedef typename Indirect_kernel::Geometric_point_2 Static_point_2;


  typedef typename Indirect_kernel::Swapable_container Points;
  typedef typename Simulation_traits::Kinetic_kernel::Point_2 Kinetic_point_2;
  typedef typename Simulation_traits::Kinetic_kernel::Motion_function Kinetic_coordinate;

  typedef typename Indirect_kernel::Current_coordinates IK_current_coordinates;

  typedef typename Simulation_traits::Simulator::NT NT;
  typedef CGAL::Interval_nt_advanced INT;

  struct Update_information: public Kinetic::Ref_counted<Update_information> {
    typedef boost::dynamic_bitset<> Active;
    Update_information(Indirect_kernel ik, 
		       typename Simulation_traits::Active_points_2_table::Handle aot):aot_(aot),
									     ik_(ik), 
									     active_(ik_.number_of_point_2s(), false),
									     start_time_(ik_.number_of_point_2s(), 1),
									     failure_time_(ik_.number_of_point_2s(), 1){
      clear_stats();
    }

    void clear_stats() {
      static_certificates_=0;
      //unsigned int num_events_=0;
      bad_edges_=0;
      num_interpolations_=0;
      constant_filtered_=0;
      interval_filtered_=0;
      interval2_filtered_=0;
      interval3_filtered_=0;
      kinetic_certificates_=0;
      unfailing_kinetic_certificates_=0;
    }


    INT interval(Point_key a, int i) const {
      if (!is_active(a)) return CGAL::to_interval(initial(a)[i]);
      INT ii= CGAL::to_interval(initial(a)[i]);
      INT fi= CGAL::to_interval(final(a)[i]);
      return INT(std::min(fi.inf(), ii.inf()),
		 std::max(fi.sup(), ii.sup()));
    }


    // need start time for each point
    INT cur_interval(Point_key a, INT ct, int i) const {
      if (is_active(a)) {

	INT ii= CGAL::to_interval(initial(a)[i]);
	INT fi= CGAL::to_interval(final(a)[i]);

	INT c[2];
	// (st-1)fi-ii+fi
	// stfi -ii
	INT stm1= (start_time(a)-1);
	c[1]=(ii-fi)/stm1;
	c[0]= (start_time(a) *fi - ii)/stm1;

	//CGAL_precondition(ct.inf() >= 0);
	//CGAL_precondition(ct.sup() <= 1);

	return c[0]+ct*c[1];
      } else {
	return CGAL::to_interval(initial(a)[i]);
      }
    }


    void write_stats(std::ostream &out) {
      out << "The triangululation had " << bad_edges_ << " bad edges" << std::endl;
      out << "Points were interpolated on " << num_interpolations_ << " occasions" << std::endl;
      out << "There were " << static_certificates_ 
	  << " static computations for "
	  << num_edges_
	  << " edges" << std::endl;
      
      out << "Constant filtered removed " << constant_filtered_ 
	  << " computations and interval removed "
	  << interval_filtered_ << " computations" 
	  << " and the second round of interval filtering remove " << interval2_filtered_
	  << " and the third " << interval3_filtered_
	  << std::endl;
      out << "There remained " << kinetic_certificates_ << " certificates that were "
	  << "computed " << unfailing_kinetic_certificates_ 
	  << " of which could have been filtered" << std::endl;
    }


    bool is_active(Point_key k) const {
      return active_[k.to_index()];
    }


    void activate(Point_key k, Kinetic_point_2 kp, INT t) {
      CGAL_precondition(active_[k.to_index()]==false);
      active_[k.to_index()]=true;
      start_time_[k.to_index()]= t;
      aot_->set(k, kp);
    }

    void set_is_editing(bool tf) {
      aot_->set_is_editing(tf);
    }

    void set_final_kernel(Indirect_kernel &fk){
      fk_=fk;
    }

    void active_clear() {
      active_.reset();
    }
    
    Static_point_2 initial(Point_key pk) const {
      return ik_.current_coordinates_object()(pk);
    }
    Static_point_2 final(Point_key pk) const {
      return fk_.current_coordinates_object()(pk);
    }

    INT start_time(Point_key k) const {
      return start_time_[k.to_index()];
    }

    /*double failure_time(Point_key k) const {
      return start_time_[k.to_index()];
      }*/
    

    template <class CNT>
    CNT  incircle(CNT ax, CNT ay, CNT bx, CNT by, 
		  CNT cx, CNT cy, CNT dx, CNT dy,
		  bool stat=false) {
      if (stat) ++static_certificates_;
      CNT qpx = bx - ax;
      CNT qpy = by - ay;
      CNT rpx = cx - ax;
      CNT rpy = cy - ay;
      CNT tpx = dx - ax;
      CNT tpy = dy - ay;
      CNT det=CGAL::det2x2_by_formula(qpx*tpy - qpy*tpx, tpx*(dx - bx) + tpy*(dy - by),
				      qpx*rpy - qpy*rpx, rpx*(cx - bx) + rpy*(cy - by));
      return det;
    }

    typename Simulation_traits::Active_points_2_table::Handle aot_;
    Indirect_kernel ik_, fk_;
    Active active_;
    std::vector<INT> start_time_;
    std::vector<double> failure_time_;

    mutable unsigned int kinetic_certificates_;
    mutable unsigned int unfailing_kinetic_certificates_;
    mutable unsigned int static_certificates_;
    //unsigned int num_events_=0;
    mutable unsigned int bad_edges_;
    mutable unsigned int num_edges_;
    mutable unsigned int num_interpolations_;
    mutable unsigned int constant_filtered_;
    mutable unsigned int interval_filtered_;
    mutable unsigned int interval2_filtered_;
    mutable unsigned int interval3_filtered_;
  };

 
    
    
  struct Traits: public Kinetic::Delaunay_triangulation_default_traits_2<Simulation_traits, Triangulation> {
    typedef Kinetic::Delaunay_triangulation_default_traits_2<Simulation_traits, Triangulation> P;

    struct Cert_tuple {
      Cert_tuple(Point_key a, Point_key b,
		 Point_key c, Point_key d) {
	k_[0]=a;
	k_[1]=b;
	k_[2]=c;
	k_[3]=d;
	std::sort(k_, k_+4);
      }
      bool operator<(const Cert_tuple o) const {
	for (unsigned int i=0; i< 4; ++i){
	  if (k_[i] < o.k_[i]) return true;
	  else if (k_[i] > o.k_[i]) return false;
	}
	return true;
      }
      bool operator==(const Cert_tuple o) const {
	for (unsigned int i=0; i< 4; ++i){
	  if (k_[i] != o.k_[i]) return false;
	}
	return true;
      }

      Point_key k_[4];
    };

    Traits(Simulation_traits tr, typename Update_information::Handle ui): P(tr), ui_(ui){}
    
    //using typename P::Certificate_pair P::certificate_failure_time(typename P::Edge e, P::Certificate_data d);
    using P::certificate_failure_time;

    typename P::Certificate_pair certificate_failure_time(typename P::Edge e) {
      if (is_hull_edge(e)) {
	return null_pair();
      } else {
	typename P::Point_key ks[4];
	P::edge_points(e, ks);
	if (filter_certificate(ks[0], ks[1], ks[2], ks[3])) {
	  return P::certificate_failure_time(e, certificate_data(ks));
	} else {
	  return null_pair();
	}
      }
    }


    typename P::Certificate_data certificate_data(typename P::Point_key ks[4]) {
      Cert_tuple ct(ks[0], ks[1], ks[2], ks[3]);
      if (exact_certs_.find(ct) == exact_certs_.end()) {
	typename P::Certificate_data  s= P::soc_(P::point(ks[0]), 
						 P::point(ks[1]),
						 P::point(ks[2]), 
						 P::point(ks[3]),
						 P::simulator_handle()->current_time(),
						 P::simulator_handle()->end_time());
	++ui_->kinetic_certificates_;
	exact_certs_[ct]= s;
	/*typename P::Certificate_pair rp= P::certificate_failure_time(e, s);
	if (rp.first >= P::simulator_handle()->end_time()) {
	  ++ui_->unfailing_kinetic_certificates_;
	}
	return rp;*/
	return s;
      } else {
	while (exact_certs_[ct].failure_time() <  P::simulator_handle()->current_time()) {
	  exact_certs_[ct].pop_failure_time();
	}
	return exact_certs_[ct];
	//return P::certificate_failure_time(e, exact_certs_[ct]);
      }
    }

    INT join(INT a, INT b) const {
      return INT(std::min(a.inf(), b.inf()),
		 std::max(a.sup(), b.sup()));
    }

    CGAL::Sign sign_at(Point_key a, Point_key b, 
		       Point_key c, Point_key d,
		       INT ct) const {
      INT det= ui_->incircle(ui_->cur_interval(a,ct,0),
			     ui_->cur_interval(a,ct,1),
			     ui_->cur_interval(b,ct,0),
			     ui_->cur_interval(b,ct,1),
			     ui_->cur_interval(c,ct,0),
			     ui_->cur_interval(c,ct,1),
			     ui_->cur_interval(d,ct,0),
			     ui_->cur_interval(d,ct,1),
			     true);
      
      if (det.sup() < 0) return CGAL::NEGATIVE;
      else if (det.inf() > 0) return CGAL::POSITIVE;
      else return CGAL::ZERO;
    }


    bool filter_certificate(Point_key a, Point_key b, 
			    Point_key c, Point_key d) const {
      if (!ui_->is_active(a) && !ui_->is_active(b)
	  && !ui_->is_active(c) && !ui_->is_active(d)) {
	++ui_->constant_filtered_;
	return false;
      }
      //return true;
      CGAL::Protect_FPU_rounding<true> prot;
      
      /*INT det= ui_->incircle(ui_->interval(a,0), ui_->interval(a,1),
			     ui_->interval(b,0), ui_->interval(b,1),
			     ui_->interval(c,0), ui_->interval(c,1),
			     ui_->interval(d,0), ui_->interval(d,1),
			     true);

      CGAL_assertion(det.sup() >0);
      bool ret= det.inf() < 0;*/
      bool ret=true;
      if (!ret) {
	++ui_->interval_filtered_;
	return false;
      } else {
	INT ct= to_interval(P::simulator_handle()->current_time());
	INT rct(ct.inf(), 1);
	bool ret= (sign_at(a,b,c,d, rct) != CGAL::POSITIVE);
	
	if (!ret) {
	  ++ui_->interval2_filtered_;
	  return ret;
	} else {
	  int NS=100;
	  double ld= ct.inf();
	  double step= (1.0-ld)/NS;
	  INT failures;
	  bool has_failures=false;
	  for (int i=0; i< NS; ++i) {
	    double nld= ld+step;
	    INT ci(ld, nld);
	    ld= nld;
	    CGAL::Sign sn= sign_at(a,b,c,d, ci);
	    if (sn != CGAL::POSITIVE) {
	      if (has_failures){ 
		if (sn == CGAL::ZERO) {
		  failures = join(failures, ci);
		}
	      } else {
		failures=ci;
		has_failures=true;
	      }
	      //return true;
	    }
	  }
	  if (has_failures) {
	    
#ifndef NDEBUG
	    typename P::Certificate_data  s= P::soc_(P::point(a), 
						     P::point(b),
						     P::point(c), 
						     P::point(d),
						     P::simulator_handle()->current_time(),
						     P::simulator_handle()->end_time());

	    std::cout << "Failures found for " << a << ", " << b 
		      << ", " << c << ", " << d << " at " << failures 
		      << " exact is " << s.failure_time() << std::endl;
#endif
	    return true;
	  } else {
	    ++ui_->interval3_filtered_;
	  }
	  
	  return false;
	}
      }
    }

    typename P::Certificate_pair null_pair() const {
      return std::make_pair(P::simulator_handle()->end_time(),
			    typename P::Certificate_data());
    }

    typename Update_information::Handle ui_;
    std::map<Cert_tuple,typename P::Certificate_data> exact_certs_;
  };

  struct Visitor: public CGAL::Kinetic::Delaunay_triangulation_visitor_base_2 {
  
    
    Visitor(Simulation_traits tr,
	    typename Update_information::Handle ui): tr_(tr), ui_(ui) {
    } 

   

    bool test_and_add(Edge e, std::vector<Point_key> &active) const {
      ++ui_->static_certificates_;
      if (!compute_ok(e, ui_->fk_)) {
	add(e.first->vertex(0)->point(), active);
	add(e.first->vertex(1)->point(), active);
	add(e.first->vertex(2)->point(), active);
	add(TDS_helper::mirror_vertex(e)->point(), active);
	return true;
      } else return false;
    }

    template <class TDS> 
    void initialize_events(const TDS &triangulation,
			   Indirect_kernel fk) {
      ui_->num_edges_=triangulation.number_of_edges();
      ui_->set_final_kernel(fk);
      std::vector<Point_key> active;
      for (typename TDS::Edge_iterator it= triangulation.edges_begin(); 
	   it != triangulation.edges_end(); ++it){
	if (test_and_add(*it, active)) {
	  ++ui_->bad_edges_;
	}
      }
      
      Interpolate_event ev(tr_,
			   ui_,
			   ui_->ik_.current_coordinates_object(),
			   ui_->fk_.current_coordinates_object(), 
			   active.begin(), 
			   active.end());
      if (!ev.empty()) {
	INT iat= CGAL::to_interval(ev.time());
	tr_.simulator_handle()->new_event(ev.time(), ev);
	/*for (unsigned int i=0; i< active.size(); ++i){
	  ui_->activate(active[i],iat);
	  }*/
      }
    }


  
    void after_flip(Edge e) {
      //++num_events_;
      // schedule a bulk set event for next rational time
      std::vector<Point_key> active;
      test_and_add(e, active);

      Edge em= TDS_helper::mirror_edge(e);
      test_and_add(Edge(e.first, (e.second+1)%3), active);
      test_and_add(Edge(e.first, (e.second+2)%3), active);
      test_and_add(Edge(em.first, (em.second+1)%3), active);
      test_and_add(Edge(em.first, (em.second+2)%3), active);
      
      Interpolate_event ev(tr_, 
			   ui_,
			   ui_->ik_.current_coordinates_object(),
			   ui_->fk_.current_coordinates_object(), 
			   active.begin(), active.end());
      if (!ev.empty()) {
	tr_.simulator_handle()->new_event(ev.time(), ev);
	//active_.set(active.begin(), active.end());
	++ui_->num_interpolations_;
	/*for (unsigned int i=0; i< active.size(); ++i){
	  ui_->activate(active[i], true);
	  }*/
      }
    }

  

   
    

    bool is_active(Point_key k) const {
      return ui_->is_active(k);
    }

    void add( Point_key k, std::vector<Point_key> &active)  const {
      if (!is_active(k)) {
	active.push_back(k);
      }
    }


    void active_clear() {
      ui_->active_clear();
    }
    
    Static_point_2 initial(Point_key pk) const {
      return ui_->ik_.current_coordinates_object()(pk);
    }
    Static_point_2 final(Point_key pk) const {
      return ui_->fk_.current_coordinates_object()(pk);
    }

    void stats_clear() {
      ui_->clear_stats();
    }
    void stats_write(std::ostream &out) {
      ui_->write_stats(out);
    }

    Simulation_traits tr_;
    typename Update_information::Handle ui_;
  
  };



  static Kinetic_point_2 interpolate_t1(NT time,
					Static_point_2 ip,
					Static_point_2 fp) {
    typedef typename Simulation_traits::Kinetic_kernel::Motion_function MF;
    typedef typename MF::NT NT;
    MF mf[2];
    for (unsigned int i=0; i< 2; ++i){
      NT c[2];
      c[1]=(NT(ip[i])-NT(fp[i]))/(time-1);
      c[0]=NT(fp[i])-c[1];
      mf[i]=MF(c, c+2);
    }
    return Kinetic_point_2(mf[0], mf[1]);
  }


  static Kinetic_point_2 interpolate_12(Static_point_2 ip,
					Static_point_2 fp) {
    typedef typename Simulation_traits::Kinetic_kernel::Motion_function MF;
    MF mf[2];
    for (unsigned int i=0; i< 2; ++i){
      NT c[2];
      c[1]=(fp[i]-ip[i]);
      c[0]=2*ip[i]-fp[i];
      mf[i]=MF(c, c+2);
    }
    return Kinetic_point_2(mf[0], mf[1]);
  }



  struct Interpolate_event: public CGAL::Kinetic::Free_event_base {
    //typedef typename Simulation_traits::Active_points_2_table Table;
    typedef typename Simulation_traits::Active_points_2_table::Key Table_key;
    typedef typename std::pair<Table_key, Kinetic_point_2>  MP;

    
  
    
    template <class It>
    Interpolate_event(Simulation_traits tr,
		      typename Update_information::Handle ui,
		      IK_current_coordinates ic,
		      IK_current_coordinates fc, 
		      It b, It e): ui_(ui) {
      time_= tr.simulator_handle()->rational_current_time();
      std::sort(b,e);
      It ne= std::unique(b,e);
      for (It c=b; c!= ne; ++c){
	if (ic(*c) != fc(*c)) {
	  motions_.push_back(MP(*c, interpolate_t1(time_, ic(*c), fc(*c))));
	}
      }
    }
    
    typename Simulation_traits::Simulator::Time time() const {
      return time_;
    }

    std::ostream & write(std::ostream&out) const {
      out << "Updating ";
      for (unsigned int i=0; i< motions_.size(); ++i){
	out << motions_[i].first << " ";
      }
      return out;
    }

    bool empty() const {
      return motions_.empty();
    }

    void process() {
      INT it= CGAL::to_interval(time_);
      ui_->set_is_editing(true);
      for (unsigned int i=0; i< motions_.size(); ++i) {
	/*out << "Setting motion of " << motions_[i].first 
	  << " to " <<  motions_[i].second 
	  << " which is currently " 
	  << motions_[i].second.x()(time_) << " " 
	  << motions_[i].second.y()(time_) 
	  << " with a current position of " 
	  << tr_.active_points_2_table_handle()->at(motions_[i].first) 
	  << std::endl;*/
	ui_->activate(motions_[i].first, motions_[i].second, it);
      }
      ui_->set_is_editing(false);
    }
    typename Update_information::Handle ui_;
    typename std::vector<MP> motions_;
    typename Simulation_traits::Simulator::NT time_;
  };

  typedef CGAL::Kinetic::Delaunay_triangulation_2<Simulation_traits, Visitor, Triangulation, Traits> KDel;
  typedef CGAL::Kinetic::Qt_widget_2<typename Simulation_traits::Simulator> Qt_gui;
  typedef CGAL::Kinetic::Qt_moving_points_2<Simulation_traits, Qt_gui> Qt_mps;
  typedef CGAL::Kinetic::Qt_triangulation_2<KDel, typename Simulation_traits::Instantaneous_kernel, Qt_gui> Qt_triangulation;
    

  struct Final_event: public CGAL::Kinetic::Free_event_base {
    typedef typename Simulation_traits::Active_points_2_table::Key Table_key;

    Final_event(Simulation_traits tr, typename KDel::Handle kdel,
		IK_current_coordinates ic,
		IK_current_coordinates fc):
      tr_(tr), kdel_(kdel),
      ic_(ic), fc_(fc){
    }

    
    std::ostream & write(std::ostream&out) const {
      out << "Final event ";
      return out;
    }

    void process() {
      
      kdel_->write_stats(std::cout);
      kdel_->visitor().stats_write(std::cout);
      
      tr_.simulator_handle()->set_interval(1,2);
      tr_.active_points_2_table_handle()->set_is_editing(true);
      
      for (typename Simulation_traits::Active_points_2_table::Key_iterator 
	     it = tr_.active_points_2_table_handle()->keys_begin(); 
	   it != tr_.active_points_2_table_handle()->keys_end(); ++it) {
	if (!kdel_->visitor().is_active(*it)) {
	  tr_.active_points_2_table_handle()->set(*it, interpolate_12(initial(*it),
								      final(*it)));
	} else {
	  /*std::cout << "Stopping point " << *it 
	    << " at " << fpoints_[it->to_index()]
	    << " ( " << tr_.active_points_2_table_handle()->at(*it).x()(1) 
	    << tr_.active_points_2_table_handle()->at(*it).y()(1)
	    << ")";*/
	  Kinetic_point_2 np(NT(final(*it).x()),
			     NT(final(*it).y()));
	  tr_.active_points_2_table_handle()->set(*it,np);
						  
	}

      }

      tr_.active_points_2_table_handle()->set_is_editing(false);
    }

    Static_point_2 initial(Point_key pk) const {
      return ic_(pk);
    }
    Static_point_2 final(Point_key pk) const {
      return fc_(pk);
    }

    Simulation_traits tr_;
    typename KDel::Handle kdel_;
    IK_current_coordinates ic_, fc_;
  };


 
  


  template <class It> 
  Updatable_Delaunay_triangulation_2(It b, It e): tr_(0,0) {
    
    typename Indirect_kernel::Key_range rg= ik_.new_point_2s(b,e);
    Triangulation tr(ik_);

    tr.insert(rg.first, rg.second);
    typename Update_information::Handle ui= new Update_information(ik_, tr_.active_points_2_table_handle());

    Traits traits(tr_, ui);

    traits.active_points_2_table_handle()->set_is_editing(true);
    for (It c=b; c!= e; ++c){
      traits.active_points_2_table_handle()->insert(Kinetic_point_2(Kinetic_coordinate(c->x()),
								 Kinetic_coordinate(c->y())));
    }
    traits.active_points_2_table_handle()->set_is_editing(false);

  
    
    kdel_= new KDel(traits, tr, Visitor(tr_, ui));
    kdel_->clear_stats();
    kdel_->visitor().stats_clear();
  }
 
  const Triangulation &triangulation() const {
    return kdel_->triangulation(0);
  }

  Triangulation &triangulation() {
    return kdel_->triangulation(0);
  }



  void update_coordinates_demo(const Points &pts) {
    audit();
    kdel_->visitor().stats_clear();
    kdel_->clear_stats();
    Indirect_kernel fk=set_up_update(pts);

    char *argv[1]={"UpdateDel"};
    typename Qt_gui::Handle qtsim= new Qt_gui(1, argv,
					      tr_.simulator_handle());
    
    typename Qt_mps::Handle qtmps= new Qt_mps(qtsim, tr_);
    //qtmps->set_point_size(10);
    typename Qt_triangulation::Handle qtdel
      = new Qt_triangulation(kdel_,
			     tr_.instantaneous_kernel_object(), 
			     qtsim);
    
    tr_.simulator_handle()->new_final_event(Final_event(tr_,
							kdel_,
							ik_.current_coordinates_object(),
							fk.current_coordinates_object()));

    std::cout << "Green edges just flipped, grey edges will not flip until"
	      << " their certificate changes and black edges will flip." << std::endl;
  
    qtsim->begin_event_loop();
    
    ik_.swap(fk);
    audit();
  }

  Indirect_kernel set_up_update(const Points &pts) {
    kdel_->visitor().active_clear();
 tr_.simulator_handle()->set_interval(0,1);
   tr_.active_points_2_table_handle()->set_is_editing(true);
   typename Indirect_kernel::Current_coordinates cc= ik_.current_coordinates_object();
   for (typename Simulation_traits::Active_points_2_table::Key_iterator
	  kit= tr_.active_points_2_table_handle()->keys_begin();
	kit !=  tr_.active_points_2_table_handle()->keys_end(); ++kit) {
     tr_.active_points_2_table_handle()->set(*kit, Kinetic_point_2(Kinetic_coordinate(cc(*kit).x()),
								   Kinetic_coordinate(cc(*kit).y())));
   }
   tr_.active_points_2_table_handle()->set_is_editing(false);


   Indirect_kernel fk;
   fk.new_point_2s(pts.begin(), pts.end());
   
   kdel_->visitor().initialize_events(kdel_->triangulation_data_structure(), fk);
   return fk;
  }

 void update_coordinates(const Points &pts) {
   Indirect_kernel fk=set_up_update(pts);
   tr_.simulator_handle()->set_current_time(1);
   ik_.swap(fk);
 }

  void write_statistics(std::ostream &out) const {
    kdel_->write_stats(std::cout);
    kdel_->visitor().stats_write(std::cout);
  }


  void audit() const {
    typename Indirect_kernel::Current_coordinates cc
      = ik_.current_coordinates_object();
    for (typename KDel::Triangulation::Edge_iterator it
	   = kdel_->triangulation().edges_begin(); 
	 it != kdel_->triangulation().edges_end(); ++it){
      if (!compute_ok(*it, ik_)) {
	std::cout << "Problem with edge " 
		  << it->first->vertex((it->second+1)%3)->point() << " "
		  << it->first->vertex((it->second+2)%3)->point() << " "
		  << cc(it->first->vertex((it->second+1)%3)->point()) << ": "
		  << cc(it->first->vertex((it->second+2)%3)->point())
		  << std::endl;
      }
    }
  }

  static void read(std::string name, Points &points) {
    std::ifstream in(name.c_str());
    if (!in) {
      std::cerr << "Error opening file " << name  << std::endl;
      exit(1);
    }
    
    while (true) {
      char ln[10000];
      in.getline(ln, 10000);
      if (!in) {
	break;
      }
      Static_point_2 p;
      std::istringstream iss(ln);
      iss >> p;
      if (!iss) {
	CGAL::Simple_cartesian<double>::Point_2 dpt;
	std::istringstream iss2(ln);
	iss2 >> dpt;
	if (!iss2) {
	  std::cerr << "Error processing line " << ln << std::endl;
	} else {
	  points.push_back(Static_point_2(dpt.x(), dpt.y()));
	}
      } else {
	points.push_back(p);
      }
    };
  }
  
  static bool is_hull_edge(const Edge &e) {
    return ! TDS_helper::mirror_vertex(e)->point().is_valid()
      || ! TDS_helper::third_vertex(e)->point().is_valid()
      || ! TDS_helper::origin(e)->point().is_valid()
      || ! TDS_helper::destination(e)->point().is_valid();
  }

  
  static bool compute_ok(const Edge &e,  Indirect_kernel sk) {
    //typename Indirect_kernel::Current_coordinates 
    //cc= sk.current_coordinates_object();
    if (is_hull_edge(e)){
      typename Indirect_kernel::Orientation_2 o2= sk.orientation_2_object();
      Point_key ks[4];
      ks[0]= TDS_helper::origin(e)->point();
      ks[1]= TDS_helper::third_vertex(e)->point();
      ks[2]= TDS_helper::destination(e)->point();
      ks[3]= TDS_helper::mirror_vertex(e)->point();
      
      bool odd_parity=false;
      bool infinity=false;
      for (unsigned int i=0; i<4; ++i) {
	if (infinity) {
	  ks[i-1]=ks[i];
	} else {
	  if (ks[i] == Point_key()) {
	    infinity=true;
	    odd_parity= ((i%2)==1);
	  }
	}
      }
      if (odd_parity) {
	std::swap(ks[0], ks[1]);
      }
      CGAL::Orientation o=o2(ks[0], ks[1], ks[2]);
      return o==CGAL::POSITIVE;
    } else {
      typename Indirect_kernel::Side_of_oriented_circle_2 soc
	= sk.side_of_oriented_circle_2_object();

      Point_key ks[4];
      ks[0]= TDS_helper::origin(e)->point();
      ks[1]= TDS_helper::third_vertex(e)->point();
      ks[2]= TDS_helper::destination(e)->point();
      ks[3]= TDS_helper::mirror_vertex(e)->point();
      
      CGAL::Oriented_side s=soc(ks[0], ks[1], ks[2], ks[3]);

      if (s== CGAL::ON_ORIENTED_BOUNDARY) {
	std::cout << "Degeneracy with edge " 
		  << ks[0] << " " << ks[2] << std::endl;
      }
      return s!= CGAL::ON_NEGATIVE_SIDE;
    }
  }

 



  static int run(int argc, char *argv[], int n, int d, int seed, std::string ifile, std::string ffile) {
    //typedef CGAL::Kinetic::Inexact_simulation_traits_2 Traits;
 

    CGAL_KINETIC_SET_LOG_LEVEL(CGAL::Kinetic::LOG_SOME);



    if (ifile.empty() || ffile.empty()) {
      std::cerr << "Need an initial and final coordinate files." << std::endl;
    }

    Points ipoints, fpoints;
    read(ifile, ipoints);
    read(ffile, fpoints);
    CGAL_assertion(ipoints.size() == fpoints.size());

    Updatable_Delaunay_triangulation_2 udt2(ipoints.begin(), ipoints.end());
    udt2.update_coordinates_demo(fpoints);
    udt2.update_coordinates_demo(ipoints);
    return 0;
  }
  Simulation_traits tr_;
  typename KDel::Handle kdel_;
  Indirect_kernel ik_;
};

CGAL_END_NAMESPACE

#endif
