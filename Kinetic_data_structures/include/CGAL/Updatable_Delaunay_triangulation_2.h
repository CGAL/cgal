#ifndef UPDATABLE_DELAUNAY_TRIANGULATION_2_H
#define UPDATABLE_DELAUNAY_TRIANGULATION_2_H
#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Interval_simulator_traits.h>
#include <CGAL/Kinetic/Simulation_traits.h>
#include <CGAL/Kinetic/Free_event_base.h>
#include <CGAL/Kinetic/Delaunay_triangulation_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_vertex_base_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_recent_edges_visitor_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_visitor_base_2.h>
#include <CGAL/Kinetic/Active_objects_update_vector.h>
#include <CGAL/Indirect_point_2_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


#include <CGAL/Kinetic/IO/Qt_moving_points_2.h>
#include <CGAL/Kinetic/IO/Qt_triangulation_2.h>
#include <CGAL/Kinetic/IO/Qt_widget_2.h>



#include <vector>
#include <boost/dynamic_bitset.hpp>

#ifdef NDEBUG
#define CGAL_UD_DEBUG(x)
#else
#define CGAL_UD_DEBUG(x) std::cout << x
#endif

CGAL_BEGIN_NAMESPACE

template <class Point_key>
struct Update_cert_tuple {
  Update_cert_tuple(){}
  Update_cert_tuple(Point_key k[4]) {
    std::copy(k, k+4, k_);
    int swaps=0;
    for (unsigned int i=0; i< 3; ++i){
      for (unsigned int j=i+1; j< 4; ++j) {
	if (k_[i] > k_[j]) {
	  std::swap(k_[i], k_[j]);
	  ++swaps;
	}
      }
    }
    
    if (swaps%2 ==1) {
      std::swap(k_[2], k_[3]);
    }
  }

  /*static std::pair<Cert_tuple, bool> make(Point_key k[4]) {
    Cert_tuple r(k);
    int swaps=0;
    for (unsigned int i=0; i< 3; ++i){
    for (unsigned int j=i+1; j< 4; ++j) {
    if (r.k_[i] > r.k_[j]) {
    std::swap(r.k_[i], r.k_[j]);
    ++swaps;
    }
    }
    }
    return std::make_pair(r, swaps%2==0);
    }*/

  Update_cert_tuple canonicalize() const {
    Update_cert_tuple ret= *this;
    if (ret.k_[2] > ret.k_[3]) {
      std::swap(ret.k_[2], ret.k_[3]);
    }
    return ret;
  }
  Update_cert_tuple opposite() const {
    Update_cert_tuple ret= *this;
    std::swap(ret.k_[2], ret.k_[3]);
    return ret;
  }

  bool operator<(const Update_cert_tuple o) const {
    for (unsigned int i=0; i< 4; ++i){
      if (k_[i] < o.k_[i]) return true;
      else if (k_[i] > o.k_[i]) return false;
    }
    return false;
  }
  bool operator==(const Update_cert_tuple o) const {
    for (unsigned int i=0; i< 4; ++i){
      if (k_[i] != o.k_[i]) return false;
    }
    return true;
  }
  Point_key operator[](int i) const {
    return k_[i];
  }
  void write(std::ostream&out) const {
    out << k_[0] << " " << k_[1] << " " << k_[2] << " " << k_[3];
  }
  Point_key k_[4];
};

template <class K>
std::ostream &operator<<(std::ostream &out, 
			 const Update_cert_tuple<K> &ct) {
  out << ct[0] << " " << ct[1] << " " << ct[2] << " " << ct[3];
  return out;
}

/*
  Where to go next:

  - compute interval bound on root and see if that is good enough (use
  compare_concurrent if they overlap). Currently I make sure that
  the event actually occurs as this is easier than trying to handle
  events disappearing. There are some difficulties with new events
  that occur close to the current event. I think if I keep a last
  event around and use its inexact or exact rep to prune any new
  events I look at then I am fine. How to get the rep of the current
  time in the traits is a bit tricky.

  - there should be a way of making sure I have a correct root since
  if there are roots close together and nothing happened in between
  I can take the last (mod 2). So can I mantain the set of
  everything before and everything after? Or do I need a whole
  graph? I know I can disambiguate everything that came before. So
  for a tuple/interval pair either I never computed the exact root,
  in which case I can take the last one with the right sign (can I
  still determine the right sign?) or I computed it, in which case I
  could store it.

  - each tuple can only be in the queue once. If it is there and has
  no exact current root then we can take the first root in the
  interval to be the root. Otherwise, it has an exact root.

  - for advancing, if I don't have an exact root, look for a disjoint
  interval (since I know none of the later roots overlap). Otherwise
  find the next root after the current exact root.

  - cached Certificate means that 


  - maintain a list of certificates which are valid at the next event
  and just check at each successive event to make sure that they can
  be advanced. This requires being able to create certificates from
  the Traits which is an easy change to Delaunay but still a change.
*/
template <class IndirectKernel>
struct Updatable_Delaunay_triangulation_2 {

  typedef typename IndirectKernel::Point_2 Point_key;
  
  typedef typename Kinetic::Suggested_exact_simulation_traits_base SimTraits_base;


  typedef Update_cert_tuple<Point_key> Cert_tuple;




  struct Simulation_traits {
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
      return Instantaneous_kernel(ap_, typename P::Static_kernel());
    }
    typedef typename Kinetic::Interval_simulator_traits<Cert_tuple> Simulator_traits;
    typedef typename Kinetic::Two_list_pointer_event_queue<Simulator_traits, 2> Queue;
    typedef typename Kinetic::Default_simulator<Simulator_traits, Queue > Simulator;

    typename Simulator::Handle simulator_handle() {
      return sim_;
    }


    typename Simulator::Handle simulator_handle() const {
      return sim_;
    }


    typedef typename Simulator::Time Time;


    Simulation_traits(const Time &lb,
		      const Time &ub): ap_(new Active_points_2_table()),
				       sim_(new Simulator(lb, ub)){}

    typedef typename P::Function_kernel Function_kernel;
    Function_kernel function_kernel_object() {
      return Function_kernel();
    }

    typedef typename P::Kinetic_kernel Kinetic_kernel;
    Kinetic_kernel kinetic_kernel_object() const {
      return Kinetic_kernel();
    }

    typedef typename P::Static_kernel Static_kernel;
    Static_kernel static_kernel_object() const {
      return Static_kernel();
    }

    typedef typename P::NT NT;
    
  protected:
    typename Active_points_2_table::Handle ap_;
    typename Simulator::Handle sim_;
  };




  typedef typename Simulation_traits::Simulator::Event_key Event_key;





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
  typedef typename Simulation_traits::Kinetic_kernel::Motion_function::NT ENT;
  typedef CGAL::Interval_nt_advanced INT;

















  struct Update_information: public Kinetic::Ref_counted<Update_information> {
    struct Coef_data {
      Coef_data(INT x_0, INT x_1, INT y_0, INT y_1){
	c_[0][0]=x_0;
	c_[0][1]=x_1;
	c_[1][0]=y_0;
	c_[1][1]=y_1;
      }
      Coef_data(){}
      const INT* operator[](int i) const {
	return c_[i];
      }
      INT c_[2][2];
    };

    typedef boost::dynamic_bitset<> Active;
    Update_information(Indirect_kernel ik, 
		       typename Simulation_traits::Active_points_2_table::Handle aot):aot_(aot),
										      ik_(ik), 
										      active_(ik_.number_of_point_2s(), false),
										      coef_cache_(ik_.number_of_point_2s()),
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
      certificate_advances_=0;
      uncertain_exact_computations_=0;
      exact_current_time_certificates_=0;
      comparison_certificates_=0;
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

	//INT c[2];
	// (st-1)fi-ii+fi
	// stfi -ii

	//CGAL_precondition(ct.inf() >= 0);
	//CGAL_precondition(ct.sup() <= 1);

	return coef_cache_[a.to_index()][i][0] + ct*coef_cache_[a.to_index()][i][1];
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
	  << "computed which were advanced " << certificate_advances_ 
	  << " times and " << unfailing_kinetic_certificates_ 
	  << " of which could have been filtered" << std::endl;
      out << "Of these " << uncertain_exact_computations_ << " were caused "
	  << " by inability to determine if the root was there and " 
	  << exact_current_time_certificates_ 
	  << " were caused to make the current time exact " 
	  << " and " << comparison_certificates_ 
	  << " were caused in order to compare." << std::endl;
    }


    bool is_active(Point_key k) const {
      return active_[k.to_index()];
    }


    void activate(Point_key k, Kinetic_point_2 kp, INT t) {
      CGAL_precondition(t.inf() == t.sup());
      //CGAL_precondition(active_[k.to_index()]==false);
      if (!active_[k.to_index()]) {
	active_[k.to_index()]=true;
	//start_time_[k.to_index()]= t.inf();
	INT stm1=(t-1);

	INT xii= CGAL::to_interval(initial(k)[0]);
	INT xfi= CGAL::to_interval(final(k)[0]);
	INT xc0= (t*xfi - xii)/stm1;
	INT xc1=(xii-xfi)/stm1;

	INT yii= CGAL::to_interval(initial(k)[1]);
	INT yfi= CGAL::to_interval(final(k)[1]);
	INT yc0= (t*yfi - yii)/stm1;
	INT yc1=(yii-yfi)/stm1;
	coef_cache_[k.to_index()]= Coef_data(xc0, xc1, yc0, yc1);
	aot_->set(k, kp);
      }
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

    /*INT start_time(Point_key k) const {
      return CGAL::to_interval(start_time_[k.to_index()]);
      }*/

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
    std::vector<Coef_data > coef_cache_;
    std::vector<double> failure_time_;

    mutable unsigned int kinetic_certificates_;
    mutable unsigned int unfailing_kinetic_certificates_;
    mutable unsigned int certificate_advances_;
    mutable unsigned int static_certificates_;
    //unsigned int num_events_=0;
    mutable unsigned int bad_edges_;
    mutable unsigned int num_edges_;
    mutable unsigned int num_interpolations_;
    mutable unsigned int constant_filtered_;
    mutable unsigned int interval_filtered_;
    mutable unsigned int interval2_filtered_;
    mutable unsigned int interval3_filtered_;
    mutable unsigned int uncertain_exact_computations_;
    mutable unsigned int exact_current_time_certificates_;
    mutable unsigned int comparison_certificates_;
  };

 
    
















  typedef  Kinetic::Delaunay_triangulation_default_traits_2<Simulation_traits, Triangulation> Default_traits;

    
  struct Traits: public Default_traits{
    typedef Default_traits  P;



  
    typedef typename Simulation_traits::Simulator::Time Time;
    struct Certificate_data{};

    typedef typename Simulation_traits::Kinetic_kernel::Positive_side_of_oriented_circle_2::result_type Exact_certificate;
    typedef typename Simulation_traits::Function_kernel::Root Exact_time;

    struct Cache_data {
      Cache_data(){}
      Cache_data(Exact_certificate cert): cert_(cert){}

      const Exact_time& failure_time() const {
	return cert_.failure_time();
      }
  
      void pop_failure_time() {
	//()) {
	cert_.pop_failure_time();
	//}
      }

      Exact_certificate cert_;
    };

    typedef std::map<Cert_tuple,Cache_data> Cache;


    Traits(Simulation_traits tr,
	   typename Update_information::Handle ui): P(tr), ui_(ui), 
						    soc_(tr.kinetic_kernel_object().positive_side_of_oriented_circle_2_object()){}
    
    typedef std::pair<Time, Certificate_data> Certificate_pair;

    Certificate_pair null_pair() const {
      return std::make_pair(Time(2), Certificate_data());
    }
    
    Cert_tuple tuple(typename P::Edge e) const {
      typename P::Point_key ks[4];
      P::edge_points(e, ks);
      return Cert_tuple(ks);
    }


    bool has_exact_failure_time(Cert_tuple ct) const {
      return cache_.find(ct) != cache_.end();
    }

    Certificate_pair return_pair(Time rt) const {
      return Certificate_pair(rt, Certificate_data());
    }

    void point_changed(Point_key k){
      std::vector<typename Cache::iterator> tk;
      for (typename Cache::iterator it = cache_.begin(); it != cache_.end(); ++it) {
	for (unsigned int i=0; i< 4; ++i){
	  if (it->first[i]==k) {
	    tk.push_back(it);
	    break;
	  }
	}
      }
      for (unsigned int i=0; i< tk.size(); ++i) {
	CGAL_UD_DEBUG("Erasing " << tk[i]->first << std::endl);
	cache_.erase(tk[i]);
      }
    }
    
    const Time &current_time() const {
      return P::simulator_handle()->current_time();
    }

    Certificate_pair certificate_failure_time(typename P::Edge e) {
      if (is_hull_edge(e)) return null_pair();

      Cert_tuple ct= tuple(e);
#ifndef NDEBUG
      Exact_time actual_failure_time;
      {
	Cert_tuple cct= tuple(e);
	if (check_.find(cct) != check_.end()) {
	  check_.erase(cct);
	}
	if (check_.find(cct.opposite()) != check_.end()) {
	  check_.erase(cct.opposite());
	}
	Exact_time ect;
	if (current_time().data() == Cert_tuple()) {
	  ect= Exact_time(CGAL::to_interval(current_time()).first);
	} else {
	  CGAL_assertion(check_.find(current_time().data()) != check_.end());
	  ect= check_[current_time().data()].failure_time();
	}
	Exact_certificate ecert=compute_exact_certificate(cct, ect);
	actual_failure_time= ecert.failure_time();
	check_[cct]= ecert;
	
      }
#endif
   
      //std::pair<Cert_tuple, bool> ctp= Cert_tuple::make(ks);
      
      if (!can_fail(ct)) {
	CGAL_postcondition(actual_failure_time > 1);
	return null_pair();
      }
      
      CGAL_UD_DEBUG("Computing failure time for " << ct << std::endl);
      CGAL_UD_DEBUG("Exact failure time is " << actual_failure_time<< std::endl);


    
      


      Time ft;
      bool isc= isolate_failure(ct, 
				CGAL::to_interval(current_time()).first,
				ft);
      if (ft== Time()) { 
	CGAL_UD_DEBUG("Not isolated" << std::endl);
	CGAL_postcondition(actual_failure_time > 1);
	return null_pair();
      }
      if (isc) {
	if (CGAL::compare(ft, current_time()) == CGAL::EQUAL) {
	  // interval overlaps current
	  Exact_time ect= exact_current_time();
	  
	  if (Exact_time(CGAL::to_interval(ft).first) <= ect) {
	    compute_exact_failure_time(ct, ect);
	    CGAL_UD_DEBUG("Exact" << exact_failure_time(ct) <<  std::endl);
	    ft= interval_from_exact_failure_time(ct);
	  }
	}
	CGAL_UD_DEBUG("Returning isolated " << ft << std::endl);
	CGAL_postcondition(actual_failure_time >= CGAL::to_interval(ft).first);
	CGAL_postcondition(actual_failure_time <= CGAL::to_interval(ft).second);
	return return_pair(ft);
      } else {
	//if (!has_exact_failure_time(ct)) {
	++ui_->uncertain_exact_computations_;
	//	}
	Exact_time ect= exact_current_time();
	compute_exact_failure_time(ct, ect);
	ft= interval_from_exact_failure_time(ct);
	CGAL_UD_DEBUG("Not isolated " << exact_failure_time(ct) <<  std::endl);
	if (ft == Time()) {
	  CGAL_UD_DEBUG("Phantom root " << ft << std::endl);
	  CGAL_postcondition(actual_failure_time > 1);
	  return null_pair();
	} else {
	  CGAL_UD_DEBUG("Returning exact " << ft << "(" 
			<< exact_failure_time(ct) << ")" << std::endl);
	  CGAL_postcondition(actual_failure_time >= CGAL::to_interval(ft).first);
	  CGAL_postcondition(actual_failure_time <= CGAL::to_interval(ft).second);
	  return return_pair(ft);
	}
      }
    }


    Certificate_pair certificate_failure_time(typename P::Edge e, 
					      Certificate_data ) {
      Cert_tuple ct= tuple(e);
   
#ifndef NDEBUG
      Exact_time actual_failure_time;
      {
	CGAL_precondition(check_.find(ct.opposite()) != check_.end());
	Exact_certificate ec= check_[ct.opposite()];
	check_.erase(ct.opposite());
	ec.pop_failure_time();
	check_[ct]= ec;
	actual_failure_time= ec.failure_time();
	CGAL_UD_DEBUG("Exact time is " << ec.failure_time() << std::endl);
      }
#endif
      //std::pair<Cert_tuple, bool> ctp= Cert_tuple::make(ks);
      
      if (!can_fail(ct)) { 
	CGAL_assertion(0);
	//std::cout << "Can't fail." << std::endl;
	CGAL_postcondition(actual_failure_time > 1);
	return null_pair();
      }

      CGAL_UD_DEBUG("Advancing failure time for " << ct << std::endl);
      CGAL_UD_DEBUG("Exact time is " << actual_failure_time << std::endl);

      Time ft;
      // this depends on being the last event of the batch whose time is computed
      Time net= P::simulator_handle()->next_event_time();

      bool isc= isolate_failure(ct, CGAL::to_interval(net).first, ft);
      
      if (ft == Time()) { 
	CGAL_UD_DEBUG("Can't isolate." << std::endl);
	CGAL_postcondition(actual_failure_time > 1);
	return null_pair();
      }
      if (isc) {
	if (CGAL::compare(ft, current_time()) == CGAL::EQUAL) {
	  //ensure_exact_failure_time(ct, current_time());
	  // this works because if there is no exact failure, we can take the first in the interval
	  ensure_exact_failure_time(ct, current_time());
	  advance_exact_failure_time(ct);
	  if (exact_failure_time(ct) > 1) {
	    ++ui_->unfailing_kinetic_certificates_;
	  }
	  ft= interval_from_exact_failure_time(ct);
	  CGAL_UD_DEBUG("Separated from current " << exact_failure_time(ct) 
			<< std::endl);
	}
	CGAL_UD_DEBUG("Returning isolated " << ft << std::endl);
	CGAL_postcondition(actual_failure_time >= CGAL::to_interval(ft).first);
	CGAL_postcondition(actual_failure_time <= CGAL::to_interval(ft).second);
	return return_pair(ft);
      } else {
	//ensure_exact_failure_time(ct, current_time());
	//if (!has_exact_failure_time(ct)) {
	++ui_->uncertain_exact_computations_;
	//}
	ensure_exact_failure_time(ct, current_time());
	advance_exact_failure_time(ct);
	if (cache_[ct].failure_time() > 1) {
	  ++ui_->unfailing_kinetic_certificates_;
	}
	ft= interval_from_exact_failure_time(ct);

	if (ft == Time()) {
	  CGAL_UD_DEBUG("Phantom root " << ft << std::endl);
	  CGAL_postcondition(actual_failure_time > 1);
	  return null_pair();
	} else {
	  CGAL_UD_DEBUG("Returning exact " << ft << "(" 
			<< exact_failure_time(ct) << ")" << std::endl);
	  CGAL_postcondition(actual_failure_time >= CGAL::to_interval(ft).first);
	  CGAL_postcondition(actual_failure_time <= CGAL::to_interval(ft).second);
	  return return_pair(ft);
	}
      }
      
    }

  

    void ensure_exact_failure_time(Cert_tuple ct, Time rt) const {
      compute_exact_failure_time(ct, Exact_time(CGAL::to_interval(rt).first));
    }

    void advance_exact_failure_time(Cert_tuple ct) const {
      CGAL_precondition(cache_.find(ct) != cache_.end());
      ++ui_->certificate_advances_;
      cache_[ct].pop_failure_time();
    }
    Time interval_from_exact_failure_time(Cert_tuple ct) const {
      CGAL_precondition(cache_.find(ct) != cache_.end());
      if (cache_[ct].failure_time() >= Exact_time(1)) {
	return Time();
      } else {
	std::pair<double,double> ip= CGAL::to_interval(cache_[ct].failure_time());
	return Time(ip.first, ip.second, ct);
      }
    }

    void compute_exact_failure_time(Cert_tuple ct, Exact_time et) const {
      CGAL_UD_DEBUG("Computing exact time for " << ct 
		    << " from " << et << "(" << cache_.size() << ")" << std::endl);
      if (cache_.find(ct.opposite()) != cache_.end()) {
	Cache_data cd = cache_[ct.opposite()];
	cache_.erase(ct.opposite());
	++ui_->certificate_advances_;
	cd.pop_failure_time();
	cache_[ct]= cd;
	CGAL_postcondition(cache_.find(ct) != cache_.end());
      }
      if (cache_.find(ct) == cache_.end()) {
	++ui_->kinetic_certificates_;
	cache_.insert(typename Cache::value_type(ct,Cache_data(compute_exact_certificate(ct, et))));
	CGAL_postcondition(cache_.find(ct) != cache_.end());
      }
      CGAL_postcondition(cache_.find(ct) != cache_.end());
      while (cache_[ct].failure_time() < et) {
	cache_[ct].pop_failure_time();
	cache_[ct].pop_failure_time();
	++ui_->certificate_advances_;
	++ui_->certificate_advances_;
	CGAL_postcondition(cache_.find(ct) != cache_.end());
      }
      CGAL_UD_DEBUG("Got " << cache_[ct].failure_time() << std::endl);
      CGAL_UD_DEBUG("Check is " << check_.find(ct)->second.failure_time() << std::endl);
      CGAL_postcondition(cache_[ct].failure_time() == check_.find(ct)->second.failure_time());
      CGAL_postcondition(cache_.find(ct) != cache_.end());
      if (exact_failure_time(ct) > 1) {
	++ui_->unfailing_kinetic_certificates_;
      }
    }

    const Exact_time& exact_current_time() const {
      static Exact_time tmp;
      //ensure_exact_failure_time(current_time().data(), current_time());
      //Cert_tuple cct= ct.canonicalize();
      CGAL_UD_DEBUG("Approximate is " << current_time() << std::endl);
      if (current_time().data() == Cert_tuple()) {
	tmp= Exact_time(CGAL::to_interval(current_time()).first);
	return tmp;
      } else {
	if (!has_exact_failure_time(current_time().data())) {
	  ++ui_->exact_current_time_certificates_;
	}
	ensure_exact_failure_time(current_time().data(),
				  current_time());
	CGAL_UD_DEBUG("Exact current time is " 
		      << cache_[current_time().data()].failure_time() << std::endl);
	
	return cache_[current_time().data()].failure_time();
      }
    }

    const Exact_time& exact_failure_time(Cert_tuple ct) const {
      CGAL_precondition(cache_.find(ct) != cache_.end());
      return cache_.find(ct)->second.failure_time();
    }

  
    Exact_certificate compute_exact_certificate(Cert_tuple ct, Exact_time b) const {
      Exact_certificate  s= P::soc_(P::point(ct[0]), 
				    P::point(ct[1]),
				    P::point(ct[2]), 
				    P::point(ct[3]),
				    b, 1);
      return s;
    }






  
   

    INT join(INT a, INT b) const {
      return INT(std::min(a.inf(), b.inf()),
		 std::max(a.sup(), b.sup()));
    }


    Exact_time check_failure_time(Cert_tuple t) const {
      CGAL_precondition(check_.find(t) != check_.end());
      return check_.find(t)->second.failure_time();
    }

    void refine_from_exact(Event_key e, Cert_tuple t) const {
      P::simulator_handle()->event_time(e).refine(CGAL::to_interval(interval_from_exact_failure_time(t)));
    }

    
    CGAL::Comparison_result compare_concurrent(Event_key a,
					       Edge ea,
					       Event_key b,
					       Edge eb) const {
      CGAL_UD_DEBUG("Perturbing " << a << " and " << b << std::endl);
      
      Time ta= P::simulator_handle()->event_time(a);
      Time tb= P::simulator_handle()->event_time(b);
      Cert_tuple tua= ta.data();
      Cert_tuple tub= tb.data();
      CGAL_UD_DEBUG(tua << ": " << ta << " and " << tub << ": " << tb << std::endl);

      double asz= CGAL::to_interval(ta).second - CGAL::to_interval(ta).first;
      double bsz= CGAL::to_interval(tb).second - CGAL::to_interval(tb).first;
      Exact_time eta;
      Exact_time etb;
      if (asz > bsz) {
	if (!has_exact_failure_time(tua)) {
	  ++ui_->comparison_certificates_;
	}
	ensure_exact_failure_time(tua, CGAL::to_interval(ta).first);
	eta= exact_failure_time(tua);
	refine_from_exact(a, tua);

	if (eta > CGAL::to_interval(tb).second) {
	  CGAL_postcondition(check_failure_time(tua) > check_failure_time(tub));
	  return CGAL::LARGER;
	} else if (eta < CGAL::to_interval(tb).first) {
	  CGAL_postcondition(check_failure_time(tua) < check_failure_time(tub));
	  return CGAL::SMALLER;
	}
	if (!has_exact_failure_time(tub)) {
	  ++ui_->comparison_certificates_;
	}
	ensure_exact_failure_time(tub, CGAL::to_interval(tb).first);
	etb= exact_failure_time(tub);
	refine_from_exact(b, tub);
      } else {
	if (!has_exact_failure_time(tub)) {
	  ++ui_->comparison_certificates_;
	}
	ensure_exact_failure_time(tub, CGAL::to_interval(tb).first);
	etb= exact_failure_time(tub);
	refine_from_exact(b, tub);

	if (etb > CGAL::to_interval(ta).second) {
	  CGAL_postcondition(check_failure_time(tua) < check_failure_time(tub));
	  return CGAL::SMALLER;
	} else if (etb < CGAL::to_interval(ta).first) {
	  CGAL_postcondition(check_failure_time(tua) > check_failure_time(tub));
	  return CGAL::LARGER;
	}
	if (!has_exact_failure_time(tua)) {
	  ++ui_->comparison_certificates_;
	}
	ensure_exact_failure_time(tua, CGAL::to_interval(ta).first);
	eta= exact_failure_time(tua);
	refine_from_exact(a, tua);
      }
 
      CGAL_UD_DEBUG("Exact are " << eta << " and " << etb << std::endl);
      CGAL_postcondition(eta == check_.find(tua)->second.failure_time());
      CGAL_postcondition(etb == check_.find(tub)->second.failure_time());
      return CGAL::compare(eta, etb);
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


    bool can_fail(Cert_tuple ct) const {
      if (!ui_->is_active(ct[0]) && !ui_->is_active(ct[1])
	  && !ui_->is_active(ct[2]) && !ui_->is_active(ct[3])) {
	++ui_->constant_filtered_;
	return false;
      }
      CGAL::Protect_FPU_rounding<true> prot;
      
    
      INT curt= to_interval(P::simulator_handle()->current_time());
      INT rct(curt.inf(), 1);
      bool ret= (sign_at(ct[0],ct[1],ct[2],ct[3], rct) != CGAL::POSITIVE);
      
      if (!ret) {
	++ui_->interval2_filtered_;
	return false;
      } else {
	return true;
      }
    }
				
    bool isolate_failure(const Cert_tuple& ct, double lb, double ub, int depth,
			  Time &ret) const {
      int NS=20;
      double ld= lb;
      double step= (ub-ld)/NS;
      INT failures;
      bool has_zero=false;
      bool has_negative=false;
      bool certain=false;
      for (int i=0; i< NS; ++i) {
	double nld= std::min(ld+step, 1.0);
	{
	  INT ci(ld, nld);
	  ld= nld;
	  CGAL::Sign sn= sign_at(ct[0],ct[1],ct[2],ct[3], ci);
	  if (sn == CGAL::NEGATIVE) {
	    CGAL_assertion(has_zero); 
	    certain=true;
	    break;
	  } else if (sn == CGAL::ZERO) {
	    if (has_zero){
	      failures = join(failures, ci);
	    } else {
	      failures=ci;
	      has_zero=true;
	    }
	  } else if (sn == CGAL::POSITIVE) {
	    if (has_negative) break;
	  }
	}
	{
	  INT ci(nld, nld);
	  CGAL::Sign sn= sign_at(ct[0],ct[1],ct[2],ct[3], ci);
	  if (sn == CGAL::NEGATIVE) {
	    CGAL_assertion(has_zero); 
	    certain=true;
	    break;
	  } else if (sn == CGAL::POSITIVE) {
	    if (has_negative) break;
	  }
	}
	
      }
      if (has_negative || has_zero) {
	if (depth == 2 || .7*(ub-lb) < (failures.sup() - failures.inf())) {
	  ret= Time(failures.inf(), failures.sup(), ct);
	  return certain;
	} else {
	  return isolate_failure(ct, failures.inf(), failures.sup(), depth+1, ret);
	}
      } else {
	++ui_->interval3_filtered_;
	ret= Time();
	return certain;
      }
     
    }
    
    bool isolate_failure(Cert_tuple ct, double lb, Time &ret) const {
      CGAL::Protect_FPU_rounding<true> prot;
      return isolate_failure(ct, lb, 1.0, 0, ret);
    }


    typename Update_information::Handle ui_;
    mutable std::map<Cert_tuple,Cache_data> cache_;
    typename Simulation_traits::Kinetic_kernel::Positive_side_of_oriented_circle_2 soc_;
    std::map<Cert_tuple, Exact_certificate> check_;
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
      time_= CGAL::to_interval(tr.simulator_handle()->current_time()).second;
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
      CGAL::Protect_FPU_rounding<true> prot;
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
	  Kinetic_point_2 np(ENT(final(*it).x()),
			     ENT(final(*it).y()));
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

    typedef CGAL::Kinetic::Qt_widget_2<typename Simulation_traits::Simulator> Qt_gui;
    typedef CGAL::Kinetic::Qt_moving_points_2<Simulation_traits, Qt_gui> Qt_mps;
    typedef CGAL::Kinetic::Qt_triangulation_2<KDel, typename Simulation_traits::Instantaneous_kernel, Qt_gui> Qt_triangulation;
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
	CGAL_UD_DEBUG("Degeneracy with edge " 
		      << ks[0] << " " << ks[2] << std::endl);
      }
      return s!= CGAL::ON_NEGATIVE_SIDE;
    }
  }

 



  static int run(int argc, char *argv[], int n, int d, 
		 int seed, std::string ifile, std::string ffile) {
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
#undef CGAL_UD_DEBUG
