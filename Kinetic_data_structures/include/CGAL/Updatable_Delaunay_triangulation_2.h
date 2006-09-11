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

  - Want to make sure isolated intervals don't touch current_time
  even if they are otherwise certain. Ahhh, could just compare.

  - compute interval bound on root and see if that is good enough (use
  compare_concurrent if they overlap). Currently I make sure that
  the event actually occurs as this is easier than trying to handle
  events disappearing. There are some difficulties with new events
  that occur close to the current event. I think if I keep a last
  event around and use its inexact or exact rep to prune any new
  events I look at then I am fine. How to get the rep of the current
  time in the traits is a bit tricky.

  - So for a tuple/interval pair either I never computed the exact
  root, in which case I can take the first one or any +2 still
  before all other events. That is complicated. If I computed it, in
  which case I could store it.

  - for advancing, if I don't have an exact root, look for a disjoint
  interval (since I know none of the later roots overlap). Otherwise
  find the next root after the current exact root.

  - get rid of safe_exact_root

  - things are very delicate right now since you can't force the
  creation of exact roots for previous times (like current_time()
  after processing is done). 

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



  typedef CGAL::Interval_nt_advanced INT;
  

  typedef CGAL::Gmpq NT;
  typedef CGAL::POLYNOMIAL::Polynomial<NT> Function;
  typedef CGAL::POLYNOMIAL::Sturm_root_stack_traits<Function> Root_stack_traits;
  typedef CGAL::POLYNOMIAL::Sturm_root_stack<Root_stack_traits> Root_stack;
  typedef CGAL::POLYNOMIAL::Kernel<Function, Root_stack> Function_kernel;


  typedef typename CGAL::Kinetic::Handle_degeneracy_function_kernel<Function_kernel, true>  KK_function_kernel;
  typedef typename CGAL::Kinetic::Cartesian_kinetic_kernel<KK_function_kernel> Kinetic_kernel;

  typedef typename Kinetic::Active_objects_update_vector<typename Kinetic_kernel::Point_2, Point_key> Active_points_2_table;


  typedef IndirectKernel Indirect_kernel;


  typedef typename Indirect_kernel::Geometric_point_2 Static_point_2;


  typedef typename Indirect_kernel::Swapable_container Points;
  typedef typename Kinetic_kernel::Point_2 Kinetic_point_2;
  typedef typename Kinetic_kernel::Motion_function Kinetic_coordinate;

  typedef typename Indirect_kernel::Current_coordinates IK_current_coordinates;



  typedef typename Kinetic_kernel::Certificate Exact_certificate;
  typedef typename Function_kernel::Root Exact_time;





  /*typedef typename Kinetic_kernel::Positive_side_of_oriented_circle_2::result_type Exact_certificate;
 
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
#ifndef NDEBUG
    Exact_certificate check_;
#endif
  };
  */





  struct Update_information: public Kinetic::Ref_counted<Update_information> {

   

    
    //typedef std::map<Cert_tuple,Cache_data> Cache;
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
    Update_information(typename Active_points_2_table::Handle aot,
		       Kinetic_kernel kk,
		       Indirect_kernel ik):aot_(aot),
					   ik_(ik), 
					   soc_(kk.positive_side_of_oriented_circle_2_object()),
					   active_(ik_.number_of_point_2s(), false),
					   coef_cache_(ik_.number_of_point_2s()),
					   failure_time_(ik_.number_of_point_2s(), 1){
      reset(false);
      clear_stats();
    }

    bool is_initializing() const {
      return is_init_;
    }

    void set_is_initializing(bool ft) const {
      is_init_=ft;
    }

    /*NT activate_time() const {
      return activate_time_;
    }

    void set_activate_time(NT t) {
      activate_time_=t;
    }

    bool has_activate_time() const {
      return activate_time_ >=0;
      }*/

    std::vector<Point_key>& activate_list() {
      return activate_list_;
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
	CGAL_precondition(ct.inf()>=0 && ct.sup() <=1);

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

    void set_is_editing(typename Active_points_2_table::Editing_state es) {
      aot_->set_is_editing(es);
    }

    void set_final_kernel(Indirect_kernel &fk){
      fk_=fk;
    }

    void reset(bool clear=true) {
      //activate_time_=-1;
      next_activation_=1.0;
      start_time_=-1;
      is_init_=false;
      if (clear) active_.reset();
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
		  bool stat=false) const {
      if (stat) ++static_certificates_;
      CNT qpx = bx - ax;
      CNT qpy = by - ay;
      CNT rpx = cx - ax;
      CNT rpy = cy - ay;
      CNT tpx = dx - ax;
      CNT tpy = dy - ay;
      CNT det=CGAL::det2x2_by_formula(qpx*tpy - qpy*tpx, tpx*(dx - bx) 
				      + tpy*(dy - by),
				      qpx*rpy - qpy*rpx, rpx*(cx - bx) 
				      + rpy*(cy - by));
      return det;
    }

   
    void point_changed(Point_key k){
      /*std::vector<typename Cache::iterator> tk;
      for (typename Cache::iterator it = cache_.begin(); it != cache_.end(); ++it) {
	if (it->first[0] == k
	    || it->first[1] == k
	    || it->first[2] == k
	    || it->first[3] == k) {
	  tk.push_back(it);
	  break;
	}
      }
      for (unsigned int i=0; i< tk.size(); ++i) {
	CGAL_UD_DEBUG("Erasing " << tk[i]->first << std::endl);
	cache_.erase(tk[i]);
	}*/
    }


    /*template <class Time>
    bool has_exact_failure_time(Time ct) const {
      return time.has_exact_failure_time();
	//cache_.find(ct) != cache_.end();
	} */   
   

    /*void ensure_exact_failure_time(Cert_tuple ct,
				   double lb) const {
      compute_exact_failure_time(ct, Exact_time(lb));
      }*/


    /*HMMMMMMM
    void advance_exact_failure_time(Cert_tuple ct) const {
      CGAL_precondition(cache_.find(ct) != cache_.end());
      ++certificate_advances_;
      cache_[ct].pop_failure_time();
      }*/


    /*std::pair<double,double> interval_from_exact_failure_time(Cert_tuple ct) const {
      CGAL_precondition(cache_.find(ct) != cache_.end());
      if (cache_[ct].failure_time() >= Exact_time(1)) {
	return std::make_pair(1.0, -1.0);
      } else {
	std::pair<double,double> ip= CGAL::to_interval(cache_[ct].failure_time());
	return ip;
      }
      }*/

    Exact_certificate compute_exact_failure_time(Cert_tuple ct, 
						 Exact_time et) const {
      CGAL_UD_DEBUG("Computing exact time for " << ct 
		    << " from " << et  << std::endl);
      ++kinetic_certificates_;
      Exact_certificate ec= soc_(point(ct[0]),
				 point(ct[1]),
				 point(ct[2]), 
				 point(ct[3]),
				 et, 1);
      if (ec.failure_time() > 1) {
	++unfailing_kinetic_certificates_;
      }
      CGAL_UD_DEBUG("Got " << ec.failure_time() << std::endl);
      return ec;
    }

    void update_exact_failure_time(Cert_tuple ct, Exact_time et,
				    Cert_tuple ot, 
				    Exact_certificate& oc) const {
      //CGAL_precondition(check_.find(ct) != check_.end());
      CGAL_UD_DEBUG("Updating exact time for  " << ct 
		    << " starting at " << et << std::endl);

      if (ot != ct) {
	CGAL_precondition(ot.opposite() == ct);
	CGAL_UD_DEBUG("Flipping exact time for " << ct 
		      << " from " << oc.failure_time()
		      << std::endl);
	++certificate_advances_;
	oc.pop_failure_time();
      }
      
      while (oc.failure_time() < et) {
	oc.pop_failure_time();
	oc.pop_failure_time();
	CGAL_UD_DEBUG("Advancing exact time for " << ct 
		      << " from " << oc.failure_time() 
		      << std::endl);
	++certificate_advances_;
	++certificate_advances_;
	if (oc.failure_time() > 1) {
	  ++unfailing_kinetic_certificates_;
	}
      }
      CGAL_UD_DEBUG("Got " << oc.failure_time() << std::endl);
      
    }
  
   

    std::pair<double,double> join(std::pair<double,double> a, INT b) const {
      return std::make_pair(std::min(a.first, b.inf()),
			    std::max(a.second, b.sup()));
    }


    CGAL::Sign sign_at(Point_key a, Point_key b, 
		       Point_key c, Point_key d,
		       INT ct) const {
      INT det= incircle(cur_interval(a,ct,0),
			cur_interval(a,ct,1),
			cur_interval(b,ct,0),
			cur_interval(b,ct,1),
			cur_interval(c,ct,0),
			cur_interval(c,ct,1),
			cur_interval(d,ct,0),
			cur_interval(d,ct,1),
			true);
      
      if (det.sup() < 0) return CGAL::NEGATIVE;
      else if (det.inf() > 0) return CGAL::POSITIVE;
      else return CGAL::ZERO;
    }


    enum Isolate_result {NO_FAILURE=-1, POSSIBLE_FAILURE=0, CERTAIN_FAILURE=1};
				
    Isolate_result isolate_failure(const Cert_tuple& ct, double lb, double ub, 
				   bool starts_positive, int rem_depth,
				   int NS,
				   std::pair<double,double> &ret) const {
      //const int NS=10;
      double ld= lb;
      const double growth=1.5;
      double step= (ub-lb)/std::pow(growth,NS-1);
      std::pair<double,double> failures(std::numeric_limits<double>::infinity(),
					-std::numeric_limits<double>::infinity());
      bool has_zero=false;
      //bool has_negative=false;
      bool has_positive = starts_positive;
      Isolate_result certain=POSSIBLE_FAILURE;
      
      /*if (!has_positive) {
	CGAL::Sign sn= sign_at(ct[0],ct[1],ct[2],ct[3], INT(lb, lb));
	if (sn== CGAL::POSITIVE) has_positive=true;
	}*/


      CGAL_UD_DEBUG( "Interval is " << lb << " to " << ub 
		     << " with initial step " << step << "(" << rem_depth << ")" 
		     << std::endl);

      for (int i=0; i< NS; ++i) {
	double nld= std::min(ld+step, ub);
	if (nld == ld) break;

	step *= growth;

	CGAL_assertion(i != NS-1 || nld == ub);

	INT ci(ld, nld);
	ld= nld;
	CGAL::Sign csn= sign_at(ct[0],ct[1],ct[2],ct[3], ci);
	CGAL_UD_DEBUG("Sign on " << ci << " is " << csn << "(" << i << ")" 
		      << std::endl);
	if (csn == CGAL::NEGATIVE) {
	  //CGAL_assertion(has_zero); 
	  // could start negative in the (large) initial interval
	  //has_negative=true;
	  if (has_positive) {
	    certain=CERTAIN_FAILURE;
	    break;
	  }
	} else if (csn == CGAL::ZERO) {
	  has_zero=true;
	  //if (has_zero){
	  failures = join(failures, ci);
	  INT ci(nld, nld);
	  CGAL::Sign sn= sign_at(ct[0],ct[1],ct[2],ct[3], ci);
	  if (sn == CGAL::NEGATIVE) {
	    CGAL_assertion(csn != CGAL::POSITIVE);
	    if (has_positive) {
	      certain=CERTAIN_FAILURE;
	      break;
	    }
	  } else if (sn == CGAL::POSITIVE) {
	    CGAL_assertion(csn != CGAL::NEGATIVE);
	    has_positive=true;
	  }
	  
	} else if (csn == CGAL::POSITIVE) {
	  has_positive=true;
	}
	//if (i==0) {
      }







      if ( has_zero) {
	if ((certain==CERTAIN_FAILURE 
	     && (lb != failures.first || starts_positive)) 
	    || rem_depth == 0
	    || lb == failures.first && ub == failures.second
	    /* || .9*(ub-lb) < (failures.second - failures.first)*/) {
	  CGAL_UD_DEBUG(  "Not recursing with " 
			  << failures.first << " " << failures.second 
			  << " because " << has_zero << has_positive 
			  << certain << std::endl);
	  ret= failures;
	  return certain;
	} else {
	  CGAL_UD_DEBUG( "Recursing with " 
			 << failures.first << " " << failures.second 
			 << " because " << has_zero << has_positive 
			 << certain << std::endl);
	  double nwid = failures.second-failures.first;
	  double wid = ub-lb;
	  if (nwid > wid/2.0) NS*=2;
	  if (nwid < wid/4.0) NS/=2;
	  return isolate_failure(ct, failures.first, failures.second,
				 has_positive,
				 NS,
				 rem_depth-1, ret);
	}
      } else {
	++interval3_filtered_;
	//ret= std::pair<double,double>(1,-1);
	return NO_FAILURE;
      }
     
    }
    
    Isolate_result isolate_failure(Cert_tuple ct, double lb, double ub, bool starts_positive,std::pair<double,double> &ret) const {
      CGAL::Protect_FPU_rounding<true> prot;
      return isolate_failure(ct, lb, ub, starts_positive, 7, 10, ret);
    }

    
    double next_activation() const {
      return next_activation_;
    }
    void set_next_activation(double f) {
      next_activation_=f;
    }
  
    double start_time() const {
      return start_time_;
    }
    void set_start_time(double f) {
      start_time_=f;
    }
  


    const typename Active_points_2_table::Data &point(Point_key pk) const {
      return aot_->at(pk);
    }
 
    /*typename Default_traits::Simulator::Handle simulator_handle() {
      return dt_.simulator_handle();
      }

      typename Default_traits::Simulator::Const_handle simulator_handle() const {
      return dt_.simulator_handle();
      }

      typename Default_traits::Instantaneous_kernel
      instantaneous_kernel_object() const {
      return dt_.instantaneous_kernel_object();
      }*/


    // Default_traits dt_;
    typename Active_points_2_table::Handle aot_;
    //typename Simulation_traits::Simulator::Handle sim_;
    Indirect_kernel ik_, fk_;
    typename Kinetic_kernel::Positive_side_of_oriented_circle_2 soc_;
    Active active_;
    std::vector<Coef_data > coef_cache_;
    std::vector<double> failure_time_;
    //mutable std::map<Cert_tuple,Cache_data> cache_;
    //std::map<Cert_tuple, Exact_certificate> check_;
    bool is_init_;
    //NT activate_time_;
    std::vector<Point_key> activate_list_;
    double next_activation_;
    double start_time_;

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




  struct Refiner {
    enum State {INTERVAL, EXACT, INVALID};

    Refiner(Cert_tuple t, typename Update_information::Handle tbl): tuple_(t), ui_(tbl), state_(INTERVAL){}
    Refiner(): state_(INVALID) {
    }
    bool operator==(const Refiner &o) const {
      if (state_ != o.state_) return false;
      CGAL_assertion(state_==INVALID && o.state_==INVALID);
      return true;
    }
    
    typedef Exact_time Exact_root;
    /*This operator-() const {
      return *this;
      }*/
    
    bool refine(std::pair<double,double> &iv) const {
      CGAL_precondition(ui_ != typename Update_information::Handle());
      if (iv.first == iv.second) return false;
      if (has_exact_root()) return false;

      double dd= iv.second-iv.first;
      if (dd > .000001) {
	CGAL::Protect_FPU_rounding<true> prot;
	std::pair<double,double> oiv=iv;
	ui_->isolate_failure(tuple(), iv.first, iv.second, 1, 20, true, iv);
	return oiv != iv;
      } else {
	// compute exact;
	CGAL_assertion(!has_exact_root());
	++ui_->comparison_certificates_;
	
	ensure_exact_root(iv);
	iv= CGAL::to_interval(exact_root());
	return true;
      }
    }

    const Exact_root& 
    exact_root() const {
      /*if (ui_ == typename Update_information::Handle()) {
	CGAL_assertion(0);
	static Exact_root er;
	return er;
	} else {*/
      return cert_.failure_time();
      //}
    }

    const Exact_root& 
    exact_root(const std::pair<double,double> &iv) const {
      /*if (ui_ == typename Update_information::Handle()) {
	CGAL_assertion(0);
	static Exact_root er;
	return er;
	} else {*/
      ensure_exact_root(iv);
      return exact_root();
      //}
    }

    void write(std::ostream &out) const {
      if (tuple_ == Cert_tuple()) {
      } else{
	out << tuple_;
	if (has_exact_root()) {
	  out << " " << exact_root();
	}
      }
    }

    const Cert_tuple tuple() const {
      return tuple_;
    }

    bool equal_description(const Refiner &o) const {
      return tuple_== o.tuple_;
    }

    bool has_exact_root() const {
      return state_== EXACT;
    }
    
    void ensure_exact_root(const std::pair<double,double> &iv) const {
      if (!has_exact_root()) {
	state_=EXACT;
	cert_= ui_->compute_exact_failure_time(tuple_, iv.first);
	CGAL_assertion(check_.failure_time() == cert_.failure_time());
      }
    }
    void ensure_exact_root(const Exact_time &et) const {
      if (!has_exact_root()) {
	cert_= ui_->compute_exact_failure_time(tuple_, et);
	CGAL_assertion(check_.failure_time() == cert_.failure_time());
      }
    }

    const Exact_certificate& exact_certificate() const {
      return cert_;
    }

    std::pair<double,double> interval_from_exact() const {
      CGAL_precondition(has_exact_root());
      if (cert_.failure_time() >= Exact_time(1)) {
	return std::make_pair(1.0, -1.0);
      } else {
	std::pair<double,double> ip= CGAL::to_interval(cert_.failure_time());
	return ip;
      }
    }

    void set_exact_certificate(const Exact_certificate& ec) {
      CGAL_precondition(!has_exact_root());
      cert_= ec;
      CGAL_postcondition(cert_.failure_time() == check_.failure_time());
    }

    Cert_tuple tuple_;
    typename Update_information::Handle ui_;
    mutable Exact_certificate cert_;
    mutable State state_;
#ifndef NDEBUG
    mutable Exact_certificate check_;
#endif
  };








  struct Simulation_traits {
    typedef SimTraits_base P;
    
    struct Sillyness {
      typedef Kinetic_kernel KK;
      typedef Active_points_2_table APT;
    };
    typedef typename Sillyness::KK Kinetic_kernel;
    
    typedef typename Sillyness::APT Active_points_2_table;
    Active_points_2_table* active_points_2_table_handle() {
      return ap_.get();
    }
    const Active_points_2_table* active_points_2_table_handle() const {
      return ap_.get();
    }
    
    typedef CGAL::Simple_cartesian<NT> Static_kernel;

    typedef typename Kinetic::Cartesian_instantaneous_kernel<Active_points_2_table, Static_kernel> Instantaneous_kernel;
    Instantaneous_kernel instantaneous_kernel_object() const {
      return Instantaneous_kernel(ap_, Static_kernel());
    }
    typedef typename Kinetic::Interval_simulator_traits<Refiner> Simulator_traits;
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

    Kinetic_kernel kinetic_kernel_object() const {
      return Kinetic_kernel();
    }

  protected:
    typename Active_points_2_table::Handle ap_;
    typename Simulator::Handle sim_;
  };


  typedef typename Simulation_traits::Simulator::Event_key Event_key;





  typedef CGAL::Delaunay_triangulation_2<IndirectKernel,
					 CGAL::Triangulation_data_structure_2<
    CGAL::Kinetic::Delaunay_triangulation_vertex_base_2<IndirectKernel>,
    CGAL::Kinetic::Delaunay_triangulation_face_base_2<Simulation_traits> > > Triangulation;
  typedef CGAL::Kinetic::internal::Triangulation_data_structure_helper_2<typename Triangulation::Triangulation_data_structure> TDS_helper;


  typedef typename Triangulation::Edge Edge;

  typedef typename Simulation_traits::Kinetic_kernel::Motion_function::NT ENT;


  typedef  Kinetic::Delaunay_triangulation_default_traits_2<Simulation_traits, Triangulation> Default_traits;



































    
  struct Traits: public Default_traits {
    typedef Default_traits  P;
  
    typedef typename Simulation_traits::Simulator::Time Time;
    struct Certificate_data{};
    typedef std::pair<Time, Certificate_data> Certificate_pair;

    typedef typename Default_traits::Triangulation Triangulation;
    typedef typename Default_traits::Point_2 Point_2;
    typedef typename Default_traits::Simulator Simulator;


    bool is_exact() const {
      return true;
    }
    /*typename Simulation_traits::Active_points_2_table::Handle
      active_points_2_table_handle() {
      return ui_->active_points_2_table_handle();
      }
    
   
      typename Default_traits::Simulator::Handle simulator_handle() {
      return ->simulator_handle();
      }

      typename Default_traits::Simulator::Const_handle simulator_handle() const {
      return ui_->simulator_handle();
      }*/

    Cert_tuple tuple(typename Default_traits::Edge e) const {
      Point_key ks[4];
      ks[0]= TDS_helper::origin(e)->point();
      ks[1]= TDS_helper::third_vertex(e)->point();
      ks[2]= TDS_helper::destination(e)->point();
      ks[3]= TDS_helper::mirror_vertex(e)->point();
      if (ks[1] == Point_key() || ks[3]==Point_key()
	  || ks[0] == Point_key() || ks[2]== Point_key()) return Cert_tuple();
      else return Cert_tuple(ks);
    }
    
    void point_changed(Point_key k){
      ui_->point_changed(k);
    }

    const Time &current_time() const {
      return P::simulator_handle()->current_time();
    }


    /*typename Default_traits::Instantaneous_kernel
      instantaneous_kernel_object() const {
      return ui_->instantaneous_kernel_object();
      }

      const Point_2 &point(Point_key pk) const {
      return ui_->point(pk);
      }*/


    Certificate_pair null_pair() const {
      return std::make_pair(Time(2), Certificate_data());
    }

    Certificate_pair return_pair(Time rt) const {
      return Certificate_pair(rt, Certificate_data());
    }
    

    double rational_current_time() const {
      return CGAL::to_interval(current_time()).first;
    }


    Traits(Simulation_traits st,
	   typename Update_information::Handle ui): P(st),
						    ui_(ui){}
    //soc_(tr.kinetic_kernel_object().positive_side_of_oriented_circle_2_object()){}
    
 

    bool can_fail(Cert_tuple ct, double end_time=1) const {
      if (!ui_->is_active(ct[0]) && !ui_->is_active(ct[1])
	  && !ui_->is_active(ct[2]) && !ui_->is_active(ct[3])) {
	++ui_->constant_filtered_;
	return false;
      }
      CGAL::Protect_FPU_rounding<true> prot;
      
    
      INT curt= to_interval(current_time());
      INT rct(curt.inf(), end_time);
      bool ret= (ui_->sign_at(ct[0],ct[1],ct[2],ct[3], rct) != CGAL::POSITIVE);
      
      if (!ret) {
	++ui_->interval2_filtered_;
	return false;
      } else {
	return true;
      }
    }


    Certificate_pair certificate_failure_time(typename Default_traits::Edge e) {

      Cert_tuple ct= tuple(e);
      if (ct == Cert_tuple()) return null_pair();

      double end_time= ui_->next_activation();
      //ui_->set_next_activation(1.0);
#ifndef NDEBUG
      Exact_time check_failure_time;
      Exact_certificate check_cert;
      {
	Exact_time ect;
	if (P::simulator_handle()->current_time().refiner() == Refiner()) {
	  ect= Exact_time(CGAL::to_interval(P::simulator_handle()->current_time()).first);
	} else {
	  ect= P::simulator_handle()->current_time().refiner().check_.failure_time();
	}
	
	check_cert=ui_->compute_exact_failure_time(ct, ect);
	check_failure_time= check_cert.failure_time();

	if (check_failure_time > end_time) check_failure_time= 2;
      }
#endif
   
      //std::pair<Cert_tuple, bool> ctp= Cert_tuple::make(ks);
      
      if (!can_fail(ct, end_time)) {
	CGAL_postcondition(check_failure_time > 1);
	return null_pair();
      }
      
      CGAL_UD_DEBUG("Computing failure time for " << ct << std::endl);
      CGAL_UD_DEBUG("Exact failure time is " << check_failure_time<< std::endl);


      double bt = ui_->start_time();
      if (bt <0) {
	bt= CGAL::to_interval(current_time()).first;
      }
      


      std::pair<double,double> ft;
      typename Update_information::Isolate_result isc= ui_->isolate_failure(ct, 
									    bt,
									    end_time,
									    true,
									    ft);


      if (isc == Update_information::NO_FAILURE) { 
	CGAL_UD_DEBUG("Not isolated" << std::endl);
	CGAL_assertion(check_failure_time >=1);
	return null_pair();
      }

      Time rett(ft.first, ft.second, Refiner(ct, ui_));
#ifndef NDEBUG
      rett.refiner().check_= check_cert;
#endif

      if (isc == Update_information::CERTAIN_FAILURE) {
	if (ft.first <= CGAL::to_interval(current_time()).second) {
	  // interval overlaps current
	  P::simulator_handle()->current_time().refiner().ensure_exact_root(CGAL::to_interval(P::simulator_handle()->current_time()));
	  const Exact_time &ect= P::simulator_handle()->current_time().refiner().exact_root();
	  
	  if (Exact_time(ft.first) <= ect) {
	    rett.refiner().ensure_exact_root(ect);
	    CGAL_UD_DEBUG("Exact" << rett.refiner().exact_root() <<  std::endl);
	    rett.set_interval(CGAL::to_interval(rett.refiner().exact_root()));
	  }
	}
	CGAL_UD_DEBUG("Returning isolated " << rett << std::endl);
	
	return return_pair(rett);
      } else {
	//if (!has_exact_failure_time(ct)) {
	++ui_->uncertain_exact_computations_;
	//	}
	CGAL_UD_DEBUG("Not isolated " <<  std::endl);
	P::simulator_handle()->current_time().refiner().ensure_exact_root(CGAL::to_interval(P::simulator_handle()->current_time()));
	const Exact_time &ect= P::simulator_handle()->current_time().refiner().exact_root();
	rett.refiner().ensure_exact_root(ect);
	CGAL_UD_DEBUG("Exact" << rett.refiner().exact_root() <<  std::endl);
	if (rett.refiner().exact_root() >= 1) {
	  CGAL_assertion(check_failure_time >=1);
	  CGAL_UD_DEBUG("Phantom root " << std::endl);
	  CGAL_postcondition(check_failure_time > 1);
	  return null_pair();
	} else {
	  rett.set_interval(CGAL::to_interval(rett.refiner().exact_root()));
	  CGAL_UD_DEBUG("Returning exact " << rett << std::endl);
	  return return_pair(rett);
	}
      }
    }




    Certificate_pair certificate_failure_time(typename Default_traits::Edge e, 
					      Certificate_data ) {
      const Time &curt= P::simulator_handle()->current_time();
      Cert_tuple ct= curt.refiner().tuple().opposite();
#ifndef NDEBUG
      Cert_tuple check_ct= tuple(e);
#endif
      CGAL_assertion(ct == check_ct);

      double end_time= ui_->next_activation();
      ui_->set_next_activation(1.0);
 
#ifndef NDEBUG
      Exact_time check_failure_time;
      Exact_certificate check_cert;
      {
	check_cert= curt.refiner().check_;
	check_cert.pop_failure_time();

	check_failure_time= check_cert.failure_time();
	if (check_failure_time > end_time) check_failure_time= 2;
	CGAL_UD_DEBUG("Exact time is " << check_cert.failure_time() << std::endl);

      }
#endif
      //std::pair<Cert_tuple, bool> ctp= Cert_tuple::make(ks);
      
      if (!can_fail(ct, end_time)) {
	CGAL_postcondition(check_failure_time > 1);
	return null_pair();
      }

      CGAL_UD_DEBUG("Advancing certificate for " << ct
		    << std::endl);

      // this depends on being the last event of the batch whose time is computed
      const Time& net= P::simulator_handle()->next_event_time();

      double lb= CGAL::to_interval(net).first;
      bool pos_start=false;
      if (lb == ui_->start_time()) {
	pos_start=true;
      } else {
	// done twiceish
	CGAL::Protect_FPU_rounding<true> prot;
	CGAL::Sign sn= ui_->sign_at(ct[0],ct[1],ct[2],ct[3], INT(lb, lb));
	if (sn != CGAL::POSITIVE) {
	  lb= CGAL::to_interval(current_time()).first;
	} else {
	  pos_start=true;
	}
      }

      std::pair<double,double> ft;
      typename Update_information::Isolate_result isc
	= ui_->isolate_failure(ct, lb, end_time, pos_start, ft);
      
  

      if (isc == Update_information::NO_FAILURE) { 
	CGAL_UD_DEBUG("No root there." << std::endl);
	CGAL_postcondition(check_failure_time > 1);
	return null_pair();
      }

      Time rett(ft.first, ft.second, Refiner(ct, ui_));
#ifndef NDEBUG
      rett.refiner().check_= check_cert;
#endif
      if (isc == Update_information::CERTAIN_FAILURE
	  && ft.first <= CGAL::to_interval(curt).second
	  || isc == Update_information::POSSIBLE_FAILURE
	  || curt.refiner().has_exact_root()) {
	if (isc == Update_information::POSSIBLE_FAILURE) {
	  ++ui_->uncertain_exact_computations_;
	}
	curt.refiner().ensure_exact_root(CGAL::to_interval(curt));
	Exact_certificate ec= curt.refiner().exact_certificate();
	ec.pop_failure_time();
	if (ec.failure_time() > 1) {
	  ++ui_->unfailing_kinetic_certificates_;
	  CGAL_UD_DEBUG("Phantom root " << std::endl);
	  CGAL_postcondition(check_failure_time > 1);
	  return null_pair();
	}
	++ui_->certificate_advances_;
	ft= CGAL::to_interval(curt.refiner().exact_root());
	rett.set_interval(ft);
	rett.refiner().set_exact_certificate(ec);
	return return_pair(rett);
      } else {
	CGAL_UD_DEBUG("Returning " << rett << std::endl);
	return return_pair(rett);
      }
    }

  

    CGAL::Comparison_result compare_concurrent(Event_key a,
					       Edge,
					       Event_key b,
					       Edge) const {
      return CGAL::compare(a,b);
    }
    

    typename Update_information::Handle ui_;
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

    bool test_and_add_one(Edge e, std::vector<Point_key> &active) const {
      ++ui_->static_certificates_;
      if (!compute_ok(e, ui_->fk_)) {
	//add(e.first->vertex(0)->point(), active);
	//add(e.first->vertex(1)->point(), active);
	//add(e.first->vertex(2)->point(), active);
	add(TDS_helper::mirror_vertex(e)->point(), active);
	return true;
      } else return false;
    }

    template <class TDS> 
    void initialize_events(const TDS &triangulation,
			   Indirect_kernel fk) {
      ui_->num_edges_=triangulation.tds().number_of_edges();
      ui_->set_final_kernel(fk);
      std::vector<Point_key> active;
      for (typename TDS::Finite_edges_iterator it= triangulation.finite_edges_begin(); 
	   it != triangulation.finite_edges_end(); ++it){
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
      //test_and_add(e, active);
      //CGAL_precondition(ui_->is_active(e.first->vertex((e.second+1)%3)->point()));
      //CGAL_precondition(ui_->is_active(e.first->vertex((e.second+2)%3)->point()));
      //CGAL_precondition(ui_->is_active(e.first->vertex(e.second)->point()));
      add(e.first->vertex(0)->point(), active);
      add(e.first->vertex(1)->point(), active);
      add(e.first->vertex(2)->point(), active);
      Edge em= TDS_helper::mirror_edge(e);
      //CGAL_precondition(ui_->is_active(em.first->vertex(em.second)->point()));
      add(em.first->vertex(em.second)->point(), active);
      //test_and_add_one(Edge(e.first, (e.second+1)%3), active);
      //test_and_add_one(Edge(e.first, (e.second+2)%3), active);
      //test_and_add_one(Edge(em.first, (em.second+1)%3), active);
      //test_and_add_one(Edge(em.first, (em.second+2)%3), active);
      if (!active.empty()) {
	Interpolate_event ev(tr_, 
			     ui_,
			     ui_->ik_.current_coordinates_object(),
			     ui_->fk_.current_coordinates_object(), 
			     active.begin(), active.end());
	if (!ev.empty()){
	  tr_.simulator_handle()->new_event(ev.time(), ev);
	  //active_.set(active.begin(), active.end());
	  ++ui_->num_interpolations_;
	  /*for (unsigned int i=0; i< active.size(); ++i){
	    ui_->activate(active[i], true);
	    }*/
	  ui_->set_next_activation(CGAL::to_interval(ev.time()).second);
	}
      }

      bool ok=true;
      double skip_to= CGAL::to_interval(tr_.simulator_handle()->next_event_time()).first;
      {
	CGAL::Protect_FPU_rounding<true> prot;
	CGAL::Sign sn= ui_->sign_at(e.first->vertex(0)->point(), e.first->vertex(1)->point(),
				    e.first->vertex(2)->point(), em.first->vertex(em.second)->point(),
				    INT(skip_to));
	if (sn == CGAL::POSITIVE) {
	  for (int i=0; i< 3; ++i) {
	    if (i != e.second) {
	      CGAL::Sign sn= ui_->sign_at(e.first->vertex(0)->point(), e.first->vertex(1)->point(),
					  e.first->vertex(2)->point(), e.first->mirror_vertex(i)->point(),
					  INT(skip_to));
	      if (sn != CGAL::POSITIVE) {
		ok=false;
		break;
	      }
	    }
	    if (i != em.second) {
	      CGAL::Sign sn= ui_->sign_at(em.first->vertex(0)->point(), em.first->vertex(1)->point(),
					  em.first->vertex(2)->point(), em.first->mirror_vertex(i)->point(),
					  INT(skip_to));
	      if (sn != CGAL::POSITIVE) {
		ok=false;
		break;
	      }
	    }
	  }
	}
      }
      if (ok) {
	CGAL_UD_DEBUG("Skipping to " << skip_to << " from " << 	tr_.simulator_handle()->current_time() << std::endl);
	ui_->set_start_time(skip_to);
      } else {
	ui_->set_start_time(-1);
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


    void reset() {
      ui_->reset();
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
      c[1]=(NT(ip[i])-NT(fp[i]))/NT(time-1);
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
      // make sure it is clear of the event to not prompt building of exact
      /*if (tr.simulator_handle()->current_time()
	== typename Simulation_traits::Simulator::NT(0)) {
	time_=0;
	} else {*/
      time_= CGAL::to_interval(tr.simulator_handle()->current_time()).second;
      //}
      std::sort(b,e);
      It ne= std::unique(b,e);
      for (It c=b; c!= ne; ++c){
	if (ic(*c) != fc(*c)) {
	  motions_.push_back(MP(*c, interpolate_t1(time_, ic(*c), fc(*c))));
	} else {
	  CGAL_UD_DEBUG("Point " << *c << " doesn't move." << std::endl);
	}
      }
      if (!empty()) {
	CGAL_UD_DEBUG("Will Interpolate at " << time());
      }
    }
    
    typename Simulation_traits::Simulator::Time time() const {
      return time_;
    }

    std::ostream & write(std::ostream&out) const {
      out << "Update ";
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
      ui_->set_is_editing(Active_points_2_table::LOGGED);
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
      ui_->set_is_editing(Active_points_2_table::NOT);
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
      tr_.active_points_2_table_handle()->set_is_editing(Active_points_2_table::LOGGED);
      
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

      tr_.active_points_2_table_handle()->set_is_editing(Active_points_2_table::NOT);
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
    typename Update_information::Handle ui
      = new Update_information(tr_.active_points_2_table_handle(),
			       tr_.kinetic_kernel_object(),
			       ik_);
    
    Traits traits(tr_,ui);

    traits.active_points_2_table_handle()->set_is_editing(Active_points_2_table::UNLOGGED);
    for (It c=b; c!= e; ++c){
      traits.active_points_2_table_handle()->insert(Kinetic_point_2(Kinetic_coordinate(c->x()),
								    Kinetic_coordinate(c->y())));
    }
    traits.active_points_2_table_handle()->set_is_editing(Active_points_2_table::NOT);

  
    
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
    kdel_->visitor().reset();
    kdel_->set_has_certificates(false);
    tr_.simulator_handle()->set_interval(0,1);
    tr_.active_points_2_table_handle()->set_is_editing(Active_points_2_table::UNLOGGED);
    typename Indirect_kernel::Current_coordinates cc= ik_.current_coordinates_object();
    for (typename Simulation_traits::Active_points_2_table::Key_iterator
	   kit= tr_.active_points_2_table_handle()->keys_begin();
	 kit !=  tr_.active_points_2_table_handle()->keys_end(); ++kit) {
      tr_.active_points_2_table_handle()->set(*kit, Kinetic_point_2(Kinetic_coordinate(cc(*kit).x()),
								    Kinetic_coordinate(cc(*kit).y())));
    }
    tr_.active_points_2_table_handle()->set_is_editing(Active_points_2_table::NOT);


    Indirect_kernel fk;
    fk.new_point_2s(pts.begin(), pts.end());
    kdel_->set_has_certificates(true, true);
    kdel_->visitor().initialize_events(kdel_->triangulation(), fk);
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
  
  /*static bool is_hull_edge(const Edge &e) {
    return ! TDS_helper::mirror_vertex(e)->point().is_valid()
      || ! TDS_helper::third_vertex(e)->point().is_valid()
      || ! TDS_helper::origin(e)->point().is_valid()
      || ! TDS_helper::destination(e)->point().is_valid();
      }*/

  
  static bool compute_ok(const Edge &e,  Indirect_kernel sk) {
    //typename Indirect_kernel::Current_coordinates 
    //cc= sk.current_coordinates_object();
    Point_key ks[4];
    ks[0]= TDS_helper::origin(e)->point();
    ks[1]= TDS_helper::third_vertex(e)->point();
    ks[2]= TDS_helper::destination(e)->point();
    ks[3]= TDS_helper::mirror_vertex(e)->point();
    for (unsigned int i=0; i< 4; ++i){
      if (ks[i]==Point_key()) {
	return true;
      }
    }
    
    typename Indirect_kernel::Side_of_oriented_circle_2 soc
      = sk.side_of_oriented_circle_2_object();
      
    CGAL::Oriented_side s=soc(ks[0], ks[1], ks[2], ks[3]);
    
    if (s== CGAL::ON_ORIENTED_BOUNDARY) {
      CGAL_UD_DEBUG("Degeneracy with edge " 
		    << ks[0] << " " << ks[2] << std::endl);
    }
    return s!= CGAL::ON_NEGATIVE_SIDE;
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
