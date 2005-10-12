#ifndef CGAL_POLYNOMIAL_ROOT_BOUND_SOLVER_CORE_H
#define CGAL_POLYNOMIAL_ROOT_BOUND_SOLVER_CORE_H

#include <CGAL/Polynomial/basic.h>
//#include <CGAL/Polynomial/Tools/Isolating_interval.h>
//#include <CGAL/Polynomial/Tools/Simple_interval_root.h>
#include <CGAL/Polynomial/internal/Descartes_root_count.h>


#define CGAL_DSPRINT(x)

CGAL_POLYNOMIAL_BEGIN_NAMESPACE
//! A Solver which uses
template <class Traits_t>
class Upper_bound_root_stack {
protected:
  typedef Upper_bound_root_stack<Traits_t> This;

  typedef internal::Descartes_root_count Root_count;
 
  struct Interval_info;
  typedef typename Traits_t::Isolating_interval Interval;

  typedef std::vector<Interval_info> Intervals;
  typedef typename Traits_t::Sturm_root_count Sturm_root_counter;
  typedef typename Traits_t::Function Polynomial;
  typedef typename Traits_t::NT NT;
public:
  typedef Traits_t Traits;
  typedef typename Traits_t::Root Root;
  /*struct Root: public Traits_t::Root {
    Root(){}
    template <class TNT>
    Root(TNT nt): Traits_t::Root(nt){}
    Root(const NT &nt, bool is_odd, Traits_t k): Traits_t::Root(nt, is_odd, k){}
    Root(const Interval &nt, bool is_odd=true, Traits_t k= Traits_t()):
      Traits_t::Root(nt, is_odd, k){}
    Root(const Interval &ii, const typename Traits_t::Function &sa, 
	 Sign slb, Sign sub, 
	 Traits_t k): Traits_t::Root(ii, sa, slb, sub, k){}
  };
  friend class Root;*/

  Upper_bound_root_stack(): cur_(Root::infinity()) {
  };
  Upper_bound_root_stack(const Polynomial &f, 
			 const Root &lb,
			 const Root &ub,
			 const Traits &tr): kernel_(tr), f_(f),
					    lb_(lb), ub_(ub),
					    rc_(kernel_.root_count_object(f_)),
					    has_ss_(false){
    initialize();

    do {
      initialize_intervals();
      cur_= make_next_root();
    } while (!(cur_ >lb_));
  };

  const Root& top() const {
    return cur_;
  }

  void pop() {
    initialize_intervals();
    cur_= make_next_root();
  }

  bool empty() const {
    return cur_==Root::infinity();
  }

protected:


  Root make_next_root(){
    if (intervals_.empty()){
      //std::cout << "Out of intervals.\n";
      return Root::infinity();
    }

    CGAL_precondition(!intervals_.back().num_roots.is_zero());
    Root ret;

    if (intervals_.back().interval.is_singular()) ret= Root(intervals_.back().interval,
							     intervals_.back().num_roots.is_odd());
    else ret= Root(intervals_.back().interval,f_, intervals_.back().left_sign,
		   intervals_.back().right_sign,  kernel_);
    CGAL_DSPRINT(std::cout << "Trying root " << ret << std::endl);
    CGAL_DSPRINT(std::cout << "Root count was " << intervals_.back().num_roots << std::endl);
    if ( !(ub_> ret)){
      CGAL_DSPRINT(std::cout << "Rejection due to >=" << ub_ << std::endl);
      intervals_.clear();
      return Root::infinity();
    } else {
      intervals_.pop_back();
      return ret;
    }
  }

  void initialize(){
    CGAL_Polynomial_expensive_precondition(lb_<ub_);
    
    if (f_.is_constant()) return;
    /*else if (f_.is_linear()){
      if (handle_linear()) return;
      }*/
    
    //std::cout << "Initializing on " << f_ << " with bounds " << lb_ << " and " << ub_ << std::endl;
    
    typename Traits::Root_bound rbe= kernel_.root_bound_object();
    NT rb= rbe(f_);
    Interval lbi, ubi;
    if (lb_ == -Root::infinity()) {
      lbi= Interval(-rb);
    } else {
      lbi= lb_.interval(); //power_of_two(lb_.interval().lower_bound());
    }
    if (ub_== Root::infinity()){
      ubi= Interval(rb);
    } else {
      ubi= ub_.interval();
    }
    Interval ii= lbi || ubi;
    // make up a sign
    //Root_count rc= compute_root_count(ii, POLYNOMIAL_NS::POSITIVE, POLYNOMIAL_NS::POSITIVE); //ii.apply_to_interval(rc_);

    intervals_.push_back(Interval_info(ii, // skip rc
				       ii.apply_to_endpoint(kernel_.sign_at_object(f_),
							    Interval::LOWER),
				       ii.apply_to_endpoint(kernel_.sign_at_object(f_), Interval::UPPER)));
    CGAL_DSPRINT(std::cout << "Initial interval is " << ii <<std::endl);
  }

  void subdivide_front(){
   
    
    //assert(!intervals_.back().is_good());
    
    Interval_info ii= intervals_.back();
    CGAL_DSPRINT(std::cout << "Investigating " << ii.interval << std::endl);
    intervals_.pop_back();      
    
    Interval uh= ii.interval.second_half();
    Interval lh= ii.interval.first_half();
    Interval mi= lh.upper_endpoint_interval();
    CGAL_POLYNOMIAL_NS::Sign sm= uh.apply_to_endpoint(kernel_.sign_at_object(f_), Interval::LOWER);
    //CGAL::Sign lms= sm;
    //CGAL::Sign ums= sm;
    //std::cout << "Splitting at " << mid << std::endl;
    bool mid_is_root=(sm == CGAL_POLYNOMIAL_NS::ZERO);
    
    /*if (mid_is_root){
      Interval mi1= lh.second_half();
      Interval mi2= uh.first_half();
      Interval mhc= mi1|| mi2;
      Root_count uhc = compute_root_count(mhc, CGAL::POSITIVE, CGAL::POSITIVE); // just make up signs
      CGAL_assertion(uhc != Root_count::zero());
      if (uhc== Root_count::one()){
      uh= uh.second_half();
      lh= lh.first_half();
      lms= lh.apply_to_endpoint(kernel_.sign_at_object(f_), Interval::UPPER);
      ums= uh.apply_to_endpoint(kernel_.sign_at_object(f_), Interval::LOWER);
      std::cout << "Skipping interval " << mhc << " around exact root " << mi << "\n";
      }
      }*/
    
    //Root_count uhc= compute_root_count(uh, sm, ii.right_sign); //uh.apply_to_interval(rc_);
    //if (!uhc.is_zero()){
    CGAL_DSPRINT(std::cout << "Produced u interval of " << uh << std::endl);
    intervals_.push_back(Interval_info(uh, sm, ii.right_sign));
    //}
    if (mid_is_root) {
      // need to check if it is an even root
	/*Root_count rc;
	  if (deg%2==0) rc= Root_count::even(); //single_even
	  else rc= Root_count::odd(); //single_odd*/
      intervals_.push_back(Interval_info(mi, CGAL_POLYNOMIAL_NS::ZERO, CGAL_POLYNOMIAL_NS::ZERO));
    }
    
    
    //Root_count lhc= compute_root_count(lh, ii.left_sign, sm); //lh.apply_to_interval(rc_);
    //if (!lhc.is_zero()) {
    //std::cout << "Produced l interval of " << lh << std::endl;
    intervals_.push_back(Interval_info(lh, ii.left_sign, sm));
    //}
    CGAL_DSPRINT(write_intervals(std::cout));
  }


  bool sturm_front() {
    if (!has_ss_){
      ss_= kernel_.Sturm_root_count_object(f_);
      has_ss_=true;
    }
    CGAL_DSPRINT(std::cout << "Computing sturm for " << intervals_.back().interval << std::endl);
    //std::cout << "Computing sturm for " << intervals_.back().interval << std::endl;
    unsigned int ct= intervals_.back().interval.apply_to_interval(ss_);
    if (ct == 0){
      intervals_.pop_back();
      return true;
    } else if (ct==1) {
	//if (intervals_.back().left_sign == intervals_.back().right_sign){
	  //  intervals_.back().num_roots =1; //single_even
	  //} else {
	  ; //single_odd
	  //}
	  intervals_.back().num_roots =1;
    } else {
      intervals_.back().num_roots=Root_count(ct);
    }
    return false;
  }


  void upperbound_front() {
    intervals_.back().num_roots=compute_root_count(intervals_.back().interval, 
						   intervals_.back().left_sign,
						   intervals_.back().right_sign); //ii.apply_to_interval(rc_);
  }

  void count_singular_front(){
    assert(intervals_.back().num_roots.is_unknown());
    int deg=intervals_.back().interval.apply_to_endpoint(kernel_.multiplicity_object(f_), Interval::LOWER);
    intervals_.back().num_roots=deg;
  }

  //! establish the invariant 
  /*!
    The front interval contains only one root. 
  */
  void initialize_intervals(){
    //write_intervals(std::cout);
    while (!intervals_.empty()){
       CGAL_assertion(intervals_.back().left_sign 
	     == intervals_.back().interval.apply_to_endpoint(kernel_.sign_at_object(f_), Interval::LOWER));
      CGAL_assertion(intervals_.back().right_sign 
	     == intervals_.back().interval.apply_to_endpoint(kernel_.sign_at_object(f_), Interval::UPPER));


      // compute a root count
      if (intervals_.back().interval.is_singular()) {
	count_singular_front();
	//std::cout << "Breaking in singular." << std::endl;
	break;
      } else if (intervals_.back().interval.approximate_width() <= min_interval_width()){ 
	if (sturm_front()) continue;
	if (intervals_.back().num_roots.is_single()
	    && intervals_.back().left_sign != ZERO 
	    && intervals_.back().right_sign != ZERO) break;
      } else  if (intervals_.back().num_roots.is_unknown()) {
	upperbound_front();

	if (intervals_.back().num_roots.is_zero()){
	  intervals_.pop_back();
	  continue;
	}
	
	if (intervals_.back().num_roots.is_single()
	    && intervals_.back().left_sign != ZERO
	    && intervals_.back().right_sign != ZERO) {
	  CGAL_assertion( intervals_.back().left_sign != intervals_.back().right_sign);
	  //std::cout << "Breaking in normal." << std::endl;
	  break;
	}
      }

      subdivide_front();
      
    } // while (the first has more than one interval);
    //std::cout << "Initialized.\n";
  }

  Root_count compute_root_count(const Interval &ii,
				CGAL_POLYNOMIAL_NS::Sign sl,
				CGAL_POLYNOMIAL_NS::Sign sr){
    return Root_count(ii.apply_to_interval(rc_, sl, sr));
  }
  
  void write_intervals(std::ostream &out){
    for (unsigned int i=0; i< intervals_.size(); ++i){
      out << "(" << intervals_[i].interval << ", " << intervals_[i].left_sign << intervals_[i].right_sign << ")";
    }
    std::cout << std::endl;
  }
  static double min_interval_width() {
    return .000001;
  }

  //! Just return small intervals if doubles are used
  template <class NTT>
  static bool is_small(const NTT &, double) {
    return false;
  }

  static bool is_small(double, double wid)  {
    return wid < .00001;
  }

  struct Interval_info{
    Root_count num_roots;
    CGAL_POLYNOMIAL_NS::Sign left_sign;
    CGAL_POLYNOMIAL_NS::Sign right_sign;
    Interval_info(Interval in, /*Root_count num,*/ CGAL_POLYNOMIAL_NS::Sign ls, 
		  CGAL_POLYNOMIAL_NS::Sign rs):
      num_roots(Root_count::UNKNOWN), left_sign(ls), right_sign(rs), interval(in){}
    Interval interval;
    /*bool is_good(){
      if (is_small(typename Traits::Function::NT(0), interval.approximate_width())) {
	assert(0);
	if (left_sign == right_sign) num_roots=Root_count::even(); //single_even
	else num_roots=Root_count::odd(); //single_odd
	interval=interval.upper_endpoint_interval();
	return true;
      }
      if (interval.is_singular()) return true;
      else if (left_sign != POLYNOMIAL_NS::ZERO && right_sign != POLYNOMIAL_NS::ZERO){
 	return num_roots.is_single();
      } else return false;
      }*/
  };  

  Traits kernel_;
  Polynomial f_;
  Intervals intervals_;
  Root lb_, ub_;
  Root cur_;
  typename Traits::Root_count rc_;
  Sturm_root_counter ss_;
  bool has_ss_;
};
CGAL_POLYNOMIAL_END_NAMESPACE

/*namespace std {
  template <class Traits>
  class numeric_limits<typename CGAL::POLYNOMIAL::Upper_bound_root_stack<Traits>::Root>: 
    public numeric_limits<typename Traits::Root> {};
    };*/

#endif
