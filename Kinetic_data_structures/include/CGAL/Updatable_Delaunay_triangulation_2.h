#ifndef UPDATABLE_DELAUNAY_TRIANGULATION_2_H
#define UPDATABLE_DELAUNAY_TRIANGULATION_2_H

#ifdef NDEBUG
#define CGAL_UD_DEBUG(x)
#else
#define CGAL_UD_DEBUG(x) std::cout << x
#endif

#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Interval_simulator_traits.h>
#include <CGAL/Kinetic/Simulation_traits.h>
#include <CGAL/Kinetic/Free_event_base.h>
#include <CGAL/Kinetic/Delaunay_triangulation_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_vertex_base_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_recent_edges_visitor_2.h>
#include <CGAL/Kinetic/Delaunay_triangulation_visitor_base_2.h>
#include <CGAL/Updatable_Delaunay_triangulation_table_2.h>
#include <CGAL/Indirect_point_2_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


#include <CGAL/Kinetic/IO/Qt_moving_points_2.h>
#include <CGAL/Kinetic/IO/Qt_triangulation_2.h>
#include <CGAL/Kinetic/IO/Qt_widget_2.h>



#include <vector>


CGAL_BEGIN_NAMESPACE

#define NEWTON



/*
  Open questions:
  - storage of certs
  - priority queue for deletions

  Where to go next:

  - optimize priority queue--insert is 20%

  - Delaunay::set is 14%, this seems a bit slow

  - storing certs in edge is 10%

  - delete certificate is 5%

  - static certificates: how to do it for wide intervals

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

  typedef IndirectKernel Indirect_kernel;

  typedef Updatable_delaunay_triangulation_table_2<Indirect_kernel, Kinetic_kernel> Update_information;





  typedef typename Indirect_kernel::Geometric_point_2 Static_point_2;


  typedef typename Indirect_kernel::Swapable_container Points;
  typedef typename Kinetic_kernel::Point_2 Kinetic_point_2;
  typedef typename Kinetic_kernel::Motion_function Kinetic_coordinate;

  typedef typename Indirect_kernel::Current_coordinates IK_current_coordinates;



  typedef typename Kinetic_kernel::Certificate Exact_certificate;
  typedef typename Function_kernel::Root Exact_time;


  struct Refiner {
    struct CS: public CGAL::Kinetic::Ref_counted<CS> {
      CS(const Exact_certificate &c): cert_(c){}
       Exact_certificate cert_;
    };
    
    // enum State {INTERVAL, EXACT, INVALID};

    Refiner(Cert_tuple t,
	    typename Update_information::Certificate_function_const_pointer certf,
#ifdef NEWTON
	    typename Update_information::Certificate_derivitive_const_pointer certpf,
#endif
	    typename Update_information::Handle tbl): tuple_(t),
						      ui_(tbl),
						      certf_(certf)
#ifdef NEWTON
						    , certpf_(certpf)
#endif
{
    }
    Refiner() {
    }
    bool is_invalid() const {
      return tuple_.is_invalid();
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

#ifndef NEWTON
	std::pair<double,double> oiv=iv;
	typename Update_information::Certificate_sign_at<5> se(certf_);
	ui_->isolate_failure(se, sep,
			     iv.first, iv.second, false, 0, 20, iv);
#ifndef NDEBUG
	if ( oiv ==iv) {
	  std::cout << "Can't subdivide more in " << iv.first << " " << iv.second << std::endl;
	}
#endif
	return oiv != iv;
#else
	INT niv=ui_->Newton_refine(certf_, certpf_, INT(iv));
	bool ret= (niv.inf() != iv.first || niv.sup() != iv.second);
	iv= std::pair<double,double>(niv.inf(), niv.sup());
	return ret;
#endif
      } else {
	// compute exact;
	CGAL_assertion(!has_exact_root());
	//++comparison_certificates_;
	
	ensure_exact_root(iv);
	iv= CGAL::to_interval(exact_root());
	return true;
      }
    }

    const Exact_root& 
    exact_root() const {
      return cert_->cert_.failure_time();
    }

    const Exact_root& 
    exact_root(std::pair<double,double> &iv) const {
      ensure_exact_root(iv);
      return exact_root();
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
      if (tuple_ == Cert_tuple()) return false;
      return tuple_== o.tuple_;
    }

    bool has_exact_root() const {
      return cert_!= typename CS::Handle() ;
    }
    
    CGAL::Sign sign_at(INT t) const {
      return Update_information::sign_at(certf_, t);
    }

    void ensure_exact_root(std::pair<double,double> &iv) const {
      if (!has_exact_root()) {
	if (ui_ != typename Update_information::Handle()) {
	  ++stat_exact_certificate_functions_;
	  ++stat_exact_certificate_functions_from_compare_;
	  CGAL_UD_DEBUG("Generating exact with interval " << iv.first << " " << iv.second << std::endl);
	  cert_= new CS(ui_->compute_exact_failure_time(tuple_, iv.first));
	  CGAL_UD_DEBUG("Got " << cert_->cert_.failure_time() << std::endl);
	  CGAL_assertion(check_.failure_time() == cert_->cert_.failure_time());
	  iv= CGAL::to_interval(exact_root());
	} else {
	  double cs[2];
	  cs[0]=iv.first;
	  cs[1]=-1.0;
	  typename Exact_certificate::Function f(cs, cs+2);
	  // hack
	  cert_= new CS(Exact_certificate(f, KK_function_kernel(), -1, 2));
	  CGAL_UD_DEBUG(f << ": " << exact_root()  << std::endl); 
	  CGAL_assertion(exact_root() == Exact_time(iv.first));
#ifndef NDEBUG
	  check_= cert_->cert_;
#endif
	}
      }
    }
    void compute_exact_root(const Exact_time &et) const {
      if (!has_exact_root()) {
	++stat_exact_certificate_functions_;
	CGAL_UD_DEBUG("Generating exact with root " << et << std::endl);
	CGAL_precondition(ui_ != typename Update_information::Handle());
	//if (ui_ != ) {
	cert_= new CS(ui_->compute_exact_failure_time(tuple_, et));
	CGAL_assertion(check_.failure_time() == exact_root());
      }
    }

    const Exact_certificate& exact_certificate() const {
      return cert_->cert_;
    }

    std::pair<double,double> interval_from_exact() const {
      CGAL_precondition(has_exact_root());
      if (exact_root() >= Exact_time(1)) {
	return std::make_pair(1.0, -1.0);
      } else {
	std::pair<double,double> ip= CGAL::to_interval(exact_root());
	return ip;
      }
    }

    void set_exact_certificate(const Exact_certificate& ec) {
      CGAL_precondition(!has_exact_root());
      cert_= new CS(ec);
      CGAL_postcondition(exact_root() == check_.failure_time());
    }

    typename Update_information::Certificate_function_const_pointer certificate_function() const {
      return certf_;
    }

#ifdef NEWTON
    typename Update_information::Certificate_derivitive_const_pointer certificate_derivitive() const {
      return certpf_;
    }
#endif

    Cert_tuple tuple_;
    typename Update_information::Handle ui_;
    mutable typename CS::Handle cert_;
    typename Update_information::Certificate_function certf_;
#ifdef NEWTON
    typename Update_information::Certificate_derivitive certpf_;
#endif
#ifndef NDEBUG
    mutable Exact_certificate check_;
#endif
  };








  struct Simulation_traits {
    typedef SimTraits_base P;
    
    struct Sillyness {
      typedef Kinetic_kernel KK;
      //typedef Active_points_2_table APT;
    };
    typedef typename Sillyness::KK Kinetic_kernel;
    
    typedef Update_information Active_points_2_table;
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
    typedef typename Kinetic::Two_list_pointer_event_queue<Simulator_traits, false, 2> Queue;
    typedef typename Kinetic::Default_simulator<Simulator_traits, Queue > Simulator;

    typename Simulator::Handle simulator_handle() {
      return sim_;
    }


    typename Simulator::Handle simulator_handle() const {
      return sim_;
    }


    typedef typename Simulator::Time Time;


    Simulation_traits(Indirect_kernel ik, 
		      double lb,
		      double ub): sim_(new Simulator(lb, ub)){
      ap_=new Active_points_2_table(ik, kinetic_kernel_object());
    }

    Simulation_traits(){}

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

    Cert_tuple tuple(Point_key ks[4]) const {
      return Cert_tuple(ks);
    }
    
#ifndef NDEBUG
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
#endif
    
    void point_changed(Point_key k){
      //ui()->point_changed(k);
    }

    const Time &current_time() const {
      return P::simulator_handle()->current_time();
    }


    Certificate_pair null_pair() const {
      return std::make_pair(Time(2.0), Certificate_data());
    }

    Certificate_pair return_pair(Time rt) const {
      return Certificate_pair(rt, Certificate_data());
    }
    

    double rational_current_time() const {
      return CGAL::to_interval(current_time()).first;
    }


    Traits(Simulation_traits st): P(st){}
    //soc_(tr.kinetic_kernel_object().positive_side_of_oriented_circle_2_object()){}
    Update_information* ui() {
      return P::active_points_2_table_handle().get();
    }
 

    bool hull_certificate_failure_time(typename Default_traits::Edge, Point_key [3], 
				       Time, Certificate_data) {
      return false;
    }



    bool internal_certificate_failure_time(typename Default_traits::Edge e, Point_key pks[4],
					   Time &rett, Certificate_data) {
      ++stat_certificate_computations_;

      Cert_tuple ct= tuple(pks);
      //bool check_if_point_is_activated;
      double end_time= 1; //i()->next_activation();
      //ui()->set_next_activation(1.0);
      CGAL_UD_DEBUG("Current time is " << P::simulator_handle()->current_time() << std::endl);

#ifndef MOVE_ALL
      bool found=false;
      for (unsigned int i=0; i< 4; ++i) {
	if (ui()->is_activating(ct[i])) {
	  found=true;
	  break;
	}
      }
      if (found) {
	CGAL_UD_DEBUG("Shortening end time to " << ui()->next_activation() << std::endl);
	end_time= ui()->next_activation();
      }
#endif

#ifndef NDEBUG
      Exact_time check_failure_time;
      Exact_certificate check_cert;
      {
	Exact_time ect;
	if (P::simulator_handle()->current_time().refiner().is_invalid()) {
	  ect= Exact_time(CGAL::to_interval(P::simulator_handle()->current_time()).first);
	} else {
	  ect= P::simulator_handle()->current_time().refiner().check_.failure_time();
	}
	
	check_cert= ui()->in_circle_object()(ui()->exact_point(ct[0]),
					     ui()->exact_point(ct[1]),
					     ui()->exact_point(ct[2]), 
					     ui()->exact_point(ct[3]),
					     ect, 1);
	check_failure_time= check_cert.failure_time();

	if (check_failure_time > end_time) check_failure_time= 2;
      }
      //CGAL_UD_DEBUG("Computing failure time for " << ct << std::endl);
      CGAL_UD_DEBUG("True failure time is " << check_failure_time<< std::endl);
      //CGAL_UD_DEBUG("Exact cert is " << check_cert << std::endl);
#endif


      double bt;// = ui()->start_time();
      //if (bt <0) {
      bt= CGAL::to_interval(current_time()).first;
	//}
      
#ifdef MOVE_ALL
      if (!ui()->can_fail(ct, bt)) {
        CGAL_postcondition(check_failure_time > 1);
        return false;
      }

#else
      if (!ui()->can_fail(ct, bt, end_time)) {
        CGAL_postcondition(check_failure_time > 1);
        return false;
      }
#endif
      
      typename Update_information::Certificate_function cf;
#ifdef NEWTON
      typename Update_information::Certificate_derivitive cfp;
#endif
      INT ft;
      typename Update_information::Isolate_result isc;
#ifndef NDEBUG
#ifndef MOVE_ALL
      if (!(ui()->is_active(ct[0]) || ui()->is_active(ct[1])
	    || ui()->is_active(ct[2]) || ui()->is_active(ct[3]))) {
	isc= Update_information::NO_FAILURE;
      } else 
#endif
#endif
        {
	  CGAL::Protect_FPU_rounding<true> prot;
	  
	  ++stat_interval_certificate_functions_;
	  ui()->certificate_function(ct[0], ct[1], ct[2], ct[3], cf);
	  //typename Update_information::Certificate_evaluator se(cf);
	  
#ifndef NEWTON
	  isc= ui()->isolate_failure(cf, 
				     bt,
				     end_time,
				     true,
				     ft);
#else
	  cfp= cf.prime();
	  typename Update_information::Certificate_acceleration cfpp= cfp.prime();
	  isc= ui()->Newton_isolate(cf, cfp, cfpp, 
				    bt, end_time, ft);
	  CGAL_UD_DEBUG("Newton isolate got " << isc << " with interval " << ft << std::endl);
	  if (isc == Update_information::CERTAIN_FAILURE) {
	    ft= ui()->Newton_refine(cf,cfp, ft);
	    CGAL_UD_DEBUG("After refinement got " << ft << std::endl);
	  }
	}
#endif

      if (isc == Update_information::NO_FAILURE) { 
        CGAL_UD_DEBUG("No root" << std::endl << std::endl);
        CGAL_assertion(check_failure_time >=1);
        return false;
      }
#ifdef NEWTON
      rett= Time(ft.inf(), ft.sup(), Refiner(ct, cf,cfp, ui()));
#else
      rett= Time(ft.inf(), ft.sup(), Refiner(ct, cf, ui()));
#endif
#ifndef NDEBUG
      rett.refiner().check_= check_cert;
#endif

      /* Check if the curren time is exact, if so check if I am positive at it*/


      if (isc == Update_information::CERTAIN_FAILURE) {
	if (CGAL::compare(rett, P::simulator_handle()->current_time()) != CGAL::LARGER) {
	  // interval overlaps current
	  //std::pair<double,double> cti= CGAL::to_interval(P::simulator_handle()->current_time());
	  //P::simulator_handle()->current_time().refiner().ensure_exact_root(cti);
	  const Exact_time &ect= P::simulator_handle()->current_time().exact_root();
	  rett.refiner().compute_exact_root(ect);
	  ++stat_exact_certificate_functions_from_compare_curt_;
	  CGAL_UD_DEBUG("Can't separate root from cur " << std::endl);
	  std::cout << "Have to generate exact certificate" << std::endl;
	  CGAL_UD_DEBUG("Exact" << rett.refiner().exact_root() <<  std::endl);
	  rett.set_interval(CGAL::to_interval(rett.refiner().exact_root()));
	} else {
#ifndef NDEBUG
	  CGAL::Protect_FPU_rounding<true> prot;
	  double owid;
	  INT nft= ui()->Newton_refine(cf,cfp, ft);
	  do {
	    owid= nft.sup()- nft.inf();
	    nft= ui()->Newton_refine(cf,cfp, nft);
	  } while (owid  > 1.01 * (nft.sup() -nft.inf()));
	  stat_total_interval_width_+= owid;
	  stat_number_of_intervals_+= 1;
#endif
	}
	CGAL_UD_DEBUG("Returning isolated " << rett << std::endl);
	CGAL_assertion(Exact_time(CGAL::to_interval(rett).first) <= check_failure_time);
	CGAL_assertion(Exact_time(CGAL::to_interval(rett).second) >= check_failure_time);
	return true;
      } else {
	//if (!has_exact_failure_time(ct)) {
	//	}
	CGAL_UD_DEBUG("Possible root " <<  std::endl);
	//std::pair<double,double> cti=CGAL::to_interval(P::simulator_handle()->current_time());
	//P::simulator_handle()->current_time().refiner().ensure_exact_root(cti);
	const Exact_time &ect= P::simulator_handle()->current_time().exact_root();
	rett.refiner().compute_exact_root(ect);
	CGAL_UD_DEBUG(std::cout << "Have to generate exact certificate" << std::endl);
	CGAL_UD_DEBUG("Exact root is" << rett.refiner().exact_root() <<  std::endl);
	CGAL_UD_DEBUG("Exact cert is " << rett.refiner().exact_certificate() << std::endl);
	CGAL_UD_DEBUG("Initial points are " << std:: endl);
#ifndef NDEBUG
	typename Indirect_kernel::Current_coordinates ic= ui()->initial_coordinates_object();
	typename Indirect_kernel::Current_coordinates fc= ui()->final_coordinates_object();
	
#endif
	{
	  CGAL::Protect_FPU_rounding<true> prot;
	CGAL_UD_DEBUG(ic(ct[0]) << ": "
		      << ic(ct[1]) << ": " 
		      << ic(ct[2]) << ": " 
		      << ic(ct[3]) << std::endl);
	CGAL_UD_DEBUG("Final points are " << std:: endl);
	CGAL_UD_DEBUG(fc(ct[0]) << ": "
		      << fc(ct[1]) << ": " 
		      << fc(ct[2]) << ": " 
		      << fc(ct[3]) << std::endl);
	CGAL_UD_DEBUG("Initial incircle is " 
		      << ui()->eval_incircle(INT(ic(ct[0])[0]),  INT(ic(ct[0])[1]),
					     INT(ic(ct[1])[0]),  INT(ic(ct[1])[1]),
					     INT(ic(ct[2])[0]),  INT(ic(ct[2])[1]),
					     INT(ic(ct[3])[0]),  INT(ic(ct[3])[1])) << std::endl);
	CGAL_UD_DEBUG("Final incircle is " 
		      << ui()->eval_incircle(INT(fc(ct[0])[0]),  INT(fc(ct[0])[1]),
					     INT(fc(ct[1])[0]),  INT(fc(ct[1])[1]),
					     INT(fc(ct[2])[0]),  INT(fc(ct[2])[1]),
					     INT(fc(ct[3])[0]),  INT(fc(ct[3])[1])) << std::endl);
	}
					     
	if (rett.refiner().exact_root() >= 1) {
	  CGAL_assertion(check_failure_time >=1);
	  CGAL_UD_DEBUG("Phantom root " << std::endl);
	  CGAL_postcondition(check_failure_time > 1);
	  ++stat_unfailing_exact_certificate_functions_;
	  return false;
	} else {
	  rett.set_interval(CGAL::to_interval(rett.refiner().exact_root()));
	  CGAL_UD_DEBUG("Returning exact " << rett << std::endl);
	  return true;
	}
      }
    }

    bool certificate_failure_time(
#ifndef NDEBUG
				  typename Default_traits::Edge e, 
#else
				  typename Default_traits::Edge, 
#endif
                                  Certificate_data , Time &rett, Certificate_data) {
      ++stat_certificate_advances_;
      const Time &curt= P::simulator_handle()->current_time();
      Cert_tuple ct= curt.refiner().tuple().opposite();
      CGAL_UD_DEBUG("Current time is " << P::simulator_handle()->current_time() << std::endl);

#ifndef NDEBUG
      Cert_tuple check_ct= tuple(e);
#endif
      CGAL_assertion(ct == check_ct);

      double end_time= 1; //ui()->next_activation();

#ifndef MOVE_ALL
      ui()->set_next_activation(1.0);
#endif
 
#ifndef NDEBUG
      Exact_time check_failure_time;
      Exact_certificate check_cert;
      {
	check_cert= curt.refiner().check_;
	check_cert.pop_failure_time();

	check_failure_time= check_cert.failure_time();
	if (check_failure_time > end_time) check_failure_time= 2;
	//GAL_UD_DEBUG("Exact time is " << check_cert.failure_time() << std::endl);
	CGAL_UD_DEBUG("True failure time is " << check_cert.failure_time() << std::endl);
      }
#endif


      // this depends on being the last event of the batch whose time is computed
      const Time& net= P::simulator_handle()->next_event_time();

     

      double lb= CGAL::to_interval(net).first;
      CGAL_UD_DEBUG("Advancing certificate for " << ct << " to " << lb
		    << std::endl);

#ifndef MOVE_ALL
      bool found=false;
      for (unsigned int i=0; i< 4; ++i) {
	if (ui()->is_activating(ct[i])) {
	  found=true;
	  break;
	}
      }
      if (found) {
	CGAL_UD_DEBUG("Shortening end time to " << ui()->next_activation() << std::endl);
	end_time= ui()->next_activation();
	if (end_time <= lb) {
	  CGAL_UD_DEBUG("And the interval is empty " << std::endl);
	  return false;
	}
      }
#endif
    
      typename Update_information::Certificate_function cf(curt.refiner().certificate_function(), true);
         
      
      bool wrap_this_into_newton;
      // done twiceish
      CGAL::Protect_FPU_rounding<true> prot;
      ++stat_point_predicate_evaluations_;
      //CGAL::Sign sn= curt.refiner().sign_at(INT(lb,lb));
      {
	INT lbv= cf(lb);
	
	//ui()->sign_at(ct[0],ct[1],ct[2],ct[3], INT(lb, lb));
	if (lbv.inf()<=0) {
	  CGAL_UD_DEBUG("Not positive at " << lb << std::endl);
	  lb= CGAL::to_interval(curt).first;
	} 
      }

      {
	++stat_interval_predicate_evaluations_;
	//CGAL::Sign sn= curt.refiner().sign_at(INT(lb,lb));
	INT lbv= cf(INT(lb, end_time));
	if (lbv.inf()>0) {
	  CGAL_UD_DEBUG("No root in interval." << std::endl);
	  CGAL_postcondition(check_failure_time > 1);
	  return false;
	} 
      }

  
#ifdef NEWTON
      typename Update_information::Certificate_derivitive cfp(curt.refiner().certificate_derivitive(), true);
#endif
      INT ft;
      typename Update_information::Isolate_result isc;
      {
	CGAL::Protect_FPU_rounding<true> prot;
	//ui()->certificate_function(ct[0], ct[1], ct[2], ct[3], cf);
;
	//typename Update_information::Certificate_evaluator se(cf);
#ifndef NEWTON
	isc
	  = ui()->isolate_failure(cf, lb, end_time, pos_start, ft);
#else
	
	typename Update_information::Certificate_acceleration cfpp(cfp.prime());
	isc= ui()->Newton_isolate(cf, cfp, cfpp, lb, end_time, ft);
	if (isc == Update_information::CERTAIN_FAILURE) {
	  ft= ui()->Newton_refine(cf, cfp, ft);
	}
	CGAL_UD_DEBUG("Newton isolate got " << isc << " with interval " << ft << std::endl);
#endif
      }
      
  

      if (isc == Update_information::NO_FAILURE) { 
	CGAL_UD_DEBUG("No root there." << std::endl);
	CGAL_postcondition(check_failure_time > 1);
	return false;
      }

#ifdef NEWTON
      rett=Time(ft.inf(), ft.sup(), Refiner(ct, cf, cfp, ui()));
#else
      rett=Time(ft.inf(), ft.sup(), Refiner(ct, cf, ui()));
#endif
#ifndef NDEBUG
      rett.refiner().check_= check_cert;
#endif
      if (isc == Update_information::CERTAIN_FAILURE
	  && ft.inf() <= CGAL::to_interval(curt).second
	  || isc == Update_information::POSSIBLE_FAILURE
	  || curt.refiner().has_exact_root()) {
	
	curt.exact_root();
	CGAL_UD_DEBUG("Have to generate exact certificate" << std::endl);
	Exact_certificate ec= curt.refiner().exact_certificate();
	CGAL_UD_DEBUG("Curt is now " << curt << std::endl);
	CGAL_UD_DEBUG("Root was " << ec.failure_time() << std::endl);
	ec.pop_failure_time();
	CGAL_UD_DEBUG("Root is " << ec.failure_time() << std::endl);
	CGAL_UD_DEBUG("Exact cert is " << ec << std::endl);
	++stat_exact_certificate_functions_from_advance_;
	if (ec.failure_time() > 1) {
	  CGAL_UD_DEBUG("Phantom root " << std::endl);
	  CGAL_postcondition(check_failure_time > 1);
	  return false;
	}
	ft= CGAL::to_interval(curt.refiner().exact_root());
	rett.set_interval(CGAL::to_interval(ft));
	rett.refiner().set_exact_certificate(ec);
	return true;
      } else {
	CGAL_assertion(Exact_time(CGAL::to_interval(rett).first) <= check_failure_time);
	CGAL_assertion(Exact_time(CGAL::to_interval(rett).second) >= check_failure_time);
	CGAL_UD_DEBUG("Returning " << rett << std::endl);
	return true;
      }
    }

  

    CGAL::Comparison_result compare_concurrent(Event_key a,
					       Edge,
					       Event_key b,
					       Edge) const {
      return CGAL::compare(a,b);
    }
  };




















  struct Visitor: public CGAL::Kinetic::Delaunay_triangulation_visitor_base_2 {
  
    
    Visitor(Simulation_traits tr): tr_(tr) {
    } 

    Update_information* ui() {
      return tr_.active_points_2_table_handle();
    }

    const Update_information* ui() const {
      return tr_.active_points_2_table_handle();
    }

#ifndef MOVE_ALL
    bool test_and_add(Edge e, std::vector<Point_key> &active) const {
      ++stat_point_predicate_evaluations_;
      if (!compute_ok(e, ui()->final_kernel_object())) {
        add(e.first->vertex(0)->point(), active);
	add(e.first->vertex(1)->point(), active);
	add(e.first->vertex(2)->point(), active);
	add(TDS_helper::mirror_vertex(e)->point(), active);
	return true;
      } else return false;
    }

    bool test_and_add_one(Edge e, std::vector<Point_key> &active) const {
      ++stat_point_predicate_evaluations_;
      if (!compute_ok(e, ui()->final_kernel_object())) {
	//add(e.first->vertex(0)->point(), active);
	//add(e.first->vertex(1)->point(), active);
	//add(e.first->vertex(2)->point(), active);
        add(TDS_helper::mirror_vertex(e)->point(), active);
        return true;
      } else return false;
    }
#endif

    template <class TDS> 
    void initialize_events(const TDS &triangulation,
                           Indirect_kernel fk) {
      stat_number_of_edges_=triangulation.tds().number_of_edges();
      ui()->set_final_kernel(fk);
#ifndef MOVE_ALL
      std::vector<Point_key> active;
      for (typename TDS::Finite_edges_iterator it= triangulation.finite_edges_begin(); 
           it != triangulation.finite_edges_end(); ++it){
        if (test_and_add(*it, active)) {
          ++stat_number_of_bad_edges_;
        }
      }
      
      Interpolate_event ev(tr_,
                           ui(),
                           ui()->initial_coordinates_object(),
                           ui()->final_coordinates_object(), 
                           active.begin(), 
                           active.end());
      ++stat_number_of_interpolations_;
      if (!ev.empty()) {
        INT iat= CGAL::to_interval(ev.time());
        tr_.simulator_handle()->new_event(ev.time(), ev);
      }
#endif
    }


  
    void after_flip(Edge e) {
      //++num_events_;
      // schedule a bulk set event for next rational time
      std::vector<Point_key> active;
      //test_and_add(e, active);
      //CGAL_precondition(ui_->is_active(e.first->vertex((e.second+1)%3)->point()));
      //CGAL_precondition(ui_->is_active(e.first->vertex((e.second+2)%3)->point()));
      //CGAL_precondition(ui_->is_active(e.first->vertex(e.second)->point()));
      Edge em= TDS_helper::mirror_edge(e);
#ifndef MOVE_ALL
      add(e.first->vertex(0)->point(), active);
      add(e.first->vertex(1)->point(), active);
      add(e.first->vertex(2)->point(), active);

      //CGAL_precondition(ui_->is_active(em.first->vertex(em.second)->point()));
      add(em.first->vertex(em.second)->point(), active);
      //test_and_add_one(Edge(e.first, (e.second+1)%3), active);
      //test_and_add_one(Edge(e.first, (e.second+2)%3), active);
      //test_and_add_one(Edge(em.first, (em.second+1)%3), active);
      //test_and_add_one(Edge(em.first, (em.second+2)%3), active);
      if (!active.empty()) {
	Interpolate_event ev(tr_, 
			     ui(),
			     ui()->initial_coordinates_object(),
			     ui()->final_coordinates_object(), 
			     active.begin(), active.end());
	if (!ev.empty()){
	  ++stat_number_of_interpolations_;
	  tr_.simulator_handle()->new_event(ev.time(), ev);
	  //active_.set(active.begin(), active.end());
	  //++num_interpolations_;
	  /*for (unsigned int i=0; i< active.size(); ++i){
	    ui_->activate(active[i], true);
	    }*/
          ui()->set_next_activation(CGAL::to_interval(ev.time()).second);
        }
      }
#endif
      /*
      bool ok=true;
      double skip_to= CGAL::to_interval(tr_.simulator_handle()->next_event_time()).first;
      {
	CGAL::Protect_FPU_rounding<true> prot;
	CGAL::Sign sn = tr_.simulator_handle()->current_time().refiner().sign_at(INT(skip_to));

	if (sn == CGAL::POSITIVE) {
	  for (int i=0; i< 3; ++i) {
	    if (i != e.second) {
	      CGAL::Sign sn= ui()->sign_at(e.first->vertex(0)->point(), 
					  e.first->vertex(1)->point(),
					  e.first->vertex(2)->point(), 
					  e.first->mirror_vertex(i)->point(),
					  INT(skip_to));
	      if (sn != CGAL::POSITIVE) {
		ok=false;
		break;
	      }
	    }
	    if (i != em.second) {
	      CGAL::Sign sn= ui()->sign_at(em.first->vertex(0)->point(),
					  em.first->vertex(1)->point(),
					  em.first->vertex(2)->point(),
					  em.first->mirror_vertex(i)->point(),
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
	ui()->set_start_time(skip_to);
      } else {
	ui()->set_start_time(-1);
      }
      */
    }

  

   
    
#ifndef MOVE_ALL
    bool is_active(Point_key k) const {
      return ui()->is_active(k);
    }

    void add( Point_key k, std::vector<Point_key> &active)  const {
      if (!is_active(k)) {
        active.push_back(k);
      }
    }
#endif

    void reset() {
      ui()->reset();
    }
    
    Static_point_2 initial(Point_key pk) const {
      return ui()->ik_.current_coordinates_object()(pk);
    }
    Static_point_2 final(Point_key pk) const {
      return ui()->fk_.current_coordinates_object()(pk);
    }

    void stats_clear() {
      ui()->clear_stats();
    }
    void stats_write(std::ostream &out) {
      ui()->write_stats(out);
    }

    Simulation_traits tr_;  
  };


  typedef CGAL::Kinetic::Delaunay_triangulation_2<Simulation_traits, Visitor, Triangulation, Traits> KDel;


#ifndef MOVE_ALL


  static Kinetic_point_2 interpolate_t1(NT time,
                                        Static_point_2 ip,
					Static_point_2 fp) {
    typedef typename Simulation_traits::Kinetic_kernel::Motion_function MF;
    typedef typename MF::NT NT;
    Protected_array<MF, 2> mf;
    for (unsigned int i=0; i< 2; ++i){
      Protected_array<NT,2> c;
      c[1]=(NT(ip[i])-NT(fp[i]))/NT(time-1);
      c[0]=NT(fp[i])-c[1];
      mf[i]=MF(c, c+2);
    }
    return Kinetic_point_2(mf[0], mf[1]);
  }


  static Kinetic_point_2 interpolate_12(Static_point_2 ip,
					Static_point_2 fp) {
    typedef typename Simulation_traits::Kinetic_kernel::Motion_function MF;
    Protected_array<MF, 2> mf;
    for (unsigned int i=0; i< 2; ++i){
      Protected_array<NT,2> c;
      c[1]=(fp[i]-ip[i]);
      c[0]=2*ip[i]-fp[i];
      mf[i]=MF(c.begin(), c.end());
    }
    return Kinetic_point_2(mf[0], mf[1]);
  }



  struct Interpolate_event: public CGAL::Kinetic::Free_event_base {
  
    
    template <class It>
    Interpolate_event(Simulation_traits tr,
		      typename Update_information::Handle ui,
		      IK_current_coordinates ic,
		      IK_current_coordinates fc, 
		      It b, It e): ui_(ui) {
      if (CGAL::to_interval(tr.simulator_handle()->current_time()).second != 0) {
	time_= nextafter(CGAL::to_interval(tr.simulator_handle()->current_time()).second,
			 std::numeric_limits<double>::max());
      } else {
	time_=0;
      }
      //}
      std::sort(b,e);
      It ne= std::unique(b,e);
      ui_->set_next_activation(time_);
      {
	CGAL::Protect_FPU_rounding<true> prot;
	for (It c=b; c!= ne; ++c){
	  if (ic(*c) != fc(*c)) {
	    motions_.push_back(*c);
	    ui_->preactivate(*c);
	  } else {
	    CGAL_UD_DEBUG("Point " << *c << " doesn't move." << std::endl);
	  }
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
	out << motions_[i] << " ";
      }
      return out;
    }

    bool empty() const {
      return motions_.empty();
    }

    void process() {
      ui_->set_next_activation(-1);
      INT it= CGAL::to_interval(time_);
      ui_->set_is_editing(Update_information::LOGGED);
      CGAL::Protect_FPU_rounding<true> prot;
      for (unsigned int i=0; i< motions_.size(); ++i) {

	ui_->activate(time_, motions_[i]);
      }
      ui_->set_is_editing(Update_information::NOT);
    }
    typename Update_information::Handle ui_;
    typename std::vector<Point_key> motions_;
    double time_;
  };


    

  struct Final_event: public CGAL::Kinetic::Free_event_base {
    typedef typename Simulation_traits::Active_points_2_table::Key Table_key;

    Final_event(Simulation_traits tr, typename KDel::Handle kdel):
      tr_(tr), kdel_(kdel){
    }

    
    std::ostream & write(std::ostream&out) const {
      out << "Final event ";
      return out;
    }

    void process() { 

      kdel_->write_stats(std::cout);
      kdel_->visitor().stats_write(std::cout);
      
      tr_.simulator_handle()->set_interval(1,2);
      tr_.active_points_2_table_handle()->set_is_editing(Update_information::LOGGED);
      
      for (typename Simulation_traits::Active_points_2_table::Key_iterator 
	     it = tr_.active_points_2_table_handle()->keys_begin(); 
	   it != tr_.active_points_2_table_handle()->keys_end(); ++it) {
	if (!kdel_->visitor().is_active(*it)) {	
  tr_.active_points_2_table_handle()->set(*it, interpolate_12(initial(*it),
								      final(*it)));
	} else {
	  Kinetic_point_2 np(ENT(final(*it).x()),
			     ENT(final(*it).y()));
	  tr_.active_points_2_table_handle()->set(*it,np);
						  
	}

      }

      tr_.active_points_2_table_handle()->set_is_editing(Update_information::NOT);
    }

    Static_point_2 initial(Point_key pk) const {
      return tr_.active_points_2_table_handle()->initial(pk);
    }
    Static_point_2 final(Point_key pk) const {
      return tr_.active_points_2_table_handle()->final(pk);
    }

    Simulation_traits tr_;
    typename KDel::Handle kdel_;
  };


 #endif
  


  template <class It> 
  Updatable_Delaunay_triangulation_2(It b, It e) {
    
    typename Indirect_kernel::Key_range rg= ik_.new_point_2s(b,e);
    tr_=  Simulation_traits(ik_, 0,0);
    Triangulation tr(ik_);

    tr.insert(rg.first, rg.second);
    
    Traits traits(tr_);
  
    
    kdel_= new KDel(traits, tr, Visitor(tr_));
    kdel_->clear_stats();
    kdel_->visitor().stats_clear();
    kdel_->set_neighbors_initialized(true);
#ifndef MOVE_ALL
    kdel_->set_has_certificates(true, true, true);
#endif
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
#ifndef MOVE_ALL
    tr_.simulator_handle()->new_final_event(Final_event(tr_,
                                                        kdel_));
#endif

    std::cout << "Green edges just flipped, grey edges will not flip until"
              << " their certificate changes and black edges will flip." << std::endl;
  
    qtsim->begin_event_loop();
    kdel_->visitor().stats_write(std::cout);
    ik_.swap(fk);
    audit();
  }

  Indirect_kernel set_up_update(const Points &pts) {
    kdel_->visitor().reset();

    kdel_->set_has_certificates(false, true);
    tr_.simulator_handle()->set_interval(0,1);

    Indirect_kernel fk;
    fk.new_point_2s(pts.begin(), pts.end());
    kdel_->visitor().initialize_events(kdel_->triangulation(), fk);
#ifdef MOVE_ALL
    kdel_->set_has_certificates(true, false, true);
#else
    kdel_->set_has_certificates(true, true, true);
#endif
    return fk;
  }

  void update_coordinates(const Points &pts) {
    Indirect_kernel fk=set_up_update(pts);
    tr_.simulator_handle()->set_current_time(1);
    ik_.swap(fk);
  }

  void write_statistics(std::ostream &out) const {
    kdel_->write_stats(out);
    kdel_->visitor().stats_write(out);
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

 



  static int run(int , char *[], int , int , 
                 int , std::string ifile, std::string ffile) {
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
  Indirect_kernel ik_;
  Simulation_traits tr_;
  typename KDel::Handle kdel_;
};

CGAL_END_NAMESPACE

#endif
#undef CGAL_UD_DEBUG
