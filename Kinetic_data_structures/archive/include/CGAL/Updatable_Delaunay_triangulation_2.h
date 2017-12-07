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


namespace CGAL {

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


  /*typedef CGAL::Gmpq NT;
  typedef CGAL::POLYNOMIAL::Polynomial<NT> Function;
  typedef CGAL::POLYNOMIAL::Sturm_root_stack_traits<Function> Root_stack_traits;
  typedef CGAL::POLYNOMIAL::Sturm_root_stack<Root_stack_traits> Root_stack;
  typedef CGAL::POLYNOMIAL::Kernel<Function, Root_stack> Function_kernel;*/
  typedef CGAL::POLYNOMIAL::CORE_kernel Function_kernel;
  typedef typename Function_kernel::Root_stack Root_stack;
  typedef typename Function_kernel::Function Function;
  typedef typename Function_kernel::FT NT;


  typedef typename CGAL::Kinetic::Handle_degeneracy_function_kernel<Function_kernel, true>  KK_function_kernel;
  typedef typename CGAL::Kinetic::Cartesian_kinetic_kernel<KK_function_kernel> Kinetic_kernel;

  typedef IndirectKernel Indirect_kernel;



  typedef Updatable_delaunay_triangulation_table_2<Indirect_kernel, Kinetic_kernel> Update_information;
  typedef typename Update_information::Simulator Simulator;
  typedef typename Simulator::Event_key Event_key;




  typedef typename Indirect_kernel::Geometric_point_2 Static_point_2;


  typedef typename Indirect_kernel::Swapable_container Points;
  typedef typename Kinetic_kernel::Point_2 Kinetic_point_2;
  typedef typename Kinetic_kernel::Motion_function Kinetic_coordinate;

  typedef typename Indirect_kernel::Current_coordinates IK_current_coordinates;



  typedef typename Kinetic_kernel::Certificate Exact_certificate;
  typedef typename Function_kernel::Root Exact_time;







  struct Simulation_traits {
    typedef SimTraits_base P;

    struct Sillyness {
      typedef Kinetic_kernel KK;
      typedef Simulator S;
      //typedef Active_points_2_table APT;
    };
    typedef typename Sillyness::KK Kinetic_kernel;
    typedef typename Sillyness::S Simulator;

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


  //  typedef typename Simulation_traits::Simulator::Event_key Event_key;





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

    void point_changed(Point_key ){
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


    double check_activation(const Cert_tuple &ct, double et) {
      bool found=false, activate=false;
      for (unsigned int i=0; i< 4; ++i) {
	if (!ui()->is_active(ct[i])) {
	  if (!found) {
	    if (!ui()->ok_at_end(ct[0], ct[1], ct[2], ct[3])) {
	      activate=true;
	      CGAL_UD_DEBUG("Shortening end time to " << ui()->activation_time(P::simulator_handle()) << std::endl);
	      et= ui()->activation_time(P::simulator_handle());
	    }
	    found=true;
	  }
	  if (activate) {
	    CGAL_UD_DEBUG("Will activate " << ct[i] << std::endl);
	    ui()->schedule_activation(ct[i], P::simulator_handle());

	  }
	}
      }
      return et;
    }

    bool internal_certificate_failure_time(typename Default_traits::Edge , Point_key pks[4],
					   Time &rett, Certificate_data) {
      ++stat_certificate_computations_;

      Cert_tuple ct= tuple(pks);
      //bool check_if_point_is_activated;
      double end_time= 1; //i()->next_activation();
      //ui()->set_next_activation(1.0);
      CGAL_UD_DEBUG("Current time is " << P::simulator_handle()->current_time() << std::endl);

#ifndef NDEBUG
      Exact_time check_failure_time;
      Exact_certificate check_cert;
      {
	Exact_time ect;
	if (current_time().refiner().is_invalid()) {
	  ect= Exact_time(CGAL::to_interval(current_time()).first);
	} else {
	  ect= current_time().refiner().check_.failure_time();
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
      if (!disable_filter_1_ && !ui()->can_fail(ct, bt)) {
        CGAL_postcondition(check_failure_time > 1);
        return false;
      }

#else
      if (ui()->is_activating()) {
	if (!ui()->mixed_ok_at_end(ct[0], ct[1], ct[2], ct[3])) {
	  bool found=false;
	  for (unsigned int i=0; i< 4; ++i) {
	    if (!ui()->is_active(ct[i])) {
	      found=true;
	      CGAL_UD_DEBUG("Will activate " << ct[i] << std::endl);
	      ui()->schedule_activation(ct[i], P::simulator_handle());
	    }
	  }
	  if (found) {
	    CGAL_UD_DEBUG("Will activate " << std::endl);
	    return false;
	  }
	}
      } else {
	end_time= check_activation(ct, end_time);
      }

      if (!disable_filter_1_ && (end_time == bt || !ui()->can_fail(ct, bt, end_time))) {
	CGAL_postcondition(check_failure_time > 1);
	return false;
      }


#endif

      typename Update_information::Certificate_function cf;
      typename Update_information::Certificate_derivitive cfp;

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
	  cfp= cf.prime();
	  if (!disable_filter_2_) {
	    bt= ui()->derivitive_filter(cf, cfp, bt, end_time);
	    if (bt>= end_time) {
	      CGAL_UD_DEBUG("No root with deriv." << std::endl);
	      return false;
	    }
	  }

	  //typename Update_information::Certificate_evaluator se(cf);

	  typename Update_information::Certificate_acceleration cfpp= cfp.prime();
	  if (!disable_filter_3_) {
	    isc= ui()->Newton_isolate(cf, cfp, cfpp,
				      bt >  CGAL::to_interval(current_time()).second,
				      bt, end_time, ft);
	  } else {
	    isc = Update_information::POSSIBLE_FAILURE;
	    ft = INT(bt, end_time);
	  }
	  CGAL_UD_DEBUG("Newton isolate got " << isc << " with interval " << ft << std::endl);
	  if (isc == Update_information::CERTAIN_FAILURE) {
	    ft= ui()->Newton_refine(cf,cfp, ft);
	    CGAL_UD_DEBUG("After refinement got " << ft << std::endl);
	  } else if (isc == Update_information::POSSIBLE_FAILURE) {
	    if (!current_time().is_point()) current_time().refine();
	      //current_time().refine();
	    CGAL_UD_DEBUG("Refining current time to " << current_time() << std::endl);
	    bt= CGAL::to_interval(current_time()).first;
	    if (!disable_filter_3_) {
	      isc= ui()->Newton_isolate(cf, cfp, cfpp,
					bt >  CGAL::to_interval(current_time()).second,
					bt, end_time, ft);
	    }
	    CGAL_UD_DEBUG("Newton isolate got " << isc << " with interval " << ft << std::endl);
	    if (isc == Update_information::CERTAIN_FAILURE) {
	      ft= ui()->Newton_refine(cf,cfp, ft);
	      CGAL_UD_DEBUG("After refinement got " << ft << std::endl);
	    }
	  }
	}

      if (isc == Update_information::NO_FAILURE) {
        CGAL_UD_DEBUG("No root" << std::endl << std::endl);
        CGAL_assertion(check_failure_time >=1);
        return false;
      }
      rett= Time(ft.inf(), ft.sup(), typename Update_information::Refiner(ct, cf,cfp, ui()));
#ifndef NDEBUG
      rett.refiner().check_= check_cert;
#endif

      /* Check if the curren time is exact, if so check if I am positive at it*/


      if (isc == Update_information::CERTAIN_FAILURE
	  && CGAL::compare(rett, current_time()) == CGAL::LARGER) {

#ifndef NDEBUG
	CGAL_UD_DEBUG("BEGIN NARROW" << std::endl);
	CGAL::Protect_FPU_rounding<true> prot;
	double owid;
	INT nft= ui()->Newton_refine(cf,cfp, ft);
	do {
	  owid= nft.sup()- nft.inf();
	  nft= ui()->Newton_refine(cf,cfp, nft);
	} while (owid  > 1.01 * (nft.sup() -nft.inf()));
	stat_total_interval_width_+= owid;
	stat_number_of_intervals_+= 1;
	CGAL_UD_DEBUG("END NARROW" << std::endl);
#endif
	CGAL_UD_DEBUG("Returning isolated " << rett << std::endl);
	CGAL_assertion(Exact_time(CGAL::to_interval(rett).first) <= check_failure_time);
	CGAL_assertion(Exact_time(CGAL::to_interval(rett).second) >= check_failure_time);
	return true;
      } else {
	CGAL_UD_DEBUG("Possible root " <<  std::endl);

	rett.refiner().compute_exact_root(CGAL::to_interval(ft));

	CGAL_UD_DEBUG(std::cout << "Have to generate exact certificate" << std::endl);
	CGAL_UD_DEBUG("Exact root is" << rett.refiner().exact_root() <<  std::endl);
	CGAL_UD_DEBUG("Exact cert is " << rett.refiner().exact_certificate() << std::endl);
	if (rett.exact_root()  <= 1) {
	  const Exact_time &ect= current_time().exact_root();
	  rett.refiner().advance_exact_root(ect);
	  if (rett.exact_root()  <= 1) {
	    rett.set_interval(CGAL::to_interval(rett.refiner().exact_root()));
	    CGAL_postcondition(CGAL::compare(rett.exact_root(), check_failure_time) == 0);
	    CGAL_UD_DEBUG("Returning exact " << rett << std::endl);
	    return true;
	  } else {
	    CGAL_assertion(check_failure_time >=1);
	    CGAL_UD_DEBUG("Very phantom root " << std::endl);
	    CGAL_postcondition(check_failure_time > 1);
	    ++stat_unfailing_exact_certificate_functions_;
	    return false;
	  }
	} else {
	  CGAL_assertion(check_failure_time >=1);
	  CGAL_UD_DEBUG("Phantom root " << std::endl);
	  CGAL_postcondition(check_failure_time > 1);
	  ++stat_unfailing_exact_certificate_functions_;
	  return false;
	}

#if 0
	CGAL_UD_DEBUG("Initial points are " << std:: endl);

	typename Indirect_kernel::Current_coordinates ic= ui()->initial_coordinates_object();
	typename Indirect_kernel::Current_coordinates fc= ui()->final_coordinates_object();

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
#endif
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
      const Time &curt= current_time();
      Cert_tuple ct= curt.refiner().tuple().opposite();
      CGAL_UD_DEBUG("Current time is " << current_time() << std::endl);

#ifndef NDEBUG
      Cert_tuple check_ct= tuple(e);
#endif
      CGAL_assertion(ct == check_ct);

      double end_time= 1; //ui()->next_activation();


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
      //end_time= check_activation(ct, end_time);
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

      typename Update_information::Certificate_derivitive cfp(curt.refiner().certificate_derivitive(), true);
      if (!disable_filter_2_) {
	lb= ui()->derivitive_filter(cf, cfp, lb, end_time);
	if (lb>= end_time) {
	  CGAL_UD_DEBUG("No root with deriv." << std::endl);
	  return false;
	}
      }


      INT ft;
      typename Update_information::Isolate_result isc;
      {
	CGAL::Protect_FPU_rounding<true> prot;
	//ui()->certificate_function(ct[0], ct[1], ct[2], ct[3], cf);
	;
	//typename Update_information::Certificate_evaluator se(cf);

	typename Update_information::Certificate_acceleration cfpp(cfp.prime());
	if (!disable_filter_3_) {
	  isc= ui()->Newton_isolate(cf, cfp, cfpp,
				    lb >  CGAL::to_interval(current_time()).second,
				    lb,  end_time, ft);
	} else {
	  isc= Update_information::POSSIBLE_FAILURE;
	  ft= INT(lb, end_time);
	}
	if (isc == Update_information::CERTAIN_FAILURE) {
	  ft= ui()->Newton_refine(cf, cfp, ft);
	}
	CGAL_UD_DEBUG("Newton isolate got " << isc << " with interval " << ft << std::endl);
      }



      if (isc == Update_information::NO_FAILURE) {
	CGAL_UD_DEBUG("No root there." << std::endl);
	CGAL_postcondition(check_failure_time > 1);
	return false;
      }

      rett=Time(ft.inf(), ft.sup(), typename Update_information::Refiner(ct, cf, cfp, ui()));

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
#if 0
    template <int PASS, class VH>
    bool test_and_add(Edge e,  std::set<Point_key> &active, std::vector<VH> &new_active) const {
      Protected_array<Point_key, 4> pts;
      Protected_array<VH, 4> vhs;
      pts[0]= TDS_helper::origin(e)->point();
      pts[1]= TDS_helper::third_vertex(e)->point();
      pts[2]= TDS_helper::destination(e)->point();
      pts[3]= TDS_helper::mirror_vertex(e)->point();
      vhs[0]= TDS_helper::origin(e);
      vhs[1]= TDS_helper::third_vertex(e);
      vhs[2]= TDS_helper::destination(e);
      vhs[3]= TDS_helper::mirror_vertex(e);
      CGAL_UD_DEBUG("Trying " << pts[0] << " " << pts[1] << " " << pts[2] << " " << pts[3] << std::endl);

      if ( active.find(pts[0]) != active.end()
	   && active.find(pts[1]) != active.end()
	   && active.find(pts[2]) != active.end()
	   && active.find(pts[3]) != active.end()) {
	return false;
      }
      if (pts[0] == Point_key() || pts[1]==Point_key() || pts[2] == Point_key() || pts[3] == Point_key()) {
	return false;
      }
      bool ok;
      if (PASS == 0 || PASS == 2) {
	Protected_array<typename Indirect_kernel::Geometric_point_2, 4> pcs;
	for (unsigned int i = 0; i< 4; ++i) {
	  if (PASS != 0 && active.find(pts[i]) == active.end()) {
	    pcs[i]= ui()->initial_coordinates_object()(pts[i]);
	  } else {
	    pcs[i]= ui()->final_coordinates_object()(pts[i]);
	  }
	}
	++stat_point_predicate_evaluations_;

	ok= (ui()->initial_kernel_object().direct_kernel_object().side_of_oriented_circle_2_object()(pcs[0], pcs[1], pcs[2], pcs[3])
	     == CGAL::ON_POSITIVE_SIDE);
      } else {
	CGAL::Protect_FPU_rounding<true> prot;
	++stat_interval_predicate_evaluations_;

	ok= ( ui()->sign_at(pts[0], pts[1], pts[2], pts[3], active) == CGAL::POSITIVE);
      }
      if (!ok) {
	for (unsigned int i=0; i< 4; ++i) {
	  if (active.find(pts[i]) == active.end()) {
	    active.insert(pts[i]);
	    new_active.push_back(vhs[i]);
	  }
	}
	return true;
      } else {
	return false;
      }
    }
#else
    // fix update
    bool test_and_add(Edge e) {
      ++stat_point_predicate_evaluations_;
      if (!compute_ok(e, ui()->final_kernel_object())) {
	ui()->schedule_activation(e.first->vertex(0)->point(), tr_.simulator_handle());
	ui()->schedule_activation(e.first->vertex(1)->point(), tr_.simulator_handle());
	ui()->schedule_activation(e.first->vertex(2)->point(), tr_.simulator_handle());
	ui()->schedule_activation(TDS_helper::mirror_vertex(e)->point(), tr_.simulator_handle());
	return true;
      } else return false;
    }

    /*bool test_and_add_one(Edge e, std::vector<Point_key> &active) const {
      ++stat_point_predicate_evaluations_;
      if (!compute_ok(e, ui()->final_kernel_object())) {
	//add(e.first->vertex(0)->point(), active);
	//add(e.first->vertex(1)->point(), active);
	//add(e.first->vertex(2)->point(), active);
	add(TDS_helper::mirror_vertex(e)->point(), active);
	return true;
      } else return false;
      }*/
#endif
#endif

    template <class TDS>
    void initialize_events(const TDS &triangulation,
			   Indirect_kernel fk) {
      if (stat_number_of_edges_==0) {
	stat_number_of_edges_=triangulation.tds().number_of_edges();
	stat_number_of_vertices_= triangulation.tds().number_of_vertices();
      }
      ui()->set_final_kernel(fk);
#ifndef MOVE_ALL
      /*#ifdef HYBRID
      // first check if all points are moving, if so skip
      // anything can only be done twice
      // take an edge, if it fails, add all 4 points and other two edges on other face
      std::set<Point_key> active;
      std::vector<typename TDS::Vertex_handle> new_active;
      std::vector<typename TDS::Edge> edges;
      for (typename TDS::Finite_edges_iterator it= triangulation.finite_edges_begin();
	   it != triangulation.finite_edges_end(); ++it){
	if (test_and_add<0>(*it, active, new_active)) {
	  ++stat_number_of_bad_edges_;
	}
      }

      while (!new_active.empty()) {
	std::vector<typename TDS::Vertex_handle> nn_active;
	for (unsigned int i=0; i < new_active.size(); ++i) {
	  typename TDS::Face_circulator c= triangulation.incident_faces(new_active[i]), n=c;
	  ++n;
	  do {
	    typename TDS::Face_handle e= *c;
	    for (unsigned int i=0; i< 3; ++i) {
	      if (e->neighbor(i) != n) {
		test_and_add<2>(Edge(e, i), active, nn_active);
	      }
	    }
	    c=n;
	    ++n;
	  } while (c != triangulation.incident_edges(new_active[i]));
	}
	std::swap(nn_active, new_active);
      }

      CGAL_UD_DEBUG( stat_point_predicate_evaluations_ << " and "
		     << stat_interval_predicate_evaluations_ << " evals on initialization" << std::endl);

		     #else*/
      //std::vector<Point_key> active;
      for (typename TDS::Finite_edges_iterator it= triangulation.finite_edges_begin();
	   it != triangulation.finite_edges_end(); ++it){
	if (test_and_add(*it)) {
	  ++stat_number_of_bad_edges_;
	}
      }
      /*std::sort(active.begin(), active.end());
	active.erase(std::unique(active.begin(), active.end()), active.end());*/
      //#endif
#endif

    }



    void after_flip(Edge ) {

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
    double activation_time_;
    Event_key activation_event_;
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
    kdel_->set_has_certificates(true, KDel::HAS_NO_FAILURES);
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

    kdel_->set_has_certificates(false, KDel::HAD_NO_FAILURES);
    tr_.simulator_handle()->set_interval(0,1);

    Indirect_kernel fk;
    fk.new_point_2s(pts.begin(), pts.end());
    kdel_->visitor().initialize_events(kdel_->triangulation(), fk);
#ifdef MOVE_ALL
    kdel_->set_has_certificates(true, KDel::HAD_NO_FAILURES | KDel::NO_STRUCTURE_CHANGES);
#else
    kdel_->set_has_certificates(true, KDel::HAD_NO_FAILURES | KDel::HAS_NO_FAILURES | KDel::NO_STRUCTURE_CHANGES);
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
      CGAL_error_msg( std::string("Error opening file ") + name );
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


    CGAL_SET_LOG_LEVEL(CGAL::Kinetic::Log::SOME);



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

} //namespace CGAL

#endif
#undef CGAL_UD_DEBUG
