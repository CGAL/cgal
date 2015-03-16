#ifndef SORT_TEST_H
#define SORT_TEST_H
#include <CGAL/Timer.h>

#include <CGAL/Kinetic/Sort.h>

template <class IComp>
struct ILess {
  ILess(IComp ic): ic_(ic){}
  template <class A>
  bool operator()(const A &a, const A &b) const {
    return ic_(a,b) == CGAL::SMALLER;
  }
  IComp ic_;
};

template <class Traits>
bool sort_test(Traits &tr, double max_events=std::numeric_limits<double>::infinity())
{

  std::string etag=":NUMERICAL ISSUE: ";
  //CGAL_exactness_assertion_code(bool fail=false);
  CGAL_exactness_assertion_code(etag="ERROR: ");
  //CGAL_exactness_assertion_code(fail=true);
  //CGAL_exactness_assertion_code(bool test_compiled_with_exact_checks;);
    
  typedef CGAL::Kinetic::Sort<Traits> Sort;
  Sort sort(tr);

  //#ifdef _MSC_VER
  ///#pragma warning(disable:1572)
  //#endif

  CGAL::Timer t;
  t.start();

  while (tr.simulator_handle()->next_event_time() != tr.simulator_handle()->end_time()) {
    //std::cout << *tr.simulator_pointer() << std::endl;
    tr.simulator_handle()->set_current_event_number(tr.simulator_handle()->current_event_number()+1);
#ifndef NDEBUG
    if (tr.simulator_handle()->current_event_number() > max_events){
      std::cerr << "ERROR too many events" << std::endl;
      ++CGAL::Kinetic::internal::get_static_audit_failures();
      std::cerr << *tr.active_points_1_table_handle() << std::endl;
    }
#endif
  }
  t.stop();
  double tpe= t.time()/tr.simulator_handle()->current_event_number();
  std::cout  << "Time per event is " << tpe*1000 << std::endl;

  //#ifdef _MSC_VER
  ///#pragma warning(enable:1572)
  //#endif

  std::cout << tr.simulator_handle()->current_event_number() << " events processed.\n";


  bool error =false;
#ifndef NDEBUG
  bool eret=false;
  CGAL_exactness_assertion_code(eret=true);
 
  typedef typename Sort::Iterator Kit;
  Kit c= sort.begin();
  Kit b=c;
  ++c;
  typename Traits::Simulator::NT ratt;
  //if (tr.simulator_handle()->next_time_representable_as_nt()) {
  ratt=tr.simulator_handle()->next_time_representable_as_nt();
  /*} else {
    std::cerr << etag << "Out of events, but the time is not rational." << std::endl;
    std::cerr << "Current time is " << tr.simulator_handle()->current_time()
    << " the end time is " << tr.simulator_handle()->end_time() << std::endl;
    CGAL::Kinetic::internal::fail__|= fail;
    ratt= CGAL::to_interval(tr.simulator_handle()->end_time()).second;
    error=eret;
    }*/
  
  std::cout << "End time is " << tr.simulator_handle()->end_time() << std::endl;

  while (c != sort.end()) {
    typename Traits::Simulator::Function_kernel::Function f= tr.active_points_1_table_handle()->at(*c).x() - tr.active_points_1_table_handle()->at(*b).x();
    if ( f(ratt) < 0 ) {
      std::cerr << etag << "Objects " << c->object() << " = " << tr.active_points_1_table_handle()->at(*c).x() << " and "
                << b->object() << " = " << tr.active_points_1_table_handle()->at(*b).x() << " out of order at end of time "
                <<ratt << std::endl;
      ++CGAL::Kinetic::internal::get_static_audit_failures();
      error=eret;
    }

    ++b;
    ++c;
  }
#endif
  return error;
}
#endif
