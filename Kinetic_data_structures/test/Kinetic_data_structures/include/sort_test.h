#ifndef SORT_TEST_H
#define SORT_TEST_H
#include <CGAL/Timer.h>

#include <CGAL/Kinetic/Sort.h>
template <class Traits>
bool sort_test(Traits &tr, double max_events=std::numeric_limits<double>::infinity())
{
  std::string etag="WARNING: ";

  CGAL_exactness_assertion_code(etag="ERROR: ");
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
      std::cerr << *tr.active_objects_table_handle() << std::endl;
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
 
  typedef typename Sort::Key_iterator Kit;
  Kit c= sort.begin();
  Kit b=c;
  ++c;
  typename Traits::NT ratt;
  if (tr.simulator_pointer()->has_rational_current_time()) {
    ratt=tr.simulator_pointer()->rational_current_time();
  } else {
    std::cerr << etag << "Out of events, but the time is not rational." << std::endl;
    std::cerr << "Current time is " << tr.simulator_pointer()->current_time()
	      << " the end time is " << tr.simulator_pointer()->end_time() << std::endl;
    ratt= CGAL::to_interval(tr.simulator_pointer()->end_time()).second;
    error=eret;
  }

  std::cout << "End time is " << tr.simulator_pointer()->end_time() << std::endl;

  while (c != sort.end()) {
    typename Traits::Simulator::Function_kernel::Function f= tr.active_objects_table_pointer()->at(*c).x() - tr.active_objects_table_pointer()->at(*b).x();
    if ( f(ratt) < 0 ) {
      std::cerr << etag << "Objects " << *c << " = " << tr.active_objects_table_pointer()->at(*c).x() << " and "
                << *b << " = " << tr.active_objects_table_pointer()->at(*b).x() << " out of order at end of time "
                <<ratt << std::endl;
      error=eret;
    }

    ++b;
    ++c;
  }
#endif
  return error;
}
#endif
