#ifndef SORT_TEST_H
#define SORT_TEST_H
#include <CGAL/KDS/Sort.h>

template <class Traits>
bool sort_test(Traits &tr, double max_events=std::numeric_limits<double>::infinity())
{
  std::string etag="WARNING: ";
  bool eret=false;
  CGAL_exactness_assertion_code(etag="ERROR: ");
  //CGAL_exactness_assertion_code(bool test_compiled_with_exact_checks;);
  CGAL_exactness_assertion_code(eret=true);
  
  typedef CGAL::KDS::Sort<Traits> Sort;
  Sort sort(tr);

  //#ifdef _MSC_VER
  ///#pragma warning(disable:1572)
  //#endif

  while (tr.simulator_pointer()->next_event_time() != tr.simulator_pointer()->end_time()) {
    //std::cout << *tr.simulator_pointer() << std::endl;
    tr.simulator_pointer()->set_current_event_number(tr.simulator_pointer()->current_event_number()+1);
    if (tr.simulator_pointer()->current_event_number() > max_events){
      std::cerr << "ERROR too many events" << std::endl;
      std::cerr << *tr.active_objects_table_pointer() << std::endl;
    }
  }

  //#ifdef _MSC_VER
  ///#pragma warning(enable:1572)
  //#endif

  std::cout << tr.simulator_pointer()->current_event_number() << " events processed.\n";

  typedef typename Sort::Key_iterator Kit;
  Kit c= sort.begin();
  Kit b=c;
  ++c;
  bool error =false;

  typename Traits::NT ratt;
  if (tr.simulator_pointer()->has_audit_time()) {
    ratt=tr.simulator_pointer()->audit_time();
  }
  else {
    std::cerr << etag << "Out of events, but the time is not rational." << std::endl;
    std::cerr << "Current time is " << tr.simulator_pointer()->current_time()
	      << " the end time is " << tr.simulator_pointer()->end_time() << std::endl;
    ratt= CGAL::to_interval(tr.simulator_pointer()->end_time()).second;
    error=eret;
  }

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
  return error;
}
#endif
