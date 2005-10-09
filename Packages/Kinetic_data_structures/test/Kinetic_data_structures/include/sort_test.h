#ifndef SORT_TEST_H
#define SORT_TEST_H
#include <CGAL/KDS/Sort.h>

template <class Traits>
bool sort_test(Traits &tr){
  typedef CGAL::KDS::Sort<Traits> Sort;
  Sort sort(tr);
  
  while (tr.simulator_pointer()->next_event_time() != tr.simulator_pointer()->end_time()){
    //std::cout << *tr.simulator_pointer() << std::endl;
    tr.simulator_pointer()->set_current_event_number(tr.simulator_pointer()->current_event_number()+1);
  }

  std::cout << tr.simulator_pointer()->current_event_number() << " events processed.\n";

  typedef typename Sort::Key_iterator Kit;
  Kit c= sort.begin();
  Kit b=c;
  ++c;
  bool error =false;

  typename Traits::NT ratt;
  if (tr.simulator_pointer()->has_rational_current_time()){
    ratt=tr.simulator_pointer()->rational_current_time();
  } else {
    std::cerr << "Out of events, but the time is not rational." << std::endl;
    std::cerr << "Current time is " << tr.simulator_pointer()->current_time() 
	      << " the end time is " << tr.simulator_pointer()->end_time() << std::endl;
    ratt= CGAL::to_interval(tr.simulator_pointer()->end_time()).second;
  }

  while (c != sort.end()){
    if (tr.moving_point_table_pointer()->at(*c).x()(ratt) 
	< tr.moving_point_table_pointer()->at(*b).x()(ratt)){
      std::cerr << "Objects " << *c << " = " << tr.moving_point_table_pointer()->at(*c).x() << " and " 
		<< *b << " = " << tr.moving_point_table_pointer()->at(*b).x() << " out of order at end of time " 
		<< tr.simulator_pointer()->rational_current_time() << std::endl;
      error=true;
    }
    
    ++b;
    ++c;
  }
  return error;
}

#endif
