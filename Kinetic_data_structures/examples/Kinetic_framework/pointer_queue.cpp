#include <CGAL/Kinetic/basic.h>
#include <CGAL/Kinetic/Heap_pointer_event_queue.h>
#include <CGAL/Kinetic/Inexact_simulation_traits.h>
#include <CGAL/Kinetic/Event_base.h>
#include <cstdlib>

typedef double Time;

Time proc_time_=-1;


class Event: public CGAL::Kinetic::Event_base<void*>
{
public:
  Event(double i): i_(i) {
  }
  Time time() const
  {
    return Time(i_);
  };
  void process(){
    //std::cout << "Event at " << i_ << "\n";
    proc_time_=i_;
    /*if (sim_->current_time() != i_) {
      std::cerr << "ERROR: Times do not match. Got " << sim_->current_time()
		<< " expected " << i_ <<std::endl;
    }
    assert(sim_->current_time()==i_);*/
  }
protected:
  Time i_;
};

std::ostream &operator<<(std::ostream &out, Event e)
{
  out << e.time();
  return out;
}


int main(int, char *[])
{
  typedef CGAL::Kinetic::Inexact_simulation_traits::Kinetic_kernel::Function_kernel FK;
  typedef CGAL::Kinetic::Heap_pointer_event_queue<FK> Q;
  Q pq(0, 10000, FK());
  typedef Q::Key Key;
  std::vector<Key>  items;

  for (unsigned int i=0; i< 10000; ++i) {
    Time t(std::rand()/100000.0);
    items.push_back(pq.insert(t,Event(t)));
  }
  for (unsigned int i=0; i< 5000; ++i) {
    pq.erase(items[i]);
  }

  Time last_time = -1;
  while (!pq.empty()) {
    Time t= pq.front_priority();
    if (t < last_time) {
      std::cerr << "ERROR: priority of next event (" << pq.front_priority()
		<< ") is before the last one (" << last_time << ").\n";
    }
    assert(t >= last_time);
    last_time=t;
    pq.process_front();
    if (t < last_time) {
      std::cerr << "ERROR: wrong event processed (" << pq.front_priority()
		<< ") instead of (" << proc_time_ << ").\n";
    }
    if (proc_time_!= t) {
      std::cerr << "ERROR: wrong event time. Expected: " << t << " got " << proc_time_ <<std::endl;
    }
    assert(proc_time_==t);
  }
  return EXIT_SUCCESS;
}
