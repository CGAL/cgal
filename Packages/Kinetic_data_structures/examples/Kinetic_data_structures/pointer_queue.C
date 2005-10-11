#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/Heap_pointer_event_queue.h>
#include <cstdlib>

typedef double Time;

Time proc_time_=-1;

class Event {
public:
  Event(double i): i_(i){
  }
  Time time() const {
    return Time(i_);
  };
  void process(double t) const {
    std::cout << i_ << "\n";
    proc_time_=i_;
    assert(t==i_);
  }
protected:
  Time i_;
};

std::ostream &operator<<(std::ostream &out, Event e){
  out << e.time();
  return out;
}

int main(int, char *[]){
  typedef CGAL::KDS::Heap_pointer_event_queue<Time> Q;
  Q pq(100);
  typedef Q::Key Key;
  std::vector<Key>  items;

  for (unsigned int i=0; i< 10000; ++i){
    Time t(std::rand());
    items.push_back(pq.insert(t,Event(t)));
  }
  for (unsigned int i=0; i< 5000; ++i){
    pq.erase(items[i]);
  }

  Time last_time = -1;
  while (!pq.empty()){
    Time t= pq.front_priority();
    assert(t >= last_time);
    last_time=t;
    pq.process_front();
    assert(proc_time_==t);
  }
  return EXIT_SUCCESS;
}
