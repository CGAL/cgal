#define CGAL_CHECK_EXACTNESS
#define CGAL_CHECK_EXPENSIVE

#include <CGAL/Kinetic/Exact_simulation_traits.h>
#include <cassert>
#include <CGAL/Kinetic/Inexact_simulation_traits.h>
#include <CGAL/Kinetic/Two_list_pointer_event_queue.h>
#include <CGAL/Kinetic/Event_base.h>
#include <CGAL/Timer.h>
#include <algorithm>
#include <map>
#include <vector>

const int num_events=200;

template <class T, class Time>
struct Test_event: public CGAL::Kinetic::Event_base<void*>
{
  Time t_;
  T *tt_;

  Test_event(Time t, T *tt): t_(t), tt_(tt){}
  const Time &time() const
  {
    return t_;
  }

  Test_event(){}

  void process() {
    if (tt_->expected_next_event_time() != tt_->current_time()) {
      std::cerr << "ERROR Wrong event occurred. Expecting time " << tt_->expected_next_event_time()
                << " got time " << tt_->current_time() << std::endl;
    }
    if (tt_->current_time() != t_) {
      std::cerr << "ERROR Event occured at wrong time. Expecting time " << t_ << " got time " << tt_->current_time() << std::endl;
    }
    tt_->process_one_event();
  }
  bool operator==(const Test_event<T, Time> &o) const
  {
    return t_== o.t_&& tt_==o.tt_;
  }
  std::ostream &write(std::ostream &out) const {
    return out << t_;
  }
};

template <class Sim>
struct Test
{
  typedef Test<Sim> This;
  typedef typename Sim::Time Time;
  typedef typename Sim::Event_key Key;
  typedef typename Sim::Function_kernel FK;

  typedef Test_event<This, Time> Ev;
  Ev ev;
  typename Sim::Handle sim;
  FK fk;

  std::map<Key, Time> event_map;
  std::map<Time, Key> time_map;
  std::vector<Key> to_delete;
  std::set<Key> inf_events_;

  Time expected_next_event_time() const
  {
    if (!time_map.empty()) {
      return time_map.begin()->first;
    } else {
      return sim->end_time();
    }
  }


  Time current_time() const
  {
    return sim->current_time();
  }

  void process_one_event() {
    if (!time_map.empty()) {
      event_map.erase(time_map.begin()->second);
      time_map.erase(time_map.begin());
    } else {
      assert(!inf_events_.empty());
      // just look at count for now
      inf_events_.erase(inf_events_.begin());
    }
  }

  typename FK::Function::NT  coef() {
    if (std::rand()%2==0)
      return typename FK::Function::NT(static_cast<double>(std::rand())/static_cast<double>(RAND_MAX));
    else
      return -typename FK::Function::NT(static_cast<double>(std::rand())/static_cast<double>(RAND_MAX));
  }

  void add_events(unsigned int num) {
    std::cout << "Adding events..."<< std::flush;
    typename FK::Construct_function cf= fk.construct_function_object();
    FK fk;
    for (unsigned int i=0; i<num-3; ++i) {
      typename FK::Function f= cf(CGAL::abs(coef()), coef(), coef(), coef(), coef());
      if (fk.sign_at_object()(f, sim->current_time()) == CGAL::NEGATIVE) {
	f=-f;
	assert(fk.sign_at_object()(f, sim->current_time()) != CGAL::NEGATIVE);
      }

      typename FK::Root_stack s= fk.root_stack_object(f, 
						      sim->current_time(),
						      sim->end_time());

      while ( !s.empty() ) {

	if (time_map.find(s.top()) == time_map.end()) {
	  Ev ev(s.top(), this);
	  Key k=sim->new_event(s.top(), ev);
	  assert(sim->template event<Ev>(k) == ev);
	  if (s.top() >= sim->end_time()) {
	    assert(k==sim->null_event());
	  }
	  if (k==sim->null_event()) {
	    assert(s.top() > sim->end_time());
	  }
	  else {
	    event_map[k]= s.top();
	    time_map[s.top()]=k;

	    if (std::rand()%3==0) {
	      to_delete.push_back(k);
	    }
	  }
	}
	s.pop();
      }
    }
    for (unsigned int i=0; i< 3; ++i) {
      Key k= sim->new_final_event( Ev(sim->end_time(), this));
      assert(k != sim->null_event());
      inf_events_.insert(k);
      
      if (std::rand()%3==0) {
	to_delete.push_back(k);
      }
      
    }
    std::cout << "done" << std::endl;
  }

  void check_event(Key k, Time t) {
    const Ev &mev = sim->template event<Ev>(k);
    if (mev.time() != t) {
      std::cerr << "ERROR Messed up event in queue: expected " << t << " got " << mev.time()
                << std::endl;
    }
  }

  void delete_events() {
    std::cout << "Deleting events..." << std::flush;
    while (!to_delete.empty()) {
      Key k= to_delete.back();
       //std::cout << sim << std::endl;
      if (event_map.find(k) != event_map.end()) {
	check_event(to_delete.back(), event_map[k]);
	time_map.erase(event_map[k]);
	event_map.erase(k);
	sim->delete_event(k);
      } else {
	inf_events_.erase(k);
	sim->delete_event(k);
      }
      to_delete.pop_back();
    }
    std::cout << "done" << std::endl;
  }

  void check_all_events() {
    std::cout << "Checking events..." << std::flush;
    for (typename std::map<Key, Time>::const_iterator it= event_map.begin();
	 it != event_map.end(); ++it) {
      check_event(it->first, it->second);
    }
    std::cout << "done" << std::endl;
  }

  Test() {
    test_to(Time(1000));
    test_to(std::numeric_limits<Time>::infinity());
  }

  void test_to(Time end) {
    sim= new Sim(0,end);
    //}

    //void run () {

    add_events(num_events);

    check_all_events();

    delete_events();

    //std::cout << sim << std::endl;

    for (int i=0; i< num_events/4; ++i) {
      sim->set_current_event_number(sim->current_event_number()+1);
    }

    check_all_events();

    add_events(num_events);
    delete_events();

    sim->set_current_time(sim->end_time());
    assert(inf_events_.size() != 0);
    while (!sim->empty()) {
      sim->set_current_event_number(sim->current_event_number()+1);
    }

    std::cout << "Done checking simulator.\n" << std::endl;

    assert(time_map.size()==0);
    assert(event_map.size()==0);
    assert(inf_events_.size()==0);
  }
};

/*
  Time on primal for 10000 events (2x) with my heap was 140s.

*/

int main(int, char *[])
{
  CGAL::Timer timer;
  CGAL_SET_LOG_LEVEL(CGAL::Log::NONE);

 {
    timer.reset();

    typedef CGAL::Kinetic::Inexact_simulation_traits::Kinetic_kernel::Function_kernel FK;
    typedef CGAL::Kinetic::Default_simulator<FK,
     CGAL::Kinetic::Two_list_pointer_event_queue<FK> > Sim2;

    timer.start();
    Test<Sim2> t2;
    timer.stop();
    //std::cout << "List time " << timer.time() << std::endl;

    //assert(two_list_remaining==0);
  }

  {
    timer.start();
    Test<CGAL::Kinetic::Exact_simulation_traits::Simulator> ts;
    //ts.run();
    timer.stop();
    //std::cout << "Bin heap time " << timer.time() << std::endl;
  }

  {
    timer.start();
    Test<CGAL::Kinetic::Inexact_simulation_traits::Simulator> ts;
    //ts.run();
    timer.stop();
    //std::cout << "Bin heap time " << timer.time() << std::endl;
  }
 
  // if (!error) return EXIT_SUCCESS;
  //else
  if (CGAL::Kinetic::internal::get_static_audit_failures() != 0 ) return EXIT_FAILURE;
  else return EXIT_SUCCESS;
}
