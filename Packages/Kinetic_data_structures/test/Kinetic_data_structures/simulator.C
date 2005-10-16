#include <CGAL/KDS/Exact_simulation_traits_2.h>
#include <CGAL/KDS/Inexact_simulation_traits_2.h>
#include <vector>
#include <algorithm>
#include <map>
#include <CGAL/KDS/Two_list_pointer_event_queue.h>
#include <CGAL/Timer.h>

const int num_events=1000;

template <class T> 
struct Test_event {
  typename T::Time t_;
  T *tt_;

  Test_event(typename T::Time t, T *tt): t_(t), tt_(tt){}
  const typename T::Time &time() const {
    return t_;
  }

  Test_event(){}

  void process(const typename T::Time& t) {
    if (tt_->expected_next_event_time() != t) {
      std::cerr << "ERROR Wrong event occurred. Expecting time " << tt_->expected_next_event_time() 
		<< " got time " << t << std::endl;
    }
    if (t != t_){
      std::cerr << "ERROR Event occured at wrong time. Expecting time " << t_ << " got time " << t << std::endl;
    }
    tt_->process_one_event();
  }
  bool operator==(const Test_event<T> &o) const {
    return t_== o.t_&& tt_==o.tt_;
  }
};

template <class T>
std::ostream& operator<<(std::ostream &out, const Test_event<T> &e) {
  return out << e.t_;
}



template <class Sim>
struct Test {
  typedef Test<Sim> This;
  typedef typename Sim::Time Time;
  typedef typename Sim::Event_key Key;
  typedef typename Sim::Function_kernel FK;
  
  typedef Test_event<Test<Sim> > Ev;
  Ev ev;
  typename Sim::Pointer sim;
  FK fk;

  std::map<Key, Time> event_map;
  std::map<Time, Key> time_map;
  std::vector<Key> to_delete;

  Time expected_next_event_time() const{
    return time_map.begin()->first;
  }

  void process_one_event() {
    event_map.erase(time_map.begin()->second);
    time_map.erase(time_map.begin());
  }

  typename FK::NT  coef(){
    if (std::rand()%2==0) 
      return typename FK::NT(static_cast<double>(std::rand())/static_cast<double>(RAND_MAX));
    else 
      return -typename FK::NT(static_cast<double>(std::rand())/static_cast<double>(RAND_MAX));
  }

  void add_events(unsigned int num){
    std::cout << "Adding events..."<< std::flush;
    typename FK::Construct_function cf= fk.construct_function_object();

    for (unsigned int i=0; i<num; ++i){
      typename FK::Function f= cf(abs(coef()), coef(), coef(), coef(), coef());
      if (sim->function_kernel_object().sign_at_object(f)(sim->current_time()) == CGAL::NEGATIVE){
	f=-f;
	assert(sim->function_kernel_object().sign_at_object(f)(sim->current_time()) != CGAL::NEGATIVE);
      }

      typename Sim::Root_stack s= sim->root_stack_object(f);
      
      while ( !s.empty() ){
	
	if (time_map.find(s.top()) == time_map.end()){
	  Ev ev(s.top(), this);
	  Key k=sim->new_event(s.top(), ev);
	  assert(sim->template event<Ev>(k) == ev);
	  if (s.top() >= sim->end_time()){
	    assert(k==sim->null_event());
	  }
	  if (k==sim->null_event()){
	    assert(s.top() >= sim->end_time());
	  } else {
	    event_map[k]= s.top();
	    time_map[s.top()]=k;
	    
	    if (std::rand()%3==0){
	      to_delete.push_back(k);
	    } 
	  }
	}
	s.pop();
      } 
    }
    std::cout << "done" << std::endl;
  }

  void check_event(Key k, Time t){
    const Ev &mev = sim->template event<Ev>(k);
    if (mev.time() != t){
      std::cerr << "ERROR Messed up event in queue: expected " << t << " got " << mev.time() 
		<< std::endl;
    }
  }

  void delete_events() {
    std::cout << "Deleting events..." << std::flush;
    while (!to_delete.empty()) {
      Key k= to_delete.back();
      check_event(to_delete.back(), event_map[k]);
      //std::cout << sim << std::endl;
      time_map.erase(event_map[k]);
      event_map.erase(k);
      sim->delete_event(k);
      to_delete.pop_back();
    }
    std::cout << "done" << std::endl;
  }


  void check_all_events() {
    std::cout << "Checking events..." << std::flush;
    for (typename std::map<Key, Time>::const_iterator it= event_map.begin();
	 it != event_map.end(); ++it){
      check_event(it->first, it->second);
    }
    std::cout << "done" << std::endl;
  }
  
  Test(){
    sim= new Sim(0,1000);
    //}

    //void run () {


    add_events(num_events);

    check_all_events();

    delete_events();

    

    //std::cout << sim << std::endl;
  
    

    for (int i=0; i< num_events/4; ++i){
      sim->set_current_event_number(sim->current_event_number()+1);
    }

    check_all_events();

   
    add_events(num_events);
    delete_events();

    sim->set_current_time(std::numeric_limits<typename Sim::Time>::infinity());

    std::cout << "Done checking simulator.\n" << std::endl;

    assert(time_map.empty());
    assert(event_map.empty());
  }
};

/*
  Time on primal for 10000 events (2x) with my heap was 140s.
  
*/

int main(int, char *[]){
  CGAL::Timer timer;

  {
    CGAL::KDS::Exact_simulation_traits_2::Simulator ts;
    ts.event<Test_event<CGAL::KDS::Exact_simulation_traits_2::Simulator> >(CGAL::KDS::Exact_simulation_traits_2::Simulator::Event_key());
  }
  {
    timer.start();
    Test<CGAL::KDS::Exact_simulation_traits_2::Simulator> ts;
    //ts.run();
    timer.stop(); 
    //std::cout << "Bin heap time " << timer.time() << std::endl;
  }
  {
    timer.reset();

    typedef CGAL::KDS::Exact_simulation_traits_2::Function_kernel FK;
    typedef CGAL::KDS::Simulator<CGAL::KDS::Exact_simulation_traits_2::Simulator::Function_kernel,
      CGAL::KDS::Two_list_pointer_event_queue<FK::Root, FK::NT> > Sim2;

    timer.start();
    Test<Sim2> t2;
    timer.stop();
    //std::cout << "List time " << timer.time() << std::endl;
  
    //assert(two_list_remaining==0);
  }

  // if (!error) return EXIT_SUCCESS;
  //else 
  return EXIT_SUCCESS;
}
