#include <CGAL/KDS/Ref_counted.h>
#include <CGAL/KDS/Exact_simulation_traits_1.h>
#include <CGAL/KDS/Notifying_table_listener_helper.h>
#include <CGAL/KDS/Simulator_kds_listener.h>
#include <iostream>
#include <set>
#include <vector>


// This must be external since operator<< has to be defined
template <class Object, class Time, class KDS>
struct Trivial_event{
  Trivial_event(){}
  template <class It, class Map>
  Trivial_event(It beg, It end, const Map &m, KDS* kds): kds_(kds) {
    for (; beg != end; ++beg){
      objects_.push_back(m[*beg]);
    }
  }
  void process(Time t) const {
    std::cout << "At time " << t 
	      << " the following objects are in the table: ";
    for (typename std::vector<Object>::const_iterator
	   cit= objects_.begin(); 
	 cit != objects_.end(); ++cit){
      std::cout << *cit << " ";
    }
    std::cout << std::endl;
    kds_->set_processed(true);
  }
    
  std::vector<Object> objects_;
  typename KDS::Pointer kds_;
};

template <class Object, class Time, class KDS>
std::ostream &operator<<(std::ostream &out,
			 const Trivial_event<Object, Time, KDS> &) {
  out << "\"An event\"";
  return out;
}

/*!  This is a trivial kinetic data structure which doesn't actually
  do anything, other than print what it should be doing and maintain 1
  certificate in the event queue which has a list of all the moving
  objects in the moving object table.
*/
template <class Traits>
struct Trivial_kds: CGAL::KDS::Ref_counted<Trivial_kds<Traits> > {
  typedef Trivial_kds<Traits> This;
  typedef typename Traits::Moving_point_table::Data Point;
  typedef typename Traits::Simulator::Time Time;
  typedef typename Traits::Moving_point_table::Key Point_key;
  typedef typename Traits::Simulator::Event_key Event_key;
  typedef CGAL::KDS::Notifying_table_listener_helper<
    typename Traits::Moving_point_table::Listener, This> Notifying_table_helper;
  typedef CGAL::KDS::Simulator_kds_listener<
    typename Traits::Simulator::Listener, This> Simulator_helper;

  typedef Trivial_event<Point, Time, This> Event;
 

  Trivial_kds(Traits tr): has_certificates_(true),
			  tr_(tr),
			  nth_(tr.moving_point_table_pointer(), this),
			  sh_(tr.simulator_pointer(), this){}

  // this method is called with the value true when the event is processed
  void set_processed(bool tf) {
    if (tf== true) {
      event_= Event_key();
      set_has_certificates(false);
      set_has_certificates(true);
    }
  }

  void audit() const {
    /* In a real KDS you could use tr_.instantaneous_kernel() to build
       a static version of the kinetic data structure and check it for
       equality with the kinetic version.
    */
    assert(event_);
    std::cout << "The structure is trivially correct at time " 
	      << tr_.simulator_pointer()->audit_time() << std::endl;
  }

  void set_has_certificates(bool tf) {
    typename Traits::Simulator::Pointer sp= tr_.simulator_pointer();
    if (has_certificates_ != tf){
      has_certificates_=tf;
      if (has_certificates_){
	assert(!event_);
	std::cout << "Time to rebuild all certificates from scratch." << std::endl;
	Time t= CGAL::to_interval(sp->current_time()).second+1;
	event_= sp->new_event(t, Event(objects_.begin(), 
				       objects_.end(),
				       *tr_.moving_point_table_pointer(), 
				       this));
	std::cout << "Created event (" << event_ << ") at time " << t << std::endl;
	assert(!event_ 
	       || sp->event(event_, Event()).objects_.size() == objects_.size());
      } else {
	if (event_) {
	  std::cout << "Time to delete all certificates." << std::endl;
	  std::cout << "Deleting event " << event_ << std::endl;
	  sp->delete_event(event_);
	  event_=Event_key();
	}
      }
    } 
  }

  bool has_certificates() const {
    // you can't use !event_ because the simulation
    // might end before the time you tried to schedule the event.
    return has_certificates_;
  }

  
  void insert(Point_key k) {
    std::cout << "Updating structure to include new object " 
	      << k << "." << std::endl;
    objects_.insert(k);
    if (has_certificates_){
      std::cout << "Updating all certificates which depend on "
		<< k << "." << std::endl;
      set_has_certificates(false);
      set_has_certificates(true);
    }
  }


  void set(Point_key k) {
    std::cout << "An object changed, but we shouldn't every update" 
	      << " structure when an object changes."<< std::endl;
    if (has_certificates_){
      std::cout << "Updating all certificates which depend on " 
		<< k << "." << std::endl;
      set_has_certificates(false);
      set_has_certificates(true);
    }
  }

  void erase(Point_key k) {
    std::cout << "An object " << k << " was removed."<< std::endl;
    objects_.erase(k);
    if (has_certificates_){
      std::cout << "Updating all certificates which depend on " 
		<< k << "." << std::endl;
      set_has_certificates(false);
      set_has_certificates(true);
    }
  }



protected:
  bool has_certificates_;
  std::set<Point_key > objects_;
  Event_key event_;
  Traits tr_;
  Notifying_table_helper nth_;
  Simulator_helper sh_;
};

int main(int, char *[]){
  typedef CGAL::KDS::Exact_simulation_traits_1 Traits;
  typedef Trivial_kds<Traits> TKDS;

  Traits tr;
  TKDS::Pointer tk= new TKDS(tr);

  Traits::Simulator::Pointer sp=tr.simulator_pointer();

  Traits::Simulator::Function_kernel::Construct_function cf
    = sp->function_kernel_object().construct_function_object();

  tk->set_has_certificates(true);

  typedef Traits::Kinetic_kernel::Point_2 Point;
  Traits::Moving_point_table::Key k
    = tr.moving_point_table_pointer()->insert(Point(cf(1), cf(0)));
  sp->set_current_event_number(sp->current_event_number()+10);
  tr.moving_point_table_pointer()
    ->set(k, Traits::Kinetic_kernel::Point_2(cf(2), cf(0)));
  sp->set_current_event_number(sp->current_event_number()+10);
  tr.moving_point_table_pointer()->erase(k); 
  sp->set_current_event_number(sp->current_event_number()+10);
  return EXIT_SUCCESS;
}
