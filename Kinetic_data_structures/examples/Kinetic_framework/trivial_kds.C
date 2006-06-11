#include <CGAL/Kinetic/Ref_counted.h>
#include <CGAL/Kinetic/Exact_simulation_traits_1.h>
#include <CGAL/Kinetic/Active_objects_listener_helper.h>
#include <CGAL/Kinetic/Simulator_kds_listener.h>
#include <iostream>
#include <set>
#include <vector>

// This must be external since operator<< has to be defined
template <class Object, class Time, class KDS>
struct Trivial_event
{
  Trivial_event(){}
  template <class It, class Map>
  Trivial_event(It beg, It end, const Map &m, KDS* kds): kds_(kds) {
    for (; beg != end; ++beg) {
      objects_.push_back(m[*beg]);
    }
  }
  void process(Time t) const
  {
    std::cout << "At time " << t
	      << " the following objects are in the table: ";
    for (typename std::vector<Object>::const_iterator
	   cit= objects_.begin();
	 cit != objects_.end(); ++cit) {
      std::cout << *cit << " ";
    }
    std::cout << std::endl;
    kds_->set_processed(true);
  }

  std::vector<Object> objects_;
  typename KDS::Handle kds_;
};

template <class Object, class Time, class KDS>
std::ostream &operator<<(std::ostream &out,
			 const Trivial_event<Object, Time, KDS> &)
{
  out << "\"An event\"";
  return out;
}


/*!  This is a trivial kinetic data structure which doesn't actually
  do anything, other than print what it should be doing and maintain 1
  certificate in the event queue which has a list of all the moving
  objects in the moving object table.
*/
template <class Traits>
struct Trivial_kds: CGAL::Kinetic::Ref_counted<Trivial_kds<Traits> >
{
  typedef Trivial_kds<Traits> This;
  typedef typename Traits::Active_points_1_table::Data Point;
  typedef typename Traits::Simulator::Time Time;
  typedef typename Traits::Active_points_1_table::Key Point_key;
  typedef typename Traits::Simulator::Event_key Event_key;
  typedef CGAL::Kinetic::Active_objects_listener_helper<
    typename Traits::Active_points_1_table::Listener, This> Active_points_1_helper;
  typedef CGAL::Kinetic::Simulator_kds_listener<
    typename Traits::Simulator::Listener, This> Simulator_helper;

  typedef Trivial_event<Point, Time, This> Event;

  Trivial_kds(Traits tr): has_certificates_(true),
			  tr_(tr),
			  nth_(tr.active_points_1_table_handle(), this),
			  sh_(tr.simulator_handle(), this){}

  // this method is called with the value true when the event is processed
  void set_processed(bool tf) {
    if (tf== true) {
      event_= Event_key();
      set_has_certificates(false);
      set_has_certificates(true);
    }
  }

  void audit() const
  {
    /* In a real KDS you could use tr_.instantaneous_kernel() to build
       a static version of the kinetic data structure and check it for
       equality with the kinetic version.
    */
    CGAL_assertion(event_.is_valid());
    std::cout << "The structure is trivially correct at time "
	      << tr_.simulator_handle()->audit_time() << std::endl;
  }

  void set_has_certificates(bool tf) {
    typename Traits::Simulator::Handle sp= tr_.simulator_handle();
    if (has_certificates_ != tf) {
      has_certificates_=tf;
      if (has_certificates_) {
	bool ev= event_.is_valid();
	CGAL_assertion(!ev);
	Time t= CGAL::to_interval(sp->current_time()).second+1;
	event_= sp->new_event(t, Event(objects_.begin(),
				       objects_.end(),
				       *tr_.active_points_1_table_handle(),
				       this));
	std::cout << "Created event (" << event_ << ") at time " << t << std::endl;
      } else {
	if (event_.is_valid()) {
	  std::cout << "Deleting event " << event_ << std::endl;
	  sp->delete_event(event_);
	  event_=Event_key();
	}
      }
    }
  }

  bool has_certificates() const
  {
    // you can't use !event_ because the simulation
    // might end before the time you tried to schedule the event.
    return has_certificates_;
  }

  void insert(Point_key k) {
    std::cout << "Updating structure to include new object "
	      << k << "." << std::endl;
    objects_.insert(k);
    if (has_certificates_) {
      std::cout << "Updating all certificates which depend on "
                << k << "." << std::endl;
      set_has_certificates(false);
      set_has_certificates(true);
    }
  }

  void set(Point_key k) {
    std::cout << "An object changed, but we shouldn't every update"
	      << " structure when an object changes."<< std::endl;
    if (has_certificates_) {
      std::cout << "Updating all certificates which depend on "
                << k << "." << std::endl;
      set_has_certificates(false);
      set_has_certificates(true);
    }
  }

  void erase(Point_key k) {
    std::cout << "An object " << k << " was removed."<< std::endl;
    objects_.erase(k);
    if (has_certificates_) {
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
  Active_points_1_helper nth_;
  Simulator_helper sh_;
};

int main(int, char *[])
{
  typedef CGAL::Kinetic::Exact_simulation_traits_1 Traits;
  typedef Trivial_kds<Traits> TKDS;

  Traits tr(1,100);
  TKDS::Handle tk= new TKDS(tr);

  Traits::Simulator::Handle sp=tr.simulator_handle();

  Traits::Simulator::Function_kernel::Construct_function cf
    = tr.kinetic_kernel_object().function_kernel_object().construct_function_object();

  tk->set_has_certificates(true);

  typedef Traits::Kinetic_kernel::Point_2 Point;
  Traits::Active_points_1_table::Key k
    = tr.active_points_1_table_handle()->insert(Point(cf(1), cf(0)));
  sp->set_current_event_number(sp->current_event_number()+10);
  tr.active_points_1_table_handle()
    ->set(k, Traits::Kinetic_kernel::Point_2(cf(2), cf(0)));
  sp->set_current_event_number(sp->current_event_number()+10);
  tr.active_points_1_table_handle()->erase(k);
  sp->set_current_event_number(sp->current_event_number()+10);
  return EXIT_SUCCESS;
}
