namespace Kinetic {

/*!
  \ingroup PkgKdsFrameworkConcepts
  \cgalConcept

  The class `Kinetic::Simulator` controls kinetic data structures by maintaining 
  a the current time and ensuring that events are processed when 
  necessary. 

  In addition, the `Kinetic::Simulator` can call on the kinetic data structures 
  to audit themselves at appropriate times. When the last event 
  processed and the next to be processed have different times, then 
  there is a rational value of time at which all kinetic data structures 
  should be non-degenerate (since there are no events at that time). At 
  such a time, kinetic data structures can easily verify their 
  correctness by checking that all the certificate predicates have the 
  correct value. When exactness checks are enabled, whenever the last 
  event processed and the next event to be processed have different 
  times, a 
  `Kinetic::Simulator::Listener::HAS_AUDIT_TIME` notification is made. Kinetic 
  data structures can listen for that event, and when it is made, they 
  can call `Kinetic::Simulator::audit_time()` to get the time value and 
  then verify that their structure is correct. 

  In addition, at such a time, the `Event::audit(Key)` is called on 
  each event. This allows kinetic data structures to check that the 
  event should be in the queue. 

  Typically, the simulator is created by the Kinetic:SimulationTraits 
  class and kinetic data structures request a handle to it from there. 

  Events that occur at or after `Kinetic::Simulator::end_time()` 
  will may not be processed. The exception are events which are 
  scheduled using the `new_final_event(Event)` call which are 
  guaranteed to occur after all other events (but have no particular 
  order amongst themselves). 

  \sa `CGAL::Kinetic::Simulator_objects_listener<Simulator_listener, KDS>`
  \sa `CGAL::Kinetic::Simulator_kds_listener<Simulator_listener, KDS>`

  \cgalHasModel `CGAL::Kinetic::Default_simulator<FunctionKernel, EventQueue>`

*/

class Simulator {
public:

  /// \name Types 
  /// @{

  /*!
    The type of the function kernel used to instantiate this `Kinetic::Simulator`. 
  */ 
  typedef unspecified_type Function_kernel; 

  /*!
    Extend this base class to listen to notifications from this `Kinetic::Simulator`. There are two types of notifications: `HAS_AUDIT_TIME` and `DIRECTION_OF_TIME`. The first is made when kinetic data structures can perform an audit. The second is made when the direction of time is changed. 
  */ 
  typedef unspecified_type Listener; 

  /*!
    \ingroup PkgKdsFrameworkOtherConcepts
    \cgalConcept

    The concept `Time` represents time in the simulator.

    Comparisons with other `Kinetic::Simulator::Time` objects are supported. 

    \cgalHasModel `double`
    \cgalHasModel `Kinetic::FunctionKernel::Root`

    \sa `Kinetic::Simulator`

  */
  class Time {
  public:
    /*!
      Construct an instance of time from a number 
      type, where NT is the number type used in the simulation. 
    */ 
    Time(NT); 

    /*!
      Write it to a stream. 
    */ 
    std::ostream& operator<<(std::ostream&, Time); 
  }; /* end Time */

/*!
  \ingroup PkgKdsFrameworkOtherConcepts
  \cgalConcept

  The concept `Event` represents a single event. Models of 
  `Event` should be passed to the `Kinetic::Simulator` when 
  scheduling events which will in turn pass them to the 
  `EventQueue`. 

  \cgalHasModel All over the place.
  \cgalHasModel `Kinetic::Event_base`

  \sa `Kinetic::EventQueue` 

  \cgalHeading{Example}

  All of the kinetic data structures provided have models of 
  `Event`. Here is the code implementing a swap event from the 
  sorting kinetic data structure. Events occuring at equal times are 
  perturbed so that the one that occurs first in the list is processed 
  first (just to illustrate the idea). 

  \code{.cpp} 

  template <class Certificate, class Id, class Root_enumerator> 
  class Swap_event { 
  typedef Swap_event<class Certificate, class Id, class Root_enumerator> This; 
  public: 
  Swap_event(Id o, Sort* sorter, 
  const Certificate &s): left_object_(o), 
  sorter_(sorter), 
  s_(s){} 
  void process(){ 
  sorter_->swap(left_object_, s_); 
  } 
  void *kds() const {return sorter_;} 
  CGAL::Comparison_result perturb_comparison(typename Sort::Event_key a, typename Sort::Event_key b) const { 
  return CGAL::compare(std::distance(sorter_->objects_begin(), left_object_), 
  std::distance(sorter_->objects_begin(), 
  sorter_->simulator_handle()->get_event<This>(b).left_object_)); 
  } 
  bool merge(typename Sort::Event_key a, typename Sort::Event_key b) { 
  return false; 
  } 
  Id left_object_; 
  Sort* sorter_; 
  Certificate s_; 
  }; 

  \endcode 

*/

  class Event {
  public:

/// \name Operations 
/// @{

/*!
  This method is called when the event 
  occurs. This method will only be called once per time this event is 
  scheduled and the event will be removed from the queue immediately 
  afterwards. 
*/ 
    void process(); 

/*!
  Return a `void *` which represents the KDS 
  which this event belongs to. The pointer is used solely to tell if 
  two events come from the same KDS for the purposes of handling 
  degeneracy. 
*/ 
    void* kds(); 

/*!
  The two events `a` and `b` occur at the same time 
  (`this` has key `a`). This method returns a 
  `CGAL::Comparison_result` which is used to order the two equal 
  events. If `CGAL::EQUAL` is returned then `merge` will be 
  called. 
*/ 
    CGAL::Comparison_result compare_concurrent(Key a, Key b) 
      const; 

/*!
  The two events 
  `a` and `b` occur at the same time (`this` has key 
  `a`) and cannot be perturbed to be unequal. This event allows 
  the KDS to merge event `b` with `a`. If it returns 
  `true` then `b` is dropped from the event queue. 
*/ 
    bool merge_concurrent(Key a, Key b); 

/*!
  Audit that this is a valid event. 
  To use this, kinetic data structure can check that this event is 
  indeed pointed to by the correct part of the combinatorial 
  structure. 
*/ 
    void audit(Key this_key); 

/*!
  Write the event to a stream. 
*/ 
    std::ostream& write(std::ostream&) const; 

/// @}

  }; /* end Event */



  /*!
    The basic number type used in computations. 
  */ 
  typedef unspecified_type NT; 

  /*!
    A reference counted pointer to be used for storing references to the object. 
  */ 
  typedef unspecified_type Handle; 

  /*!
    A reference counted pointer to be used for storing references to the object. 
  */ 
  typedef unspecified_type Const_handle; 

  /// @} 

  /// \name Creation 
  /// @{

  /*!
    Construct a `Kinetic::Simulator` which will process events between times start and end (events outside this window will be discarded). 
  */ 
  Simulator(const Time start=Time(0), const Time end= Time::infinity()); 

  /// @} 

  /// \name Operations 
  /// @{

  /*!
    Access the `Function_kernel` object used by the `Kinetic::Simulator`. 
  */ 
  Function_kernel function_kernel_object() const; 

  /*!
    Return the current time. 
  */ 
  Time current_time(); 

  /*!
    Set the current time to `t`, which cannot be less than `current_time`. Any events in the queue before time `t` are processed. 
  */ 
  void set_current_time(Time t); 

  /*!
    This function returns a time which can be represented using an instance of type `NT` which is shortly after the current time. You can then advance the current time to this one and act on the data structure using the return nt. 
  */ 
  NT next_time_reprsentable_as_nt() const; 

  /*!
    Return true if there is a rational number which is equivalent to the current time. Equivalent means that it has the same ordering relation to all previous and scheduled events. 
  */ 
  bool has_current_time_as_nt() const; 

  /*!

    Returns true if the current time is a rational number and there are no events at the current time. This means that the simulation can be audited at this time. 
  */ 
  bool has_audit_time() const; 

  /*!
    Return the time of the next event in the queue. 
  */ 
  Time next_event_time() const; 

  /*!
    Return the time the simulation will end. If time is running backwards, then this returns `Time::infinity()`. 
  */ 
  Time end_time() const; 

  /*!
    Set the current 
    time to \f$ t_{cur}\f$ and the end time to \f$ t_{end}\f$. The event queue 
    must be empty. Use this method if you want to reset or extend the 
    simulation. 
  */ 
  void set_interval(Time t_cur, Time t_end); 

  /*!
    Schedule a new event at time `t`. The object 
    `event` must implement the concept `Event`. The 
    `Event_key` returned can be used to access or deschedule the 
    event. 
  */ 
  template <class Event> Event_key new_event(Time t, const 
                                             Event event); 

  /*!
    Schedule a new event that will occur at the end of 
    the simulation. This type of event is useful if, for example, you 
    want to run for a while, change all motions, and then run some 
    more. 
  */ 
  template <class Event> Event_key new_final_event(const 
                                                   Event event); 

  /*!
    This method returns an 
    `Event_key` which is guaranteed never to be assigned to any real 
    event. This is a very useful placeholder for events which are known 
    never to occur (and allows data structures to differentiate between 
    uninitialized and never failing). 
  */ 
  Event_key null_event() const; 

  /*!
    Performs some checks as 
    to whether the key corresponds to a valid event. Generally, a event 
    is valid if it is not defaultly constructed and either is in the 
    queue or is the `null_event()`. 
  */ 
  void audit_event(Event_key) const; 

  /*!
    Return true if there are no more events. 
  */ 
  bool empty() const; 

  /*!
    Remove the event 
    referenced by `k` from the event queue. 
  */ 
  void delete_event(const Event_key k); 

  /*!
    This method returns a pointer to an event, which can be 
    used for recoving data, such as cached solvers, from that event. The 
    second argument really shouldn't be there, but gcc seems to 
    sometimes have issues if you try to specify the template value 
    directly. 
  */ 
  template <class Ev> typename 
  Queue::Event_handle<Ev>::Handle event(const Event_key k, const Ev 
                                        e) const; 

  /*!
    Return the time at which the event referenced by `k` occurs. 
  */ 
  Time event_time(Event_key k) const; 

  /*!
    Set the event referenced by key `k` to `ev`, for example if you want to change what happens when that event occurs. A new event key is returned. 
  */ 
  template <class Ev> 
  Event_key set_event(Event_key k, const Ev ev); 

  /*!
    Return `POSITIVE` if time is running forwards or `NEGATIVE` if it is running backwards. 
  */ 
  Sign direction_of_time() const; 

  /*!
    Set which direction time is running. 
  */ 
  void set_direction_of_time(Sign dir) const; 

  /*!
    Return the number of events which have been processed. 
  */ 
  unsigned int current_event_number() const; 

  /*!
    Process all events up to the ith event. `i` cannot be less than `current_event_number`. 
  */ 
  void set_current_event_number(unsigned int i) const; 

  /// @}


}; /* end Kinetic::Simulator */


/*!
  Return a double approximation of the time value. 
  \relates Simulator::Time
*/ 
double to_double(Simulator::Time); 

/*!
  Return an interval containing the time value. 
  \relates Simulator::Time
*/ 
std::pair<double, double> to_interval(Simulator::Time); 

} /* end namespace Kinetic */

