#ifndef CGAL_PMWX_SWEEP_LINE_FUNCTORS_H
#define CGAL_PMWX_SWEEP_LINE_FUNCTORS_H

CGAL_BEGIN_NAMESPACE

template <class Event, class SweepLineTraits_2>
class Event_less_functor 
{
public:
  typedef SweepLineTraits_2 Traits;
  
  Event_less_functor( Traits *traits) : m_traits(traits) {}
  
  bool operator()(  Event* e1, Event* e2) const 
  { 
    return (m_traits->compare_xy(e1->get_point(), e2->get_point()) == SMALLER);
  }

private:

  /*! a pointer to a traits class */
  Traits *m_traits;
};

CGAL_END_NAMESPACE

#endif
