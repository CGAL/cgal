#ifndef CGAL_KDS_FREE_EVENT_BASE_H
#define CGAL_KDS_FREE_EVENT_BASE_H
#include <CGAL/Kinetic/basic.h>
CGAL_KINETIC_BEGIN_NAMESPACE

class Free_event_base {
public:

  Free_event_base(){}
	
  void* kds() const {
    return NULL;
  }
  
  template <class Key>
  CGAL::Comparison_result compare_concurrent(Key a,
					     Key b) const {
    if (a < b) return CGAL::SMALLER;
    else if (b < a) return CGAL::LARGER;
    else return CGAL::EQUAL;
  }	
  template <class Key>
  bool merge_concurrent(Key,
			Key) {
    return false;
  }
  template <class Key>
  void audit(Key)const{}
};

CGAL_KINETIC_END_NAMESPACE
#endif
