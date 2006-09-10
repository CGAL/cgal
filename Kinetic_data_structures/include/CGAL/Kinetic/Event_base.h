#ifndef CGAL_KDS_EVENT_BASE_H
#define CGAL_KDS_EVENT_BASE_H
#include <CGAL/Kinetic/basic.h>
CGAL_KINETIC_BEGIN_NAMESPACE

template <class KDS_ptr>
class Event_base {
public:
  typedef KDS_ptr KDS_handle;
  Event_base(KDS_ptr pt): kds_(pt){}
  Event_base(): kds_(NULL){}
	
  KDS_handle kds() const {
    return kds_;
  }
  KDS_handle kds() {
    return kds_;
  }
  template <class Key>
  CGAL::Comparison_result compare_concurrent(Key a,
					     Key b) const {
    return CGAL::compare(a,b);
  }	
  template <class Key>
  bool merge_concurrent(Key a,
			Key b) {
    return false;
  }
  template <class K>
  void audit(K k) const {}
protected:
  KDS_handle kds_;
};

CGAL_KINETIC_END_NAMESPACE
#endif
