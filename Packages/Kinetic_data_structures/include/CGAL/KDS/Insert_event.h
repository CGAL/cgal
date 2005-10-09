#ifndef CGAL_KDS_INSERT_EVENT_H_
#define CGAL_KDS_INSERT_EVENT_H_
#include <CGAL/KDS/basic.h>
#include <CGAL/KDS/internal/To_static.h>

CGAL_KDS_BEGIN_NAMESPACE;


//! An event to insert a single object into a MovingObjectTable

template <class MOT>
class Insert_event {
  typedef typename MOT::Pointer Pointer;
  typedef typename MOT::Data Object;
public:
  //! Construct it with the time, the object and where to put it
  Insert_event(const Object &obj, 
	       Pointer mot):mot_(mot),
			    obj_(obj){}
  template <class T>
  void process(const T&) {
    CGAL_KDS_LOG(LOG_SOME, "Inserting object.\n");
    mot_->insert(obj_);
  }

   void write(std::ostream &out) const {
     out << " I" << obj_ ;
  }
protected:
  Pointer mot_;
  Object obj_;
};

template <class MH>
std::ostream &operator<<(std::ostream &out, const Insert_event< MH> &moi){
  moi.write(out);
  return out;
}



CGAL_KDS_END_NAMESPACE;
#endif
