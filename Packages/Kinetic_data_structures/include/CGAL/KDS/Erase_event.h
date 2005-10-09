#ifndef CGAL_KDS_MOVING_OBJECT_ERASER_H_
#define CGAL_KDS_MOVING_OBJECT_ERASER_H_
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE;


//! Delete a single moving object from the MOT at a particular time. 
/*!
  Note that this class has not been used. 
*/
template <class MOT>
class Erase_event{
  typedef typename MOT::Pointer Pointer;
  typedef typename MOT::Key Key;
public:
  Erase_event(Key k, 
			Pointer mot):mot_(mot),
				     k_(k){}
  template <class T>
  void process(const T&) {
    CGAL_KDS_LOG(LOG_SOME,"Deleting object.\n");
    mot_->erase(k_);
  }
  void write(std::ostream &out) const {
    out << "E" << k_;
  }
protected:
  Pointer mot_;
  Key k_;
};

template <class MH>
std::ostream &operator<<(std::ostream &out, const Erase_event<MH> &moi){
  moi.write(out);
  return out;
}



CGAL_KDS_END_NAMESPACE;
#endif
