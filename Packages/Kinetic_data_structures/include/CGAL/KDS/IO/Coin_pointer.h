#ifndef CGAL_KDS_IO_COIN_POINTER_H
#define CGAL_KDS_IO_COIN_POINTER_H
#include <CGAL/basic.h>
#include <boost/intrusive_ptr.hpp>
#include <Inventor/nodes/SoNode.h>

//class SoNode;
void intrusive_ptr_add_ref(SoNode *n) {
  n->ref();
}
 
void intrusive_ptr_release(SoNode *n) {
  n->unref();
}

CGAL_BEGIN_NAMESPACE;


 
//! A reference counting pointer for storing pointers to Inventor objects.
/*!  Inventor objects already have reference counts built in, so I
  have to use the existing reference count.
*/
template <class T>
class Coin_pointer: public boost::intrusive_ptr<T> {
private:
  typedef boost::intrusive_ptr<T> P;
public:
  //! Pointer constructor
  Coin_pointer(T* t): P(t){}
  //! default constructor
  Coin_pointer(): P(){}
};

CGAL_END_NAMESPACE
#endif
