#ifndef CGAL_GENERIC_HANDLE_MAP_H
#define CGAL_GENERIC_HANDLE_MAP_H

#include <CGAL/Unique_hash_map.h>

CGAL_BEGIN_NAMESPACE

struct Void_handle_hash_function {
    std::size_t operator() (void* h) const { 
        return std::size_t(h);
    }
};


template <class I>
class Generic_handle_map : public 
  Unique_hash_map<void*,I,Void_handle_hash_function>
{ typedef Unique_hash_map<void*,I,Void_handle_hash_function> Base;
public:
  Generic_handle_map() : Base() {}
  Generic_handle_map(I i) : Base(i) {}

  template <class H>
  const I& operator[](H h) const 
  { return Base::operator[](&*h); }

  template <class H>
  I& operator[](H h) 
  { return Base::operator[](&*h); }

};

CGAL_END_NAMESPACE
#endif //CGAL_GENERIC_HANDLE_MAP_H
