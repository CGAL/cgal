#ifndef CGAL_VORONOI_DIAGRAM_2_DUMMY_ITERATOR_H
#define CGAL_VORONOI_DIAGRAM_2_DUMMY_ITERATOR_H 1

#include <CGAL/Voronoi_diagram_adaptor_2/basic.h>
#include <CGAL/iterator.h>

CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

template<class Value_t>
class Dummy_iterator : public Emptyset_iterator
{
 private:
  typedef Dummy_iterator<Value_t>  Self;

 public:
  typedef Value_t      value_type;
  typedef value_type&  reference;
  typedef value_type*  pointer;

  Dummy_iterator() {}
  Dummy_iterator(const Dummy_iterator&) {}

  template< class T >
  Self& operator=(const T&) { return *this; }

  Self& operator++()        { return *this; }
  Self& operator++(int)     { return *this; }

  reference operator*()              { return *dummy_handle(); }
  pointer   operator->()             { return dummy_handle(); }
  const reference operator*() const  { return *dummy_handle(); }
  const pointer   operator->() const { return dummy_handle(); }

  bool operator==(const Self&) const { return true; }
  bool operator!=(const Self&) const { return false; }

 private:
  static value_type* dummy_handle() {
    static value_type dummy_handle_static;
    return &dummy_handle_static;
  }
};

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_DIAGRAM_2_DUMMY_ITERATOR_H
