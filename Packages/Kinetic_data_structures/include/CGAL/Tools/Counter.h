#ifndef DSR_COUNTER_H
#define DSR_COUNTER_H
#include <iterator>
#include <CGAL/basic.h>
#include <CGAL/Tools/Label.h>

CGAL_BEGIN_NAMESPACE;
//! An integer iterator
/*!  It looks like an iterator but just wraps an integer. All the
  expected iterator operations are supported (namely operator* and
  operator->, which just return the value).
  
*/
template <class T=int, class W=T>
class Counter{
  typedef Counter<T, W> This;
public:
  Counter(T t=0):t_(t){
  }
  This &operator++(){
    ++t_;
    return *this;
  }
  This operator++(int){
    This ret= *this;
    ++t_;
    return ret;
  }
  //! Return the value
  W operator*() const {
    return W(t_);
  }
  //! Return a pointer to the value.
  W* operator->() const {
    return &W(t_);
  }
  bool operator==(const This &o) const {
    return t_ == o.t_;
  }
  bool operator!=(const This &o) const {
    return t_ != o.t_;
  }
  bool operator<(const This &o) const {
    return t_ < o.t_;
  }
  T operator-(const This &o) const {
    return t_-o.t_;
  }
  /*const This& operator+=(T &v) {
    t_+=v;
    return *this;
    }*/
  const This& operator+=(unsigned int v) {
    t_+=v;
    return *this;
  }
  const This& operator-=(T &v) {
    t_-=v;
    return *this;
  }
  //! Cast it to something
  template <class J>
  void operator()(const J&){
    ++*this;
  }
  typedef std::random_access_iterator_tag iterator_category;
  typedef W value_type;
  typedef T difference_type;
  typedef W* pointer;
  typedef W& reference;
protected:
  T t_;
};


//! The Counter is specialized for int and Label.
/*!
  I don't remember why I did this. 
*/
template <>
template <class LT>
class Counter<int, Label<LT> >{
  typedef Label<LT> Label;
  typedef Counter<int, Label > This;
public:

  Counter(int t=0):t_(t){
  }
  This &operator++(){
    t_= Label(t_.index()+1);
    return *this;
  }
  This operator++(int){
    This ret= *this;
    t_= Label(t_.index()+1);;
    return ret;
  }
  Label operator*() const {
    return t_;
  }
  const Label* operator->() const {
    return &t_;
  }
  bool operator==(const This &o) const {
    return t_ == o.t_;
  }
  bool operator!=(const This &o) const {
    return t_ != o.t_;
  }
  bool operator<(const This &o) const {
    return t_ < o.t_;
  }
  int operator-(const This &o) const {
    return t_.index()-o.t_.index();
  }
  /*const This& operator+=(T &v) {
    t_+=v;
    return *this;
    }*/
  const This& operator+=(unsigned int v) {
    t_= Label(t_.index()+v);
    return *this;
  }
  const This& operator-=(int &v) {
    t_ = Label(t_.index()-v);
    return *this;
  }
  template <class J>
  void operator()(const J&){
    ++*this;
  }
  typedef std::random_access_iterator_tag iterator_category;
  typedef Label value_type;
  typedef int difference_type;
  typedef Label* pointer;
  typedef Label& reference;
protected:
  Label t_;
};

template <class T>
Counter<T> counter(T t){
  return Counter<T>(t);
}
CGAL_END_NAMESPACE;
#endif
