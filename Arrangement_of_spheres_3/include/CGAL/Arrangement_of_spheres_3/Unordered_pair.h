#ifndef DSR_GEOMETRY_IO_EDGE_H
#define DSR_GEOMETRY_IO_EDGE_H
#include <iostream>

template <class T>
class Unordered_pair {
public:
  typedef T value_type;
  Unordered_pair(){}
  Unordered_pair(T a, T b): a_(std::min(a,b)), b_(std::max(a,b)){}
  Unordered_pair(const std::pair<T,T> &p): a_(std::min(p.first,p.second)), b_(std::max(p.first,p.second)){}
  T first() const {
    return a_;
  }
  T second() const {
    return b_;
  }
  bool operator<(const Unordered_pair &o) const {
    if (a_ < o.a_) return true;
    else if (a_ > o.a_) return false;
    else return b_ < o.b_;
  }
  bool operator>(const Unordered_pair &o) const {
    if (a_ > o.a_) return true;
    else if (a_ < o.a_) return false;
    else return b_ > o.b_;
  }
  bool operator==(const Unordered_pair &o) const {
    return a_==o.a_ && b_==o.b_;
  }
  bool operator!=(const Unordered_pair &o) const {
    return a_!=o.a_ || b_!=o.b_;
  }
private:
  T a_, b_;
};



template <class T>
inline std::ostream &operator<<(std::ostream &o, Unordered_pair<T> e){
  o << "{" << e.first() << ", " << e.second() << "}";
  return o;
}

template <class T>
inline std::istream &operator>>(std::istream &i, Unordered_pair<T> &e){
  char c;
  i >> c;
  if (c != '{') {
    std::cerr << "Error reading pair, found " << c << " when expecting '{'." << std::endl;
    return i;
  }
  T f,s;
  i >> f;
  i >> c;
  if (c != ',') {
    std::cerr << "Error reading pair, found " << c << " when expecting ','." << std::endl;
    return i;
  }
  i >> s;
  i >> c;
  if (c != '}') {
    std::cerr << "Error reading pair, found " << c << " when expecting '}'." << std::endl;
    return i;
  }
  e= Unordered_pair<T>(f,s);
  return i;
}
#endif
