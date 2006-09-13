#ifndef DSR_GEOMETRY_IO_TRIPLE_H
#define DSR_GEOMETRY_IO_TRIPLE_H
#include <iostream>
template <class T>
class Unordered_triple {
public:
  typedef T value_type;
  Unordered_triple(){}
  Unordered_triple(T a, T b, T c){
    a_[0]=a; a_[1]=b; a_[2]=c;
    if (a_[0] > a_[1]) std::swap(a_[0], a_[1]);
    if (a_[1] > a_[2]) std::swap(a_[1], a_[2]);
    if (a_[0] > a_[1]) std::swap(a_[0], a_[1]);
  }
  T first() const {
    return a_[0];
  }
  T second() const {
    return a_[1];
  }
  T third() const {
    return a_[2];
  }
  bool operator<(const Unordered_triple &o) const {
    for (unsigned int i=0; i< 3; ++i){
      if (a_[i] < o.a_[i]) return true;
      else if (a_[i] > o.a_[i]) return false;
    }
    return false;
  }
  bool operator>(const Unordered_triple &o) const {
    for (unsigned int i=0; i< 3; ++i){
      if (a_[i] > o.a_[i]) return true;
      else if (a_[i] < o.a_[i]) return false;
    }
    return false;
  }
  bool operator==(const Unordered_triple &o) const {
    return a_[0]==o.a_[0] && a_[1]==o.a_[1] && a_[2]== o.a_[2];
  }
  bool operator!=(const Unordered_triple &o) const {
    return a_[0]!=o.a_[0] || a_[1]!=o.a_[1] || a_[2]!= o.a_[2];
  }
private:
  T a_[3];
};



template <class T>
inline std::ostream &operator<<(std::ostream &o, Unordered_triple<T> e){
  o << "{" << e.first() << ", " << e.second() << ", " << e.third() << "}";
  return o;
}

template <class T>
inline std::istream &operator>>(std::istream &i, Unordered_triple<T> &e){
  char c;
  i >> c;
  if (c != '{') {
    std::cerr << "Error reading triple, found " << c << " when expecting '{'." << std::endl;
    return i;
  }
  T f,s,t;
  i >> f;
  i >> c;
  if (c != ',') {
    std::cerr << "Error reading triple, found " << c << " when expecting ','." << std::endl;
    return i;
  }
  i >> s;
  i >> c;
  if (c != ',') {
    std::cerr << "Error reading triple, found " << c << " when expecting ','." << std::endl;
    return i;
  }
  i >> t;
  i >> c;
  if (c != '}') {
    std::cerr << "Error reading triple, found " << c << " when expecting '}'." << std::endl;
    return i;
  }
  e= Unordered_triple<T>(f,s,t);
  return i;
}
#endif
