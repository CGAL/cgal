#ifndef CGAL_LABEL_H_
#define CGAL_LABEL_H_
#include <CGAL/basic.h>
#include <iostream>
#include <sstream>


CGAL_BEGIN_NAMESPACE

//! A type which provides opaque, typed identifiers.
/*!
	Basically this is just an int with a type associated with it
	and one value (-1) picked out as the null value. It implements
	comparisons and things.
*/
template <class A_Type>
class Label {
protected:
  typedef Label<A_Type> This;
  int id_;
    
public:
  //! Construct it from an int
  Label(int i):id_(i){}
  //! Construct a default (null) label
  Label():id_(-1){
    if(0) print(); // make sure it is compiled
  }
  //! Make the next label
  This  next_label(){
    return Label(id_+1);
  }
  //! Return the null label
  operator bool() const  {
    return id_ != -1;
  }
  //! Convert to an index value, -1 is invalid
  int index() const {
    return id_;
  }
  bool operator<(const This &o) const {
    return id_ < o.id_;
  }
  bool operator>(const This &o) const {
    return id_ > o.id_;
  }
  bool operator==(const This &o) const {
    return id_ ==o.id_;
  }
  bool operator!=(const This &o) const {
    return id_ !=o.id_;
  }

  template <class OS>
  void write(OS &out) const {
    if (id_==-1) out << "N";
    else out << id_;
  }
  void print() const {
    write(std::cout);
  }
  //! Convert to a string
  std::string string() const {
    std::ostringstream os;
    write(os);
    return os.str();
  }
};


template <class A_type>
std::ostream& operator<<(std::ostream &out, const Label<A_type> &label){
  out << "(";
  label.write(out);
  out << ")";
  return out;
}
CGAL_END_NAMESPACE
#endif
