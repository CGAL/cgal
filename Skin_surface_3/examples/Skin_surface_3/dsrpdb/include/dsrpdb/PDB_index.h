#ifndef DSR_PDB_PDB_INDEX_H
#define DSR_PDB_PDB_INDEX_H
#include <cassert>
#include <ostream>
namespace dsrpdb {

  //! A class for representing indices in a pdb.
  /*!  This class is here to keep such indices straight from norm 1
    based C++ indices. PDB indices start at 1 and can contain gaps.
  */
  template <class Type>
  class PDB_index {
    typedef PDB_index<Type> This;
  public:
    //! Construct an index from a non-zero integer.
    explicit PDB_index(unsigned int v): v_(v){}
    //! Construct an invalid index.
    PDB_index():v_(-1){}
    //! Return the integer.
    operator unsigned int() const {
      assert(*this);
      return v_;
    }
    bool operator==(This o) const {
      if (!o || !*this) return false;
      return v_== o.v_;
    }
    bool operator!=(This o) const {
      if (!o || !*this) return false;
      return v_!= o.v_;
    }
    bool operator<(This o) const {
      if (!o || !*this) return false;
      return v_ < o.v_;
    }
    bool operator>(This o) const {
      if (!o || !*this) return false;
      return v_ > o.v_;
    }
    bool operator<=(This o) const {
      if (!o || !*this) return false;
      return v_ <= o.v_;
    }
    bool operator>=(This o) const {
      if (!o || !*this) return false;
      return v_ >= o.v_;
    }
    operator bool() const {
      return v_ != -1;
    }
    std::ostream &write(std::ostream &o) const {
      if (operator bool()) {
	o << "(" << v_ << ")";
      } else {
	o << "(null)";
      }
      return o;
    }
  protected:
    int v_;
  };

  template <class T>
  std::ostream &operator<<(std::ostream &o, PDB_index<T> i) {
    return i.write(o);
  }
}
#endif
