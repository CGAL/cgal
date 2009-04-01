#ifndef CGAL_PDB_NESTED_ITERATOR
#define CGAL_PDB_NESTED_ITERATOR
#include <CGAL/PDB/basic.h>
#include <boost/tuple/tuple.hpp>
#include <CGAL/assertions.h>

namespace CGAL { namespace PDB {
namespace internal {
  //! An iterator through the atoms of a CGAL::PDB::Chain.
  template <class T>
  class Nested_iterator {
    typedef Nested_iterator<T> This;
    typedef typename T::Make_value Make_value;
    typedef typename T::Get_inner Get_inner;
    typedef typename T::Inner Inner;
    typedef typename T::Outer Outer;
  public:
    typedef typename T::value_type value_type;
    typedef std::forward_iterator_tag iterator_category;
    typedef std::size_t difference_type;
    typedef value_type reference;
    typedef value_type* pointer;

    reference operator*() {
      return ret_;
    }
    const reference operator*() const {
      return ret_;
    }
    pointer operator->() {
      return &ret_;
    }
    const This& operator++() {
      CGAL_assertion(ait_ != aend_);
      ++ait_;
      while (ait_== aend_) {
	++rit_;
	if (rit_== rend_) {
          break;
        }
	Inner in= Get_inner()(rit_);
        ait_=in.begin();
        aend_= in.end();
      }
      if (rit_ != rend_) {
	ret_=Make_value()(rit_, ait_);
      }
      return *this;
    }
    
    bool operator==(const This& o) const {
      if (rit_ == rend_) return rit_==o.rit_;
      else return rit_== o.rit_ && ait_ == o.ait_;
    }
    bool operator!=(const This& o) const {
      if (rit_== rend_) return rit_!= o.rit_;
      else return rit_!= o.rit_ || ait_ != o.ait_;
    }
    Nested_iterator(){}
    template <class OT>
    Nested_iterator(const Nested_iterator<OT> &o) {
      copy_from(o);
    }
    template <class OT>
    This& operator=(const Nested_iterator<OT> &o) {
      copy_from(o);
      return *this;
    }
    Nested_iterator(const This &o) {
      copy_from(o);
    }
    This& operator=(const This &o) {
      copy_from(o);
      return *this;
    }

    template <class R>
    Nested_iterator(R r): rit_(r.begin()), rend_(r.end()){
      if (!r.empty()) {
	Inner in= Get_inner()(rit_);
        ait_=in.begin();
        aend_= in.end();
	//typename boost::tuple_element<1, value_type>::type st=ait_->second;
	ret_= Make_value()(rit_, ait_);
      }
    }

  protected:
    template <class OT>
    void copy_from(const Nested_iterator<OT> &o) {
      rit_= o.rit_;
      rend_= o.rend_;
      if (rit_!= rend_) {
	ait_=o.ait_;
	aend_=o.aend_;
	ret_= o.ret_;
      }
    }

  
    
    typename T::Outer::iterator rit_, rend_;
    typename T::Inner::iterator ait_, aend_;
    value_type ret_; 
  };
}
}}
#endif
