#ifndef CGAL_POLYNOMIAL_INTERNAL_DERIVATIVE_H
#define CGAL_POLYNOMIAL_INTERNAL_DERIVATIVE_H

#include <CGAL/Polynomial/basic.h>
#include <vector>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <class Fn>
struct Derivative{
private:
  typedef typename Fn::iterator Cit;
  struct It {
    It(Cit it, int i): i_(i), cit_(it) {}
    It(){}
    
    typedef typename Cit::iterator_category iterator_category;
    typedef typename Cit::value_type value_type;
    typedef typename Cit::pointer pointer;
    typedef typename Cit::reference reference;
    typedef typename Cit::difference_type difference_type;

    value_type operator*() const {
      return value_type(i_)**cit_;
    }
    pointer operator->() const {
      return value_type(i_)**cit_;
    }
    It operator++(){
      ++i_;
      ++cit_;
      return *this;
    }
    It operator++(int){
      It t=*this;
      ++i_;
      ++cit_;
      return t;
    }
    It operator--(){
      --i_;
      --cit_;
      return *this;
    }
    It operator--(int){
      It t=*this;
      --i_;
      --cit_;
      return t;
    }
    bool operator==(It o) const {
      return cit_ == o.cit_;
    }
    bool operator!=(It o) const {
      return cit_ != o.cit_;
    }
    difference_type operator-(It o) const {
      return cit_-o.cit_;
    }
    It operator+=(difference_type i) {
      cit_+= i;
      i_+= i;
      return *this;
    }
  protected:
    int i_;
    Cit cit_;
  };
public:
  typedef Fn result_type;
  typedef Fn argument_type;

  result_type operator()(const argument_type &o) const {
    if (o.is_constant()) { return result_type(typename result_type::NT(0)); }
    else {
      It b(o.begin(),0);
      ++b;
      return result_type(b, It(o.end(), o.degree()+1));
    }
  }

  void write(std::ostream &out) const {
    out << "diff";
  }
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
