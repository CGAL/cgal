#ifndef CGAL_POLYNOMIAL_POLYNOMIAL_H_
#define CGAL_POLYNOMIAL_POLYNOMIAL_H_
#include <CGAL/Polynomial/basic.h>

#include <CGAL/Polynomial/internal/interval_arithmetic.h>
#include <CGAL/Polynomial/internal/Polynomial_impl.h>
//#include <utility>
#include <sstream>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE


//! A basic polynomial class
/*!  The implementation is proved by internal::Polynomial_impl. This
  strips leading 0s. When debugging is on, a string representation of
  the polynomial is stored. This is fairly key for debugging, but
  rather slow.
*/
template <class NTT>
class Polynomial: public internal::Polynomial_impl<Polynomial<NTT>, NTT> {
  typedef internal::Polynomial_impl<Polynomial<NTT>, NTT>  Parent;
  typedef Polynomial<NTT> This;

  friend class internal::Polynomial_impl<Polynomial<NTT>, NTT>; // NOT SO CLEAN

#ifndef NDEBUG
  typedef std::string Approximation;
  void generate_approximation() const  {
    std::ostringstream s;
    this->write(s);
    approximation_= s.str();
    //return s.str();
  }

#endif

  void approximate() const {
#ifndef NDEBUG
    generate_approximation();
#endif
  }
public:
  //================
  // CONSTRUCTORS
  //================

  //! Default
  Polynomial(){
#ifndef NDEBUG
    approximation_="Not initialized.";
#endif
  }

  //! Make a constant polynomial
  Polynomial(const NTT& c): Parent(c){
#ifndef NDEBUG
    approximate();
#endif
    strip_leading_zeros();
  }

  //! Make a polynomial from an iterator range
  template<typename Iterator>
  Polynomial(Iterator first, Iterator beyond)
    : Parent(first,beyond) {
#ifndef NDEBUG
    approximate();
#endif
    strip_leading_zeros();
  }

private:
  Polynomial(const Parent &p): Parent(p){
#ifndef NDEBUG
    approximate();
#endif
    strip_leading_zeros();
  }

protected:

  void strip_leading_zeros() {
    if ( this->is_zero() ) { return; }

    do {
      Sign s = sign( this->coefs_[this->degree()] );
      if ( s == ZERO ) {
	Polynomial_assertion( this->coefs_.size() > 0 );
	this->coefs_.resize(this->coefs_.size() - 1);
      } else {
	break;
      }
    } while ( !this->is_zero() );
  }


  void finalize() {
    strip_leading_zeros();
  }

private:
#ifndef NDEBUG
  /*! A string represneting the approximation of the polynomial. For
    inspection in the debugger.*/
  mutable Approximation approximation_;
#endif
};


CGAL_POLYNOMIAL_END_NAMESPACE

#endif
