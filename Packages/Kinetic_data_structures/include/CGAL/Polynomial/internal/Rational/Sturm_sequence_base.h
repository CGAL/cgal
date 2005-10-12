#ifndef CGAL_POLYNOMIAL_INTERNAL_STURM_SEQUENCE_BASE_H
#define CGAL_POLYNOMIAL_INTERNAL_STURM_SEQUENCE_BASE_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Sign_variations_counter.h>
/*!
  \file Sturm_sequence.h A non-filtered Sturm sequence class.
*/

#ifndef POLYNOMIAL_NO_CGAL
#include <CGAL/MP_Float.h>
#endif

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE
template<class Kernel_t>
class Sturm_sequence_base
{
 public:
  typedef Kernel_t                    Kernel;
  typedef typename Kernel::Function   Polynomial;

 protected:
  typedef typename Kernel::Sign_at    Sign_at;
  typedef typename Polynomial::NT     NT;

  typedef std::vector<Polynomial>            Container;

  typedef CGAL_POLYNOMIAL_NS::Sign                Sign;
  typedef CGAL_POLYNOMIAL_NS::Comparison_result   Comparison_result;

  template<class Iterator>
  static
  unsigned int sign_variations(const Iterator& first,
			       const Iterator& beyond)
  {
    return Sign_variations_counter::sign_variations(first, beyond);
  }


  template<class NTRep>
  unsigned int sign_variations_base(const NTRep& x) const
  {
    Sign s0 = k_.sign_at_object( seq_[0] )(x);

    CGAL_exactness_precondition( s0 != CGAL::ZERO );

    std::vector<Sign> signs(size_);
    signs[0] = s0;

    for (unsigned int i = 1; i < size_; i++) {
      signs[i] = k_.sign_at_object( seq_[i] )(x);
    }

    return sign_variations(signs.begin(), signs.end());
  }


  void add(const Polynomial& f)
  {
    seq_.push_back(f);
  }


  template<class T>
  void normalize(Polynomial&, const T&)
  {
  }

#ifdef POLYNOMIAL_USE_CGAL  
  void normalize(Polynomial& r, const CGAL::MP_Float&)
  {
    // THE FOLLOWING HACK HAS BEEN DONE SO THAT MP_Float HOPEFULLY
    // DOES NOT RUN OUT OF EXPONENT BITS WHEN THE STURM SEQUENCE IS
    // COMPUTED
    NT half(0.5);
    while ( CGAL::abs(r[r.degree()]) > NT(1) ) {
      r = r * half;
    }

    NT two(2);
    while ( CGAL::abs(r[r.degree()]) < NT(2) ) {
      r = r * two;
    }
  }
#endif

public:

  Sturm_sequence_base() : size_(0) {}

  Sturm_sequence_base(const Polynomial& p, const Polynomial& q,
		      const Kernel &k)
    : size_(0), k_(k) {if (0) {Polynomial pq=p; pq=q;}}

  unsigned int size() const { return size_; }

  Polynomial operator[](unsigned int i) const
  {
    if (i >= seq_.size()) return zero_poly();
    else return Polynomial(seq_[i]);
  }

  // These are redundant
  Sign sign_at(const NT& x, unsigned int i) const
  {
    if ( i > size_ ) { return CGAL::ZERO; }
    return k_.sign_at_object(seq_[i])(x);
  }

  Sign sign_at_gcd(const NT& x) const
  {
    return k_.sign_at_object(seq_[size_-1])(x);
  }

  template<class T>
  unsigned int sign_variations(const T& x) const
  {
    return sign_variations_base(x);
  }

  void set_size(size_t sz) {
    seq_.resize(sz);
  }

protected:

  static const Polynomial &zero_poly(){
    static Polynomial zero(NT(0));
    return zero;
  }
  
  unsigned int size_;
  Container    seq_;
  Kernel k_;
};


CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif // CGAL_POLYNOMIAL_INTERNAL_STURM_SEQUENCE_BASE_H
