#ifndef CGAL_POLYNOMIAL_KERNEL_NEGATE_VARIABLE_H
#define CGAL_POLYNOMIAL_KERNEL_NEGATE_VARIABLE_H

#include <CGAL/Polynomial/basic.h>
#include <vector>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE
//------------------------------------------------------------------


template <class Polynomial>
class Negate_variable {
public:
  typedef typename Polynomial::NT NT;
  typedef NT          argument_type;
  typedef Polynomial  result_type;

  Negate_variable() {}

  void write(std::ostream &out) const {
    out << "negate_var";
  }

  result_type operator()(const Polynomial &f) const {
    int size = f.degree() + 1;
    std::vector<NT> coefs(size);
    
    for (int i = 0; i < size; i++) {
      if (i%2 == 1) {
	coefs[i]= -f[i];
      } else {
	coefs[i]= f[i];
      }
    }

    Polynomial ret(coefs.begin(), coefs.end());

    Polynomial_assertion(ret.degree() == f.degree());

    return ret;
  }
};


CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif // CGAL_POLYNOMIAL_KERNEL_NEGATE_VARIABLE_H
