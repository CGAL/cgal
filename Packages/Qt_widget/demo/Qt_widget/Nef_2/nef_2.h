#include <CGAL/Quotient.h>

template <class NT>
CGAL::Quotient<NT>
double_to_quotient(double x)
{ 
    NT num = 0; 
    NT den = 1;

    if (x != 0.0)
    { int neg = (x < 0);
      if (neg) x = -x;

      const unsigned shift = 15;   // a safe shift per step
      const unsigned int shift_pow = 32768; // = 2^shift
      const double width = 32768;  // = 2^shift
      const int maxiter = 20;      // ought not be necessary, but just in case,
                                   // max 300 bits of precision
      int expt;
      double mantissa = frexp(x, &expt);
      long exponent = expt;
      double intpart;
      int k = 0;
      
      while (mantissa != 0.0 && k++ < maxiter)

      { mantissa *= width; // shift double mantissa
        mantissa = modf(mantissa, &intpart);
        num *= shift_pow;
        num += (long)intpart;
        exponent -= shift;
      }
      int expsign = (exponent>0 ? +1 : (exponent<0 ? -1 : 0));
      exponent *= expsign;
      NT twopot(2);
      NT exppot(1);
      while (exponent!=0) {
        if (exponent & 1)
          exppot *= twopot;
        exponent >>= 1;
        twopot *= twopot;
      }

      if (expsign > 0)
        num *= exppot;
      else if (expsign < 0)
        den *= exppot;
      if (neg)
        num = -num;
    }
    CGAL::Quotient<NT> q(num,den);
    q.normalize();
    return q;
}