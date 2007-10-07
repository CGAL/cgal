#ifndef CGAL_SPHERICAL_KERNEL_CONSTANT_H
#define CGAL_SPHERICAL_KERNEL_CONSTANT_H

namespace CGAL{
  typedef float HQ_NT;//type to represent the index of one hquadrant
  enum Fct_type{TAN, COT, FIXED, TAG_M2};

  inline Fct_type auto_ftype(const HQ_NT& hquad){
  if (hquad >7 || hquad<2 || (hquad>3 && hquad<6))
    return TAN;
  return COT;
  };    
}

#endif
