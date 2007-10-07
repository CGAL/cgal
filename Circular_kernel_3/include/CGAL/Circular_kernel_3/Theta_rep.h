#ifndef CGAL_SPHERICAL_KERNEL_THETA_REP_H
#define CGAL_SPHERICAL_KERNEL_THETA_REP_H

namespace CGAL {
  namespace CGALi {

template <class SK>
class Theta_rep{
  typedef typename SK::Algebraic_kernel::Root_of_2 Root_of_2;
  typedef typename SK::FT FT;    
  typedef std::pair<HQ_NT,Root_of_2> Rep;
  
  //
  typename SK::template Handle<Rep>::type  base;
  public:
  Theta_rep(const HQ_NT& hq,const Root_of_2& r):base(std::pair<HQ_NT,Root_of_2>(hq,r)){}
  Theta_rep(){}
  const HQ_NT& hq() const {return get(base).first;}
  const Root_of_2& ftheta() const {return CGAL::get(base).second;}
};    

  }
}

#endif
