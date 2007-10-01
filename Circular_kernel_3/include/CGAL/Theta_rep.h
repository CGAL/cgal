#ifndef CGAL_THETA_REP_H
#define CGAL_THETA_REP_H

namespace CGAL {

template < typename SphericalKernel >
class Theta_rep
  : public SphericalKernel::Kernel_base::Theta_rep
{
  typedef typename SphericalKernel::Kernel_base::Theta_rep 
                                           RTheta_rep;

  typedef typename SphericalKernel::Root_of_2             Root_of_2;
  typedef typename SphericalKernel::FT                    FT;
  typedef typename SphericalKernel::Algebraic_kernel      AK;
  typedef typename SphericalKernel::HQ_NT                 HQ_NT;
  typedef typename SphericalKernel::Root_for_spheres_2_3  Root_for_spheres_2_3;  

public:
  typedef SphericalKernel   R; 
  typedef RTheta_rep Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  Theta_rep(const RTheta_rep& p)
  :RTheta_rep(p){}    
    
    
  Theta_rep()
  : RTheta_rep(
    typename R::Construct_theta_rep()())
  {}

  Theta_rep(const HQ_NT& hq,const Root_of_2& r)
  : RTheta_rep(
    typename R::Construct_theta_rep()(hq,r))
  {}

    
  typename Qualified_result_of<typename R::Compute_theta_ftheta_3,Theta_rep>::type
  ftheta() const
  { return typename R::Compute_theta_ftheta_3()(*this);}
  
  typename Qualified_result_of<typename R::Compute_theta_hq_3,Theta_rep>::type
  hq() const
  { return typename R::Compute_theta_hq_3()(*this);}
  
    
};

}
#endif

