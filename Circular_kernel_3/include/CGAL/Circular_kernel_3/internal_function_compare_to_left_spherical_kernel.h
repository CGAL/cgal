#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_COMPARE_TO_LEFT_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_COMPARE_TO_LEFT_H

namespace CGAL {
  namespace SphericalFunctors {

//Compare two Half_circle_on_reference_sphere_3 to the left of a common point
template <class SK>
struct trait_for_cmp_tgt{
  const typename SK::Algebraic_kernel::Root_for_spheres_2_3& Pt;
  typedef typename SK::Algebraic_kernel::Root_of_2 Tk_type;
  
  trait_for_cmp_tgt(const typename SK::Algebraic_kernel::Root_for_spheres_2_3& Pt):Pt(Pt){}

  typename SK::Algebraic_kernel::Root_of_2 unsigned_tkz_coeff_normal(const typename SK::Point_3& C,
                                                                                                                  const typename SK::FT& rk,
                                                                                                                  const typename SK::FT& R,
                                                                                                                  typename SK::FT (&gamma_k)(const typename SK::Point_3&,const typename SK::FT&,
                                                                                                                                                                            const typename SK::FT&) ) const
  {
    typename SK::FT gk=gamma_k(C,rk,R);
    return CGAL::sign(gk)*(C.x()*Pt.y()-C.y()*Pt.x());
  };

  typename SK::Algebraic_kernel::Root_of_2 unsigned_tkz_coeff_threaded(const typename SK::Point_3& C) const{
    return C.x()*Pt.y()-C.y()*Pt.x();
  };
    
};

// Special case: traits use to compare at theta = 0 
template <class SK>
struct trait_for_cmp_tgt_theta_0{
  typedef typename SK::Algebraic_kernel::Root_of_2 Tk_type;
  
  typename SK::FT unsigned_tkz_coeff_normal(const typename SK::Point_3& C,
                                                                                const typename SK::FT& rk,
                                                                                const typename SK::FT& R,
                                                                                typename SK::FT (&gamma_k)(const typename SK::Point_3&,const typename SK::FT&,const typename SK::FT&) ) const{
    typename SK::FT gk=gamma_k(C,rk,R);
    return -CGAL::sign(gk)*C.y();
  }

  typename SK::FT unsigned_tkz_coeff_threaded(const typename SK::Point_3& C) const{
    return -C.y();
  }
  
};


template<class SK,class Traits>
struct compare_to_hcircle_to_left{
  Traits traits;
  
  compare_to_hcircle_to_left(){}
  compare_to_hcircle_to_left(const typename SK::Algebraic_kernel::Root_for_spheres_2_3& Pt):traits(Pt){};
    
  
  typename SK::FT 
  gamma_k ( const typename  SK::Point_3& c,
                    const typename SK::FT& rk,
                    const typename SK::FT& R,
                    typename SK::FT& ak2) const
  {
    ak2=c.x()*c.x()+c.y()*c.y()+c.z()*c.z();
    return (ak2+R-rk)/(2.*ak2);
  };

    
  static typename SK::FT
  gamma_k_ (const typename SK::Point_3& c,
                    const typename SK::FT& rk,
                    const typename SK::FT& R)
  {
    typename SK::FT ak2=c.x()*c.x()+c.y()*c.y()+c.z()*c.z();
    return (ak2+R-rk)/(2.*ak2);
  }
  
  
  
  typename SK::FT square_norm_tk_normal (const typename SK::Point_3& c,const typename SK::FT& rk,const typename SK::FT& R) const;
  typename SK::FT square_norm_tk_threaded (const typename SK::Point_3& c,const typename SK::FT& rk,const typename SK::FT& R) const;
  void fill_tzk_n(const typename SK::Half_circle_on_reference_sphere_3& ES,typename Traits::Tk_type& tz,typename SK::FT& n) const;
  int sign_of_delta(const typename SK::Half_circle_on_reference_sphere_3& ES1,const typename SK::Half_circle_on_reference_sphere_3& ES2) const;
  typename SK::FT give_rk(const typename SK::Half_circle_on_reference_sphere_3& ES) const;
  int compare_for_delta_eq_0_threaded(const typename SK::Half_circle_on_reference_sphere_3& EST,const typename SK::Half_circle_on_reference_sphere_3& ES) const;
  int compare_for_delta_eq_0(const typename SK::Half_circle_on_reference_sphere_3& ES1,const typename SK::Half_circle_on_reference_sphere_3& ES2) const;
  CGAL::Comparison_result operator()(const typename SK::Half_circle_on_reference_sphere_3& ES1,const typename SK::Half_circle_on_reference_sphere_3& ES2) const; 
};


template<class SK,class Traits>
typename SK::FT compare_to_hcircle_to_left<SK,Traits>::square_norm_tk_normal(const typename SK::Point_3& c,const typename SK::FT& rk,const typename SK::FT& R) const{
  typename SK::FT ak2;
  typename SK::FT gk=gamma_k(c,rk,R,ak2);
  return ak2 * (R - gk * gk * ak2);
};

template<class SK,class Traits>
typename SK::FT compare_to_hcircle_to_left<SK,Traits>::square_norm_tk_threaded(const typename SK::Point_3& c,const typename SK::FT& rk,const typename SK::FT& R) const{
  typename SK::FT ak2;
  typename SK::FT gk=gamma_k(c,rk,R,ak2);
  return ak2 * (R - gk * gk * ak2);
};

template<class SK,class Traits>
void compare_to_hcircle_to_left<SK,Traits>::fill_tzk_n (const typename SK::Half_circle_on_reference_sphere_3& ES,typename Traits::Tk_type& tz,typename SK::FT& n) const{
  typedef typename SK::FT P_NT;
  if (ES.supporting_circle().type_of_circle_on_reference_sphere()==CGAL::THREADED){
    tz=traits.unsigned_tkz_coeff_threaded(ES.supporting_circle().supporting_sphere_center());
    tz*=(CGAL::pole_covered_by_supporting_sphere<SK>(ES.supporting_circle())==CGAL::NPOLE)?(1):(-1);
    n=square_norm_tk_threaded(ES.supporting_circle().supporting_sphere_center(),ES.supporting_circle().supporting_sphere_squared_radius(),ES.supporting_circle().reference_sphere().squared_radius());
  }
  else{
    tz=traits.unsigned_tkz_coeff_normal(ES.supporting_circle().supporting_sphere_center(),ES.supporting_circle().supporting_sphere_squared_radius(),ES.supporting_circle().reference_sphere().squared_radius(),compare_to_hcircle_to_left<SK,Traits>::gamma_k_);
    tz*=(ES.get_position()==CGAL::UPPER)?(-1):(1);
    n=square_norm_tk_normal(ES.supporting_circle().supporting_sphere_center(),ES.supporting_circle().supporting_sphere_squared_radius(),ES.supporting_circle().reference_sphere().squared_radius());    
  }
};


template<class SK,class Traits>
int compare_to_hcircle_to_left<SK,Traits>::sign_of_delta (const typename SK::Half_circle_on_reference_sphere_3& ES1,
  const typename SK::Half_circle_on_reference_sphere_3& ES2) const{
  typename Traits::Tk_type tz1,tz2;
  typename SK::FT n1,n2;
  fill_tzk_n(ES1,tz1,n1);
  fill_tzk_n(ES2,tz2,n2);
  if (tz1==0 && tz2==0)
    return 0;
  if (tz1==0)
    return CGAL::sign(tz2);
  if (tz2==0)
    return -CGAL::sign(tz1);
  
  if (CGAL::sign(tz1)!=CGAL::sign(tz2))
    return CGAL::sign(tz2);
  return CGAL::sign(tz2)*CGAL::sign(n1*CGAL::square(tz2)-n2*CGAL::square(tz1));
};


template<class SK,class Traits>
typename SK::FT compare_to_hcircle_to_left<SK,Traits>::give_rk(const typename SK::Half_circle_on_reference_sphere_3& ES) const{
  return ((ES.get_position()==CGAL::UPPER)?(1):(-1))*ES.supporting_circle().squared_radius();  
};

template<class SK,class Traits>
int compare_to_hcircle_to_left<SK,Traits>::compare_for_delta_eq_0_threaded(const typename SK::Half_circle_on_reference_sphere_3& EST,
  const typename SK::Half_circle_on_reference_sphere_3& ES) const{
  if (ES.supporting_circle().type_of_circle_on_reference_sphere()!=CGAL::THREADED)
    return (ES.get_position()==CGAL::UPPER)?(-1):(1);
  return (-CGAL::sign(EST.supporting_circle().center().z()-ES.supporting_circle().center().z()));
};
  
template<class SK,class Traits>
int compare_to_hcircle_to_left<SK,Traits>::compare_for_delta_eq_0(const typename SK::Half_circle_on_reference_sphere_3& ES1,
  const typename SK::Half_circle_on_reference_sphere_3& ES2) const{
  if (ES1.supporting_circle().type_of_circle_on_reference_sphere()==CGAL::THREADED)
    return compare_for_delta_eq_0_threaded(ES1,ES2);
  if (ES2.supporting_circle().type_of_circle_on_reference_sphere()==CGAL::THREADED)
    return -compare_for_delta_eq_0_threaded(ES2,ES1);
  typename SK::FT rc1=give_rk(ES1);
  typename SK::FT rc2=give_rk(ES2);
  CGAL_precondition(rc1!=0 && rc2!=0);
  if(CGAL::sign(rc1)*CGAL::sign(rc2)<0)
    return (CGAL::sign(rc1)>0)?(1):(-1);
  return (rc1<rc2)?(1):(-1);
};


template<class SK,class Traits>
CGAL::Comparison_result compare_to_hcircle_to_left<SK,Traits>::operator()(const typename SK::Half_circle_on_reference_sphere_3& ES1,
  const typename SK::Half_circle_on_reference_sphere_3& ES2) const{
  int res=-sign_of_delta(ES1,ES2);
  if (res==0)
    res=compare_for_delta_eq_0(ES1,ES2);
  CGAL_precondition(res!=0);
  if (res<0) return CGAL::LARGER;
  return CGAL::SMALLER;
};

}
}

#endif
