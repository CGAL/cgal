// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
//
// $URL: svn+ssh://sloriot@scm.gforge.inria.fr/svn/cgal/trunk/Circular_kernel_3/include/CGAL/Circular_kernel_3/Circular_arc_3.h $
// $Id: Circular_arc_3.h 40627 2007-10-16 15:00:59Z sloriot $
//
// Author(s) : Loriot Sebastien <Sebastien.Loriot@sophia.inria.fr>
//                   Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                   Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//                   Pedro Machado    <tashimir@gmail.com>
//                   Julien Hazebrouck
//                   Damien Leroy

#ifndef CGAL_SPHERICAL_KERNEL_CONSTANT_H
#define CGAL_SPHERICAL_KERNEL_CONSTANT_H

namespace CGAL{
  typedef float HQ_NT;//type to represent the index of one hquadrant
  enum Fct_type{TAN, COT, FIXED, TAG_M2};
  enum Circle_type{NORMAL,THREADED,POLAR,BIPOLAR};
  enum Hcircle_type{UPPER, LOWER, SENT_NPOLE=3, SENT_SPOLE,UNDEF};
  enum Pole_type{NOTAPOLE=-1,NPOLE=3,SPOLE=4};
  enum EvtPt_num{TANGENCY_PT,FIRST_PT,SECOND_PT,NONE,ALL,START_PT=1,END_PT,NORTH,SOUTH,UNKNOW,EVTPT_INI=-1};
  
  inline Fct_type auto_ftype(const HQ_NT& hquad){
  if (hquad >7 || hquad<2 || (hquad>3 && hquad<6))
    return TAN;
  return COT;
  };    

  template<class SK>
  typename CGAL::HQ_NT hquadrant(const typename SK::Circular_arc_point_3& R)
  {
    int x=CGAL::sign(R.x());
    int y=CGAL::sign(R.y());
    if (y>0){
      if (x>0)  switch (CGAL::sign(R.x()-R.y())){case -1: return 2; break;  case 0: return 1.5; break; case 1: return 1; break; }; 
      if (x<0) switch (CGAL::opposite(CGAL::sign(R.x()+R.y()))){case -1: return 3; break;  case 0: return 3.5; break; case 1: return 4; break; };
      return 2.5;//OPTI : we have more information here by x_y
    }
    else{
      if (y<0){
        if (x>0) switch (CGAL::sign(R.x()+R.y())){case -1: return 7; break;  case 0: return 7.5; break; case 1: return 8; break; };
        if (x<0) switch (CGAL::opposite(CGAL::sign(R.x()-R.y()))){case -1: return 6; break; case 0: return 5.5; break; case 1: return 5; break; };
        return 6.5;//OPTI : we have more information here by x_y
      }
      else{
        if (x>0)  return 0.5;
        if (x<0) return 4.5;
        return 0;
      }
    }
  }

  
  struct Inter_alg_info{
    CGAL::HQ_NT qF;
    CGAL::HQ_NT qS;
    int F_index;
    int S_index;
    bool is_polar;
    Inter_alg_info()
      :qF(-1),qS(-1),F_index(-1),S_index(-1),is_polar(false){};
    void print() const{
      std::cout << "(" << qF <<","<< qS <<") ("<< F_index <<","<< S_index <<") " << std::endl;
    }
  };

  //bool operator==(const Inter_alg_info& IA1,const Inter_alg_info& IA2){
  //  return (IA1.qF==IA2.qF) && (IA1.qS==IA2.qS) && (IA1.F_index==IA2.F_index) && (IA1.S_index==IA2.S_index);
  //}  
  
  
  template <class T>
  inline static const T& get_point(const T& p){return p;};
  
  template <class T>
  inline static const T& get_point(const std::pair<T,unsigned>& p){return p.first;}  
  
  template<class SK,class pt_container>
  static void init_indices(CGAL::Inter_alg_info& IA,const pt_container& Ipts,bool is_tangency){
    
    if (is_tangency){ //OPTI here adapt to have square free polynomial for tangency or do "if Root_of.poly().size()<3"
      IA.qF=hquadrant<SK>(get_point(Ipts[0]));
      IA.F_index=0;
      IA.S_index=0;
      IA.qS=IA.qF;      
      return;
    }
    
    CGAL::HQ_NT q0=hquadrant<SK>(get_point(Ipts[0]));
    CGAL::HQ_NT q1=hquadrant<SK>(get_point(Ipts[1]));

    
    CGAL::HQ_NT res=q0-q1;
    if (res == 0){//same quadrant
      res=CGAL::sign(get_point(Ipts[0]).y()*get_point(Ipts[1]).x()-get_point(Ipts[1]).y()*get_point(Ipts[0]).x());
      if (res==0){//same f(theta) value
        res=typename SK::Compare_z_3() (get_point(Ipts[1]),get_point(Ipts[0]));
        CGAL_precondition(res!=0);//DEBUG
      }
    }
    if (res <0){
      IA.qF=(q0==0)?(10):(q0);//hande 2 polar by same pole
      IA.qS=q1;
      IA.F_index=0;
      IA.S_index=1;
    }
    else{
      IA.qS=q0;
      IA.qF=(q1==0)?(10):(q1);//hande  2 polar by same pole
      IA.F_index=1;
      IA.S_index=0;
    }
  };
  
  template<class G>
  inline void exchange(G& v1,G& v2)
  {
    G tmp=v2;
    v2=v1;
    v1=tmp;
  }  
  
  //Never more than Pi between two intersection points
  template<class pt_container>
  static void set_inter_pt_conv(CGAL::Inter_alg_info& IA,const pt_container& Ipts){
    if (IA.qS-IA.qF>4){
      exchange(IA.F_index,IA.S_index);
      exchange(IA.qF,IA.qS);
    }
    else{
      if (!(IA.qS-IA.qF<4)){//cad =4  |  conflict who is the first?
        int s=CGAL::sign(get_point(Ipts[IA.F_index]).y()*get_point(Ipts[IA.S_index]).x()-get_point(Ipts[IA.S_index]).y()*get_point(Ipts[IA.F_index]).x());
        if (s>0){
          exchange(IA.F_index,IA.S_index);
          exchange(IA.qF,IA.qS);
        }
      }
    }
  }

  template<class SK,class pt_container>
  static void set_IA(CGAL::Inter_alg_info& IA,const pt_container& Ipts,bool is_tangency){
    init_indices<SK>(IA,Ipts,is_tangency);
    if (!is_tangency) set_inter_pt_conv(IA,Ipts);
  }  
  
  template<class FT,class Point_3>
  inline FT compute_a(const Point_3& c,const FT& R2,const FT& squared_radius){
    return c.x() * c.x() + c.y() * c.y()+ c.z() * c.z() + R2 - squared_radius;
  };
  
  template<class FT>
  FT circle_center_coefficent(const FT& x,const FT& y,const FT& z,const FT& r2,const FT& R2){
    return ((FT)(0.5) +  (R2 - r2)/(FT)(2* (x*x +y*y +z*z))) ;
  }
  
  template<class SK, class circle_on_sphere>
  CGAL::Circle_type classify_one_circle(const circle_on_sphere& C){
    if (C.supporting_sphere().center().z()==0){
      typename SK::Point_3 Pt=C.center();
      if (Pt.z()==0 && Pt.y()==0 && Pt.x()==0)
      return CGAL::BIPOLAR;
    }
    std::vector<CGAL::Object> cont;
    typename SK::Plane_3 Pl=SK().construct_plane_3_object()(typename SK::Algebraic_kernel::Polynomial_1_3(0,1,0,0));
    typename SK::Intersect_3()(C.reference_sphere(),C.supporting_sphere(),Pl,std::back_inserter(cont));
    
    switch (cont.size()){
      case 0: 
        return CGAL::NORMAL;
      case 2:{ 
        std::pair<typename SK::Circular_arc_point_3,unsigned> p1,p2;
        CGAL::assign(p1,cont[0]);CGAL::assign(p2,cont[1]);
        CGAL::Sign s1=CGAL::sign(p1.first.x());
        CGAL::Sign s2=CGAL::sign(p2.first.x());
        if (s1==CGAL::opposite(s2))
          return CGAL::THREADED;
        else
          if (s1!=s2) return CGAL::POLAR;
        }
        break;
      
      case 1:
        if (CGAL::abs(C.extremal_point_z())==C.reference_sphere().radius()) return CGAL::POLAR;
    }
    return CGAL::NORMAL;
  }  

  template <class SK>
  int Sign_power_of_pole(const typename SK::Circle_on_reference_sphere_3& C,CGAL::Pole_type P){
    typename SK::Point_3 center=C.supporting_sphere().center();
    typename SK::FT Part1=CGAL::squared_distance(typename SK::Point_3(0.,0.,0.),center)+C.reference_sphere().squared_radius()-C.supporting_sphere_squared_radius();
    typename SK::FT Part2=2*(P==CGAL::NPOLE?-1:1)*center.z();
    
    if (Part1<0 && Part2<0)
      return -1;
    if (Part1>0 && Part2>0)
      return 1;
    return (CGAL::sign(Part2)>0?-1:1)*CGAL::sign(CGAL::square(Part1)-CGAL::square(Part2)*C.reference_sphere().squared_radius());
  }    
  
  template <class SK>
  static CGAL::Pole_type pole_covered_by_supporting_sphere(const typename SK::Circle_on_reference_sphere_3& C){
    CGAL_precondition(C.type_of_circle_on_reference_sphere()==CGAL::THREADED);
    typename SK::FT NP=CGAL::Sign_power_of_pole<SK>(C,CGAL::NPOLE);
    if (NP < 0)
      return CGAL::NPOLE;
    else{
      if(CGAL::Sign_power_of_pole<SK>(C,CGAL::SPOLE) < 0 )
        return CGAL::SPOLE;
      return (NP==0)?(CGAL::NPOLE):(CGAL::SPOLE);
    }
  }

  template <class SK>
  CGAL::Pole_type pole_of_polar_circle(const typename SK::Circle_on_reference_sphere_3& C){
    CGAL_precondition(C.type_of_circle_on_reference_sphere()==CGAL::POLAR);
    return (CGAL::Sign_power_of_pole<SK>(C,CGAL::NPOLE)==0)?(CGAL::NPOLE):(CGAL::SPOLE);
  }    
  
    //indicate the direction of the tangent to a polar or bipolar circle according to a Start or end tag
  template <class SK>
  CGAL::Point_3<SK> get_polar_coordinate(const typename SK::Circle_on_reference_sphere_3& C, CGAL::EvtPt_num num){
    CGAL_precondition(C.type_of_circle_on_reference_sphere()==CGAL::POLAR || C.type_of_circle_on_reference_sphere()==CGAL::BIPOLAR);
    typename SK::FT s= ( C.circle_center_coefficient() <0 )?-1:1;
    typename SK::FT x=((num==CGAL::START_PT)?(1):(-1))*s*C.supporting_sphere_center().y();
    typename SK::FT y=((num==CGAL::END_PT)?(1):(-1))*s*C.supporting_sphere_center().x();
    return typename SK::Point_3(x,y,typename SK::FT(0));
  }
  
  //version with the radius known
  //~ template<class SK>
  //~ static int Sign_power_of_pole(const typename SK::Circle_on_reference_sphere_3& C,CGAL::Pole_type P){
    //~ typename SK::Point_3 center=C.supporting_sphere().center();
    //~ return CGAL::sign(my_pow(center.x(),2)+my_pow(center.y(),2)+my_pow(center.z()+ (P==CGAL::NPOLE?-1:1)*get_S_0_radius<P_NT>(),2)
                                    //~ -C.supporting_sphere_squared_radius());
  //~ }  
  
}

#endif
