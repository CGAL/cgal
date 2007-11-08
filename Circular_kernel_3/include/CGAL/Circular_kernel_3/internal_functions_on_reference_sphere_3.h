#ifndef CGAL_SPHERICAL_KERNEL_PREDICATES_ON_REFERENCE_SPHERE_H
#define CGAL_SPHERICAL_KERNEL_PREDICATES_ON_REFERENCE_SPHERE_H

namespace CGAL {
  namespace SphericalFunctors {

//NOTE THAT THIS IMPLEMENTATION PROVIDE ONLY INTERSECTION POINT THAT ARE NOT A POLE
    
    
//private fonction
template <class SK, class OutputIterator>
OutputIterator intersect_3( const typename SK::Circle_on_reference_sphere_3& C1,
                                         const typename SK::Circle_on_reference_sphere_3& C2,
                                         CGAL::Inter_alg_info& IA,
                                         OutputIterator out)
{
  if (C1.type_of_circle_on_reference_sphere()==CGAL::POLAR && C2.type_of_circle_on_reference_sphere()==CGAL::POLAR &&
      CGAL::pole_of_polar_circle<SK>(C1)==CGAL::pole_of_polar_circle<SK>(C2)){//to handle 2 polar circles
    if (CGAL::collinear(CGAL::get_polar_coordinate<SK>(C1,CGAL::START_PT),CGAL::get_polar_coordinate<SK>(C2,CGAL::START_PT),typename SK::Point_3(0,0,0)))
      return out; //tangency case
    IA.is_polar=true;
  }
  
  std::vector<CGAL::Object> cont;
  typename SK::Intersect_3() (C1.supporting_sphere(),C2.supporting_sphere(),
    C1.reference_sphere(),std::back_inserter(cont));
  
  unsigned res=cont.size();
  if ( res<1 ) 
    return out;
  
  CGAL_precondition(!IA.is_polar || res!=1  );

  typename SK::Circular_arc_point_3 Cir_pts[2];
  std::pair<typename SK::Circular_arc_point_3,unsigned> p;
  CGAL::assign(p,cont[0]);
  Cir_pts[0]=p.first;
  if (res==2){
    CGAL::assign(p,cont[1]);
    Cir_pts[1]=p.first;
  }
  
  if (IA.qF==-1)//IA not initialized
    CGAL::set_IA <SK>(IA,Cir_pts,res==1);
  
  if (res==1 || IA.is_polar){//handle polar same pole intersecting
    if (IA.is_polar)
      *out++=typename SK::Circular_arc_point_on_reference_sphere_3(IA.qS,Cir_pts[IA.S_index]);
    else
        *out++=typename SK::Circular_arc_point_on_reference_sphere_3(IA.qF,Cir_pts[IA.F_index]);
  }
  else{
    *out++=typename SK::Circular_arc_point_on_reference_sphere_3(IA.qF,Cir_pts[IA.F_index]);
    *out++=typename SK::Circular_arc_point_on_reference_sphere_3(IA.qS,Cir_pts[IA.S_index]);
  }

  return out;
}

//return intersection point of two circles as Circular_arc_point_on_reference_sphere_3, the first correspond to FIRST_PT, second to SECOND_PT 
template <class SK, class OutputIterator>
inline OutputIterator intersect_3(const typename SK::Circle_on_reference_sphere_3& C1,
                                                   const typename SK::Circle_on_reference_sphere_3& C2,
                                                   OutputIterator out)
{
    CGAL::Inter_alg_info IA;
    return intersect_3<SK>(C1,C2,IA,out);
}

  
template <class SK, class OutputIterator>
OutputIterator intersect_3(const typename SK::Circle_on_reference_sphere_3& C1,
                                        const typename SK::Circle_on_reference_sphere_3& C2,
                                        CGAL::Inter_alg_info &IA,
                                        OutputIterator out, CGAL::EvtPt_num num)
{
  typename SK::Circular_arc_point_on_reference_sphere_3 Pts[2];
  intersect_3<SK>(C1,C2,IA,Pts);
  if (num!=CGAL::SECOND_PT || IA.is_polar)
    *out++=Pts[0];
  else
    *out++=Pts[1];
  if (num==CGAL::ALL)
    *out++=Pts[1];
  return out;
}


template<class SK>
bool has_on_ref_sphere(const typename SK::Half_circle_on_reference_sphere_3& H,const typename SK::Circular_arc_point_on_reference_sphere_3& P,bool point_is_on_circle){
  if (!point_is_on_circle){
    if (!SK().has_on_3_object()(H.supporting_circle(),P))
      return false;
  }
  int i=(H.get_position()==CGAL::UPPER)?(1):(-1);
  int tr=(H.supporting_circle().type_of_circle_on_reference_sphere()==CGAL::NORMAL)?(CGAL::compare(P.z(),H.supporting_circle().extremal_point_z() )):(i);//to handle threaded and polar
  return (i*tr>0);
}

template <class SK, class OutputIterator>
OutputIterator intersect_3( const typename SK::Half_circle_on_reference_sphere_3& H1,
                                        const typename SK::Half_circle_on_reference_sphere_3& H2,
                                        OutputIterator out){
  std::vector<typename SK::Circular_arc_point_on_reference_sphere_3> res;
  intersect_3<SK>(H1.supporting_circle(),H2.supporting_circle(),std::back_inserter(res));
  int i=res.size();
  if ( i>0 && has_on_ref_sphere<SK>(H1,res[0],true) &&  has_on_ref_sphere<SK>(H2,res[0],true) ) 
    *out++=typename SK::Circular_arc_point_on_reference_sphere_3(CGAL::hquadrant<SK>(res[0]),res[0]);
  if (i>1 && has_on_ref_sphere<SK>(H1,res[1],true) &&  has_on_ref_sphere<SK>(H2,res[1],true) )
    *out++=typename SK::Circular_arc_point_on_reference_sphere_3(CGAL::hquadrant<SK>(res[1]),res[1]);      
  return out;
}

template <class SK, class OutputIterator>
OutputIterator intersect_3( const typename SK::Circle_on_reference_sphere_3&        C,
                                         const typename SK::Half_circle_on_reference_sphere_3& H,
                                         OutputIterator out){
  std::vector<typename SK::Circular_arc_point_on_reference_sphere_3> res;
  intersect_3<SK>(C,H.supporting_circle(),std::back_inserter(res));
  int i=res.size();
  if ( i>0 && has_on_ref_sphere<SK>(H,res[0],true) ) 
    *out++=typename SK::Circular_arc_point_on_reference_sphere_3(CGAL::hquadrant<SK>(res[0]),res[0]);
  if (i>1 && has_on_ref_sphere<SK>(H,res[1],true) )
    *out++=typename SK::Circular_arc_point_on_reference_sphere_3(CGAL::hquadrant<SK>(res[1]),res[1]);      
  return out;
}

  
  }
}


#endif
