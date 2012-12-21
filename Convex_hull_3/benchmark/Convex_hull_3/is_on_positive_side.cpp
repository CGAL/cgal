#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/random_selection.h>


#include <iostream>

namespace CGAL{

template <class Kernel,int alternatives>
struct Is_on_positive_side_of_plane_3;

template <class Kernel>
struct Is_on_positive_side_of_plane_3<Kernel,0>{
  typedef typename Kernel::Point_3 Point_3;
  Point_3 p,q,r;
  typename Kernel::Orientation_3 orientation;
public:  
  Is_on_positive_side_of_plane_3(const Kernel& kernel,const Point_3& p_,const Point_3& q_,const Point_3& r_)
  :p(p_),q(q_),r(r_),orientation(kernel.orientation_3_object()) {}
    
  bool operator() (const Point_3& s) const 
  {
    return orientation(p,q,r,s) == CGAL::POSITIVE;
  }
};
  
template <class Kernel>
struct Is_on_positive_side_of_plane_3<Kernel,1>{
  typedef typename Kernel::Point_3                                              Point_3;
  double m10,m20,m21,Maxx,Maxy,Maxz;
  const Point_3& p_;
  
  int static_filtered(double psx,double psy, double psz) const{
    using std::fabs;

    // Then semi-static filter.
    double maxx = Maxx;
    if (maxx < fabs(psx)) maxx = fabs(psx);
    double maxy = Maxy;
    if (maxy < fabs(psy)) maxy = fabs(psy);
    double maxz = Maxz;
    if (maxz < fabs(psz)) maxz = fabs(psz);
    double det =  psx*m10 - m20*psy + m21*psz;
    
    // Sort maxx < maxy < maxz.
    if (maxx > maxz)
        std::swap(maxx, maxz);
    if (maxy > maxz)
        std::swap(maxy, maxz);
    else if (maxy < maxx)
        std::swap(maxx, maxy);

    // Protect against underflow in the computation of eps.
    if (maxx < 1e-97) /* cbrt(min_double/eps) */ {
      if (maxx == 0)
        return 0;
    }
    // Protect against overflow in the computation of det.
    else if (maxz < 1e102) /* cbrt(max_double [hadamard]/4) */ {
      double eps = 5.1107127829973299e-15 * maxx * maxy * maxz;
      if (det > eps)  return 1;
      if (det < -eps) return -1;
    }
    return 555;
  }
public:
  Is_on_positive_side_of_plane_3(const Kernel&,const Point_3& p,const Point_3& q,const Point_3& r):p_(p)
  {
    double pqx = q.x() - p.x();
    double pqy = q.y() - p.y();
    double pqz = q.z() - p.z();
    double prx = r.x() - p.x();
    double pry = r.y() - p.y();
    double prz = r.z() - p.z();   

    m10 = pqy*prz - pry*pqz;
    m20 = pqx*prz - prx*pqz;
    m21 = pqx*pry - prx*pqy;
    
    Maxx = fabs(pqx);
    if (Maxx < fabs(prx)) Maxx = fabs(prx);
    Maxy = fabs(pqy);
    if (Maxy < fabs(pry)) Maxy = fabs(pry);
    Maxz = fabs(pqz);
    if (Maxz < fabs(prz)) Maxz = fabs(prz);
  }
  
  bool operator() (const Point_3& s) const 
  {
    double psx = s.x() - p_.x();
    double psy = s.y() - p_.y();
    double psz = s.z() - p_.z(); 
    
    int static_res = static_filtered(psx,psy,psz);
    if (static_res != 555)
      return static_res == 1;
    
    std::cerr << "ERROR static predicate failure!!!\n";
    exit(EXIT_FAILURE);
  }
};


template <class Kernel>
struct Is_on_positive_side_of_plane_3<Kernel,2>{
  typedef Simple_cartesian<Interval_nt_advanced >                               CK;  
  typedef typename Kernel::Point_3                                              Point_3;
  
  Cartesian_converter<Kernel,CK>                        to_CK;
  typename CK::Plane_3 ck_plane;

  Is_on_positive_side_of_plane_3(const Kernel&,const Point_3& p,const Point_3& q,const Point_3& r):
    ck_plane(to_CK(p),to_CK(q),to_CK(r))
  {}
    
  bool operator() (const Point_3& s) const 
  {
    try{
      return ck_plane.has_on_positive_side(to_CK(s));
    }
    catch (Uncertain_conversion_exception){
      std::cerr << "ERROR Interval filtering failure\n";
      exit(EXIT_FAILURE);
    }
  }
};

}//namespace CGAL


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Is_on_positive_side_of_plane_3<K,0>  Pred_using_kernel;
typedef CGAL::Is_on_positive_side_of_plane_3<K,1>  Pred_using_optimized_static;
typedef CGAL::Is_on_positive_side_of_plane_3<K,2>  Pred_using_intervals;


template <class Predicate,class Outputiterator>
void run(const K::Point_3& p,const K::Point_3& q,const K::Point_3& r,const std::vector<K::Point_3>& queries,Outputiterator out){
  CGAL::Timer time; time.start();
  Predicate predicate(K(),p,q,r);
  for (std::vector<K::Point_3>::const_iterator it=queries.begin();it!=queries.end();++it)
    *out++=predicate(*it);
  time.stop();
  std::cout << time.time() << std::endl;
}

int main()
{
  std::size_t nb_pts=20000000;
  typedef CGAL::Creator_uniform_3<double,K::Point_3>  Creator;
  CGAL::Random_points_in_sphere_3<K::Point_3,Creator>     gen(1);
  
  std::vector<K::Point_3> points;
  points.reserve(nb_pts);
  
  K::Point_3 p=*gen++,q=*gen++,r=*gen++;
  CGAL::cpp11::copy_n(gen,nb_pts,std::back_inserter(points));
  
  std::vector<bool> res0; res0.reserve(nb_pts);
  std::vector<bool> res1; res1.reserve(nb_pts);
  std::vector<bool> res2; res2.reserve(nb_pts);
  
  std::cout << "Running kernel predicates: "; 
  run<Pred_using_kernel>(p,q,r,points,std::back_inserter(res0));
  std::cout << "Running static optimized predicates: ";
  run<Pred_using_optimized_static>(p,q,r,points,std::back_inserter(res1));
  std::cout << "Running predicates with intervals: ";
  run<Pred_using_intervals>(p,q,r,points,std::back_inserter(res2));
  
  for (std::size_t k=0;k<nb_pts;++k)
    if(res0[k]!=res1[k] || res1[k]!=res2[k]){
      std::cerr << "ERROR results different\n";
      exit(EXIT_FAILURE);
    }
      
  
  return 0;
}

