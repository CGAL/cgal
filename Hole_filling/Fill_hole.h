
#ifndef CGAL_FILL_HOLE_H
#define CGAL_FILL_HOLE_H

#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <CGAL/utility.h>
#include <vector>
#include <limits>

namespace CGAL {

template <typename K>
class Fill_hole {
public:
  typedef typename K::Point_3 Point_3;
  typedef std::vector<Point_3> Polyline_3;

private:
  struct Weight {
    
    std::pair<double,double> w;
    
    Weight()
      : w(std::make_pair<double>(0, 0))
    {}
    
    Weight(double angle, double area)
      : w(angle, area)
    {}

  
  Weight operator+(const Weight& w2) const
  {
    return Weight((std::max)(w.first, w2.w.first), w.second + w2.w.second);
  }
  
  bool operator<(const Weight& w2) const
  {
    if(w.first < w2.w.first){
      return true;
    }else if(w.first == w2.w.first){
      return w.second < w2.w.second;
    }
    return false;
  }
  };


  Weight
  sigma(const Polyline_3& P, const Polyline_3& Q, int i, int j, int k)
  {
    assert(i+1 ==j);
    assert(j+1 ==k);
    const Point_3& p = P[i];
    const Point_3& q = P[j];
    const Point_3& r = P[k];
    // The CGAL::dihedral angle is measured between the oriented triangles, that is it goes from [-pi, pi]
    // What we need is the angle between the normals of the triangles between [0, pi]
    double ang1 = 0, ang2 = 0;
    
    if(!Q.empty()){
      ang1 = 180 - CGAL::abs( CGAL::Mesh_3::dihedral_angle(p,q,r,Q[i]));
      ang2 = 180 - CGAL::abs( CGAL::Mesh_3::dihedral_angle(q,r,p, Q[j]));
    }
    return Weight((std::max)(ang1, ang2), sqrt(CGAL::squared_area(p,q,r)));
  }
  
  
  Weight
  sigma(const Polyline_3& P, const Polyline_3& Q, int i, int m, int k, const std::vector<int>& lambda)
  {
    assert(i < m);
    assert(m < k);
    int n = P.size() -1; // because the first and last point are equal
    const Point_3& pi = P[i];
    const Point_3& pm = P[m];
    const Point_3& pk = P[k];
    double ang1=0, ang2=0;
    
    if(!Q.empty()){
      if(lambda[i*n+m] != -1){
        const Point_3& pim = P[lambda[i*n+m]];
        ang1 = 180 - CGAL::abs(CGAL::Mesh_3::dihedral_angle(pi,pm,pk, pim));
      }
      if(lambda[m*n+k] != -1){
        const Point_3& pmk = P[lambda[m*n+k]];
        ang2 = 180 - CGAL::abs( CGAL::Mesh_3::dihedral_angle(pm,pk,pi, pmk));
      }
    }
    return Weight((std::max)(ang1, ang2), sqrt(CGAL::squared_area(pi,pm,pk)));
  }
  
  template <typename OutputIterator>
  OutputIterator
  trace(int n,
        const std::vector<int>& lambda, 
        int i, 
        int k, 
        OutputIterator out)
  {            
    if(i+2 == k){
      *out++ = CGAL::Triple<int,int,int>(i%n, (i+1)%n, k%n);
    } else {
      int la = lambda[i*n + k];
      if(la != i+1){
        out = trace(n, lambda, i, la, out);
      }
      *out++ = CGAL::Triple<int,int,int>(i%n, la%n, k%n);
    if(la != k-1){
      out = trace(n, lambda, la, k, out);
    }
    }
    return out;
  }
  
public:

  template <typename OutputIterator>
  OutputIterator 
  operator()(const Polyline_3& P, 
             const Polyline_3& Q,
             OutputIterator out)
  {
    assert(P.front() == P.back());
    assert(Q.empty() || (Q.front() == Q.back()));
    assert(Q.empty() || (P.size() == Q.size()));
    
    int n = P.size() - 1; // because the first and last point are equal
    std::vector<Weight> W(n*n,Weight(0,0));
    std::vector<int> lambda(n*n,-1);
    
    for(int i=0; i < n-2; ++i){
      W[i*n + (i+2)] = sigma(P, Q, i, i+1, i+2);
      lambda[i*n + (i+2)] = i+1;
    }
    
    for(int j = 3; j< n; ++j){
      for(int i=0; i<n-j; ++i){
        int k = i+j;
        int m_min = 0;
        Weight w_min((std::numeric_limits<double>::max)(), 
                     (std::numeric_limits<double>::max)());
        for(int m = i+1; m<k; ++m){
          Weight w = W[i*n + m] + W[m*n + k] + sigma(P,Q,i,m,k, lambda);
          if(w < w_min){
          w_min = w;
          m_min = m;
          }
        }
        W[i*n+k] = w_min;
        lambda[i*n+k] = m_min;
      }
    }
    return  trace(n, lambda, 0, n-1, out);
    
  }

};


template <typename InputIterator, typename OutputIterator>
OutputIterator
fill_hole(InputIterator pbegin, InputIterator pend, 
          InputIterator qbegin, InputIterator qend, 
          OutputIterator out)
{
  typedef typename CGAL::Kernel_traits< typename std::iterator_traits<InputIterator>::value_type>::Kernel Kernel;
  typedef Fill_hole<Kernel> Fill;
  typename Fill::Polyline_3 P(pbegin, pend);
  typename Fill::Polyline_3 Q(qbegin, qend);
  if(P.front() != P.back()){
    P.push_back(P.front());
  }
  if(! Q.empty() && (Q.front() != Q.back())){
    Q.push_back(Q.front());
  }
  Fill fill;
  return fill(P,Q,out);
}

} // namespace CGAL

#endif // CGAL_FILL_HOLE_H
