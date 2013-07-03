#ifndef CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYLINE_H
#define CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYLINE_H

#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <CGAL/utility.h>
#include <vector>
#include <limits>
#include <CGAL/value_type_traits.h>

namespace CGAL {
namespace internal {

template <typename K>
class Triangulate_hole_polyline {
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
    CGAL_assertion(i+1 ==j);
    CGAL_assertion(j+1 ==k);
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
    return Weight((std::max)(ang1, ang2), std::sqrt(CGAL::squared_area(p,q,r)));
  }
  
  
  Weight
  sigma(const Polyline_3& P, const Polyline_3& Q, int i, int m, int k, const std::vector<int>& lambda)
  {
    CGAL_assertion(i < m);
    CGAL_assertion(m < k);
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
    return Weight((std::max)(ang1, ang2), std::sqrt(CGAL::squared_area(pi,pm,pk)));
  }
  
  template <typename OutputIteratorValueType, typename OutputIterator>
  OutputIterator
  trace(int n,
        const std::vector<int>& lambda, 
        int i, 
        int k, 
        OutputIterator out)
  {
    if(i + 1 == k) { return out; }

    if(i+2 == k){
      *out++ = OutputIteratorValueType(i%n, (i+1)%n, k%n);
    } 
    else {
      int la = lambda[i*n + k];
      out = trace<OutputIteratorValueType>(n, lambda, i, la, out);
      *out++ = OutputIteratorValueType(i%n, la%n, k%n);
      out = trace<OutputIteratorValueType>(n, lambda, la, k, out);
    }
    return out;
  }
  
public:

  template <typename OutputIteratorValueType, typename OutputIterator>
  OutputIterator 
  triangulate(const Polyline_3& P, 
             const Polyline_3& Q,
             OutputIterator out)
  {
    CGAL_assertion(P.front() == P.back());
    CGAL_assertion(Q.empty() || (Q.front() == Q.back()));
    CGAL_assertion(Q.empty() || (P.size() == Q.size()));
    
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
    return  trace<OutputIteratorValueType>(n, lambda, 0, n-1, out);
    
  }

};
} // namespace internal

/*!
Creates triangles to fill the hole defined by points in the range `(pbegin,pend)`.
The range `(qbegin,qend)` indicate for each pair of consecutive points in the aforementioned range,
the third point of the facet this segment is incident to. Triangles are put into `out`
using the indices of the input points in the range `(pbegin,pend)`.
@tparam OutputIteratorValueType value type of OutputIterator having a constructor `OutputIteratorValueType(int p0, int p1, int p2)` available. 
        It is default to value_type_traits<OutputIterator>::type, and can be omitted when the default is fine
@tparam InputIterator iterator over input points
@tparam OutputIterator iterator over patch triangles
@param pbegin first iterator of the range of points
@param pend past-the-end iterator of the range of points
@param qbegin first iterator of the range of third points
@param qend past-the-end iterator of the range of third points
@param out iterator over output patch triangles
*/
template <typename OutputIteratorValueType, typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
          InputIterator qbegin, InputIterator qend, 
          OutputIterator out)
{
  typedef typename CGAL::Kernel_traits< typename std::iterator_traits<InputIterator>::value_type>::Kernel Kernel;
  typedef CGAL::internal::Triangulate_hole_polyline<Kernel> Fill;
  typename Fill::Polyline_3 P(pbegin, pend);
  typename Fill::Polyline_3 Q(qbegin, qend);
  if(P.front() != P.back()){
    P.push_back(P.front());
  }
  if(! Q.empty() && (Q.front() != Q.back())){
    Q.push_back(Q.front());
  }
  Fill fill;
  return fill.template triangulate<OutputIteratorValueType>(P,Q,out);
}

// overload for OutputIteratorValueType
template <typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
          InputIterator qbegin, InputIterator qend, 
          OutputIterator out)
{
  return triangulate_hole_polyline<typename value_type_traits<OutputIterator>::type>
    (pbegin, pend, qbegin, qend, out);
}


/*!
Creates triangles to fill the hole defined by points in the range `(pbegin,pend)`.
Triangles are put into `out` using the indices of the input points in the range `(pbegin,pend)`.
@tparam OutputIteratorValueType value type of OutputIterator having a constructor `OutputIteratorValueType(int p0, int p1, int p2)` available. 
        It is default to value_type_traits<OutputIterator>::type, and can be omitted when the default is fine
@tparam InputIterator iterator over input points
@tparam OutputIterator iterator over output patch triangles
@param pbegin first iterator of the range of points
@param pend past-the-end iterator of the range of points
@param out iterator over output patch triangles
*/
template <typename OutputIteratorValueType, typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
          OutputIterator out)
{
  typedef typename CGAL::Kernel_traits< typename std::iterator_traits<InputIterator>::value_type>::Kernel Kernel;
  typedef CGAL::internal::Triangulate_hole_polyline<Kernel> Fill;
  typename Fill::Polyline_3 P(pbegin, pend);
  typename Fill::Polyline_3 Q;
  if(P.front() != P.back()){
    P.push_back(P.front());
  }
  Fill fill;
  return fill.template triangulate<OutputIteratorValueType>(P,Q,out);
}

// overload for OutputIteratorValueType
template <typename InputIterator, typename OutputIterator>
OutputIterator
triangulate_hole_polyline(InputIterator pbegin, InputIterator pend, 
          OutputIterator out)
{
  return triangulate_hole_polyline<typename value_type_traits<OutputIterator>::type>
    (pbegin, pend, out);
}

} // namespace CGAL

#endif // CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYLINE_H
