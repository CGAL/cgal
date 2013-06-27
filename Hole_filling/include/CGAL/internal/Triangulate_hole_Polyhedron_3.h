#ifndef CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H
#define CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H

#include <limits>
#include <map>
#include <set>
#include <vector>

#include <CGAL/assertions.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>

namespace CGAL {
namespace internal {

template<class Polyhedron>
class Triangulate_hole_Polyhedron_3{
// typedefs
  typedef typename Polyhedron::Traits::Point_3 Point_3;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
  typedef typename Polyhedron::Facet_handle Facet_handle;
  typedef typename Polyhedron::Halfedge_around_facet_circulator  Halfedge_around_facet_circulator;
  typedef std::vector<Halfedge_handle> Polyline_3;

  struct Weight {

    std::pair<double,double> w;

    Weight()
      : w(std::make_pair<double>(0, 0))
    {}

    Weight(double angle, double area)
      : w(angle, area)
    {}

    friend Weight operator+(const Weight& w1, const Weight& w2)
    {
      return Weight((std::max)(w1.w.first, w2.w.first), w1.w.second + w2.w.second);
    }

    friend bool operator<(const Weight& w1, const Weight& w2)
    {
      if(w1.w.first < w2.w.first){
        return true;
      }else if(w1.w.first == w2.w.first){
        return w1.w.second < w2.w.second;
      }
      return false;
    }
  };

  Weight
  sigma(const Polyline_3& P, int i, int j, int k)
  {
    assert(i+1 ==j);
    assert(j+1 ==k);
    const Point_3& p = P[i]->opposite()->vertex()->point();
    const Point_3& q = P[j]->opposite()->vertex()->point();
    const Point_3& r = P[k]->opposite()->vertex()->point();
    // The CGAL::dihedral angle is measured between the oriented triangles, that is it goes from [-pi, pi]
    // What we need is the angle between the normals of the triangles between [0, pi]
    double ang1 = 180 - CGAL::abs( CGAL::Mesh_3::dihedral_angle(p,q,r,P[i]->opposite()->next()->vertex()->point()));
    double ang2 = 180 - CGAL::abs( CGAL::Mesh_3::dihedral_angle(q,r,p, P[j]->opposite()->next()->vertex()->point()));
    return Weight((std::max)(ang1, ang2), std::sqrt(CGAL::squared_area(p,q,r)));
  }

  Weight
  sigma(const Polyline_3& P, int i, int m, int k, const std::vector<int>& lambda)
  {
    assert(i < m);
    assert(m < k);
    int n = P.size() -1; // because the first and last point are equal
    const Point_3& pi = P[i]->opposite()->vertex()->point();
    const Point_3& pm = P[m]->opposite()->vertex()->point();
    const Point_3& pk = P[k]->opposite()->vertex()->point();
    double ang1=0, ang2=0;
    if(lambda[i*n+m] != -1){
      const Point_3& pim = P[lambda[i*n+m]]->opposite()->vertex()->point();
      ang1 = 180 - CGAL::abs(CGAL::Mesh_3::dihedral_angle(pi,pm,pk, pim));
    }
    if(lambda[m*n+k] != -1){
      const Point_3& pmk = P[lambda[m*n+k]]->opposite()->vertex()->point();
      ang2 = 180 - CGAL::abs( CGAL::Mesh_3::dihedral_angle(pm,pk,pi, pmk));
    }
    return Weight((std::max)(ang1, ang2), std::sqrt(CGAL::squared_area(pi,pm,pk)));
  }

  template<class OutputIterator>
  Halfedge_handle add_facets(const Polyline_3& P, const std::vector<int>& lambda, 
    int i, int k, Polyhedron& poly, OutputIterator& out, bool last = true)
  {
    if(i + 1 == k) { return P[i]; }

    Halfedge_handle h, g;
    int n = P.size() -1; // because the first and last point are equal
    if(i+2 == k){
      if(last)
      { h = poly.fill_hole(P[i+1]); }
      else 
      { h = poly.add_facet_to_border(P[i]->prev(), P[i+1]); }
      
      assert(h->facet() != Facet_handle());
      *out++ = h->facet();
      return h->opposite();
    } 
    else 
    {
      int la = lambda[i*n + k];
      h = add_facets(P, lambda, i, la, poly, out, false);
      g = add_facets(P, lambda, la, k, poly, out, false);

      if(last)
      { h = poly.fill_hole(g); }
      else 
      { h = poly.add_facet_to_border(h->prev(), g); }

      assert(h->facet() != Facet_handle());
      *out++ = h->facet();
      return h->opposite();
    }
  }

  void compute_lambda(const Polyline_3& P, std::vector<int>& lambda) {
    int n = P.size() - 1; // because the first and last point are equal
    lambda = std::vector<int>(n*n,-1);
    std::vector<Weight> W(n*n,Weight(0,0));

    for(int i=0; i < n-2; ++i){
      W[i*n + (i+2)] = sigma(P, i, i+1, i+2);
      lambda[i*n + (i+2)] = i+1;
    }

    for(int j = 3; j< n; ++j){
      for(int i=0; i<n-j; ++i){
        int k = i+j;
        int m_min = 0;
        Weight w_min((std::numeric_limits<double>::max)(), (std::numeric_limits<double>::max)());
        for(int m = i+1; m<k; ++m){
          Weight w = W[i*n + m] + W[m*n + k] + sigma(P,i,m,k, lambda);
          if(w < w_min){
            w_min = w;
            m_min = m;
          }
        }
        W[i*n+k] = w_min;
        lambda[i*n+k] = m_min;
      }
    }
  }

public:
  template<class OutputIterator>
  void operator()(Polyhedron& poly, Halfedge_handle it, OutputIterator& out) {
    Polyline_3 P;
    Halfedge_around_facet_circulator circ(it), done(circ);
    do{ 
      P.push_back(circ);
    } while (++circ != done);
    P.push_back(circ);

    std::vector<int> lambda;
    compute_lambda(P, lambda);

    int n = P.size() - 1; // because the first and last point are equal
    add_facets(P, lambda, 0, n-1, poly, out);
  }
};

}//namespace internal

/**
 * @brief Function triangulating a hole in surface mesh.
 * 
 * @tparam Polyhedron a %CGAL polyhedron
 * @tparam OutputIterator iterator holding 'Polyhedron::Facet_handle' for patch facets.
 *
 * @param[in, out] polyhedron surface mesh which has the hole
 * @param border_halfedge a border halfedge incident to the hole
 * @param[out] output iterator over patch facets.
 * 
 */
template<class Polyhedron, class OutputIterator>
OutputIterator 
triangulate_hole(Polyhedron& polyhedron, 
  typename Polyhedron::Halfedge_handle border_halfedge, 
  OutputIterator output
  )
{
  internal::Triangulate_hole_Polyhedron_3<Polyhedron>()(polyhedron, border_halfedge, output);
  return output;
}

}//namespace CGAL
#endif //CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H