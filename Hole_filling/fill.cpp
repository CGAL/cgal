#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef Polyhedron::Vertex_handle Vertex_handle;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator  Halfedge_around_facet_circulator;

typedef std::vector<Halfedge_handle> Polyline_3;


struct Weight {

  std::pair<double,double> w;

  Weight()
    : w(std::make_pair<double>(0, 0))
  {}

  Weight(double angle, double area)
    : w(angle, area)
  {}
};

Weight operator+(const Weight& w1, const Weight& w2)
{
  return Weight((std::max)(w1.w.first, w2.w.first), w1.w.second + w2.w.second);
}

bool operator<(const Weight& w1, const Weight& w2)
{
  if(w1.w.first < w2.w.first){
    return true;
  }else if(w1.w.first == w2.w.first){
    return w1.w.second < w2.w.second;
  }
  return false;
}



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
  return Weight((std::max)(ang1, ang2), sqrt(CGAL::squared_area(p,q,r)));
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
  return Weight((std::max)(ang1, ang2), sqrt(CGAL::squared_area(pi,pm,pk)));
}


void trace(const Polyline_3& P, const std::vector<int>& lambda, int i, int k, std::vector<int>& facets)
{
  int n = P.size() -1; // because the first and last point are equal
  if(i+2 == k){
    facets.push_back(i%n);
    facets.push_back((i+1)%n);
    facets.push_back(k%n);
  } else {
    int la = lambda[i*n + k];
    if(la != i+1){
      trace(P, lambda, i, la, facets);
    }
    if(la != k-1){
      trace(P, lambda, la, k, facets);
    }

    facets.push_back(i%n);
    facets.push_back(la%n);
    facets.push_back(k%n);
  }
}


void trace(const Polyline_3& P, const std::vector<int>& lambda)
{
  std::cout.precision(20);
  int n = P.size() -1; // because the first and last point are equal
  std::vector<int> facets;
  trace(P, lambda, 0, n-1, facets);

  std::cout << "OFF\n" << n << " " << facets.size()/3 << " 0" << std::endl;
  for (int i =0; i < n; i++){
    std::cout << P[i]->opposite()->vertex()->point() << std::endl;
  }
  for(int i = 0; i < facets.size();){
    std::cout << "3 " << facets[i] << " " << facets[i+1] << " " << facets[i+2] << std::endl;
    i+=3;
  }
}



Halfedge_handle trace(const Polyline_3& P, const std::vector<int>& lambda, int i, int k, Polyhedron& poly, bool last = true)
{

  Halfedge_handle h, g;
  int n = P.size() -1; // because the first and last point are equal
  if(i+2 == k){
    if(last){
      return poly.fill_hole(P[i+1]);
    } else {
      return poly.add_facet_to_border(P[i]->prev(), P[i+1])->opposite();
    }
  } else {
    int la = lambda[i*n + k];
    if(la != i+1){
      h = trace(P, lambda, i, la, poly, false);
    } else {
      h = P[i];
    }
    if(la != k-1){
      g = trace(P, lambda, la, k, poly, false);
    } else {
      g = P[la];
    }
    if(last){
      return poly.fill_hole(g);
    } else {
      return poly.add_facet_to_border(h->prev(), g)->opposite();
    }
  }
}




void fill(Polyhedron& poly, Halfedge_handle it)
{
  if(! it->is_border()){
    return;
  }
  Polyline_3 P;
  Halfedge_around_facet_circulator circ(it), done(circ);
  do{
    P.push_back(circ);
  } while (++circ != done);
  P.push_back(circ);

  int n = P.size() - 1; // because the first and last point are equal
  std::vector<Weight> W(n*n,Weight(0,0));
  std::vector<int> lambda(n*n,-1);

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
  trace(P, lambda, 0, n-1, poly);
}



int main()
{
  Polyhedron poly;
  std::cin >> poly;

  int H = 20;
  for(Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ++it){
    if(it->is_border()){
      fill(poly, it);
    }
  }
  std::cout.precision(20);
  std::cout << poly << std::endl;

  std::cerr << "done" << std::endl;
  return 0;
}
