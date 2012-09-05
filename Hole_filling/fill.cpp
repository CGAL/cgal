#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
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
  return Weight(0, sqrt(CGAL::squared_area(P[i]->vertex()->point(), P[j]->vertex()->point(), P[k]->vertex()->point())));
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
    facets.push_back(i%n);
    facets.push_back(la%n);
    facets.push_back(k%n);
    if(la != k-1){
      trace(P, lambda, la, k, facets);
    }
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
    std::cout << P[i]->vertex()->point() << std::endl;
  }
  for(int i = 0; i < facets.size();){
    std::cout << "3 " << facets[i] << " " << facets[i+1] << " " << facets[i+2] << std::endl;
    i+=3;
  }
}


void fill(const Polyline_3& P)
{
  //  std::cerr << "n = " << P.size() << std::endl;
  int n = P.size() - 1; // because the first and last point are equal
  std::vector<Weight> W(n*n,Weight(0,0));
  std::vector<int> lambda(n*n,-1);

  for(int i=0; i < n-2; ++i){
      W[i*n + (i+2)] = sigma(P, i, i+1, i+2);
  }

  for(int j = 3; j< n; ++j){ // in the papers it is j < n-1
    for(int i=0; i<n-j; ++i){
      int k = i+j;
      //      std::cerr << "i = " << i <<  "  j = " << j << "  k = " << k << std::endl;
      int m_min = 0;
      Weight w_min((std::numeric_limits<double>::max)(), (std::numeric_limits<double>::max)());
      for(int m = i+1; m<k; ++m){
        Weight w = W[i*n + m] + W[m*n + k] + sigma(P,i,m,k);
        if(w < w_min){
          w_min = w;
          m_min = m;
          //          std::cerr << "w_min = " << w_min;
          //          std::cerr << "   m_min = " << m_min << std::endl;
        }
      }
      W[i*n+k] = w_min;
      lambda[i*n+k] = m_min;
    }
  }
  //  std::cerr << "weight: " <<  W[n-1] << std::endl;
  trace(P, lambda);
}



int main()
{
  Polyline_3 P;
  Polyhedron poly;
  std::cin >> poly;

  for(Halfedge_iterator it = poly.halfedges_begin(); it != poly.halfedges_end(); ++it){
    if(it->is_border()){
      Halfedge_around_facet_circulator circ(it), done(circ);
      do{
        P.push_back(circ);
      }while (++circ != done);
      P.push_back(circ);
      break;
    }
  }
  fill(P);

  return 0;
}
