#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include <vector>
#include <limits>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef std::vector<Point_3> Polyline_3;

double
sigma(const Polyline_3& P, int i, int j, int k)
{
  return sqrt(CGAL::squared_area(P[i], P[j], P[k]));
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
    std::cout << P[i] << std::endl;
  }
  for(int i = 0; i < facets.size();){
    std::cout << "3 " << facets[i] << " " << facets[i+1] << " " << facets[i+2] << std::endl;
    i+=3;
  }
}


void fill(const Polyline_3& P)
{
  //  std::cerr << "n = " << P.size() << std::endl;
  int n = P.size() -1; // because the first and last point are equal
  std::vector<double> W(n*n,0);
  std::vector<int> lambda(n*n,-1);

  for(int i=0; i < n-2; ++i){
      W[i*n + (i+2)] = sigma(P, i, i+1, i+2);
  }

  for(int j = 3; j< n; ++j){ // in the papers it is j < n-1
    for(int i=0; i<n-j; ++i){
      int k = i+j;
      //      std::cerr << "i = " << i <<  "  j = " << j << "  k = " << k << std::endl;
      int m_min = 0;
      double w_min = (std::numeric_limits<double>::max)();
      for(int m = i+1; m<k; ++m){
        double w = W[i*n + m] + W[m*n + k] + sigma(P,i,m,k);
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
  int n;
  Polyline_3 P;
  Point_3 p;
  std::cin >> n;
  P.reserve(n);
  while(n--){
    std::cin >> p;
    P.push_back(p);
  }
  fill(P);

  return 0;
}
