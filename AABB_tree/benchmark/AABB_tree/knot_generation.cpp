
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/polygon_soup_io.h>

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;
using Triangle = std::array<int, 3>;

int main(int argc, char* argv[])
{
  const std::string filename = argc > 1 ? argv[1] : "knot.off";

  const int p = argc > 2 ? std::atoi(argv[2]) : 1; // Nb of major turns
  const int q = argc > 3 ? std::atoi(argv[3]) : 6; // Nb of knot turns

  const double R = argc > 4 ? std::atof(argv[4]) : 2.0; // major radius
  const double r = argc > 5 ? std::atof(argv[5]) : 0.9; // knot radius
  const double tube_radius = argc > 6 ? std::atof(argv[6]) : 0.4;

  const int tubular_segments = argc > 7 ? std::atoi(argv[7]) : 120;
  const int radial_segments  = argc > 8 ? std::atoi(argv[8]) : 60;

  std::vector<Point> vertices;
  std::vector<Triangle> faces;

  // Generate knot polyline
  auto knot = [&](double t){
    double cq = std::cos(q * t);
    double sq = std::sin(q * t);
    double factor = R + r * cq;
    return Point(factor * std::cos(p * t),
                 factor * std::sin(p * t),
                 r * sq);
  };
  auto knot_derivative = [&](double t){
    double cq = std::cos(q * t);
    double sq = std::sin(q * t);
    double cp = std::cos(p * t);
    double sp = std::sin(p * t);

    double factor = R + r*cq;
    double dfactor = -r*q*sq;

    return Vector(dfactor*cp - factor*p*sp,
                  dfactor*sp + factor*p*cp,
                  r*q*cq);
  };

  const double dt = 2.0 * M_PI / tubular_segments;
  const double dr = 2.0 * M_PI / radial_segments;

  // A normal to the knot axis
  Vector N(1,0,0);
  for(int i=0; i<tubular_segments; ++i)
  {
    double t = i * dt;

    // knot point and its tangent
    Point C  = knot(t);
    Vector T = knot_derivative(t);
    // Compute bases of the plane on C
    Vector B =  CGAL::cross_product(N, T);
    N = -CGAL::cross_product(B, T);
    N /= CGAL::approximate_sqrt(N.squared_length());
    B /= CGAL::approximate_sqrt(B.squared_length());

    for(int j=0; j<radial_segments; ++j){
      double phi = j * dr;
      vertices.push_back(C + tube_radius*(std::cos(phi)*N + std::sin(phi)*B));
    }
  }

  auto idx = [&](int i, int j){
    i %= tubular_segments;
    j %= radial_segments;
    return i * radial_segments + j;
  };

  for(int i=0; i<tubular_segments-1; ++i)
    for(int j=0; j<radial_segments; ++j){
      int v00 = idx(i  , j  );
      int v10 = idx(i+1, j  );
      int v01 = idx(i  , j+1);
      int v11 = idx(i+1, j+1);

      faces.push_back({v00,v10,v11});
      faces.push_back({v00,v11,v01});
    }

  // To avoid a twist, we compute the offset that minimize the distance between a point at tubular index N-1 and tubular index 0
  int n = tubular_segments-1;
  int offset = 0;
  double min = CGAL::squared_distance(vertices[idx(n, 0)], vertices[idx(0,0)]);
  for(int j=1; j<radial_segments; ++j){
    double sq = CGAL::squared_distance(vertices[idx(n, 0)], vertices[idx(0,j)]);
    if(min >sq){
      min = sq;
      offset = j;
    }
  }
  for(int j=0; j<radial_segments; ++j){
    int v00 = idx(n, j  );
    int v10 = idx(0, j+offset);
    int v01 = idx(n, j+1);
    int v11 = idx(0, j+offset+1);

    faces.push_back({v00,v10,v11});
    faces.push_back({v00,v11,v01});
  }

  CGAL::IO::write_polygon_soup(filename, vertices, faces);

  return 0;
}