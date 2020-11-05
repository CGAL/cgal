#include <CGAL/Simple_cartesian.h>
#include <CGAL/box_intersection_d.h>
#include <vector>
#include <fstream>

typedef CGAL::Simple_cartesian<float>             Kernel;
typedef Kernel::Point_3                           Point_3;

std::vector<Point_3>  points;
std::vector<Point_3*> boxes;     // boxes are just pointers to points
const float           eps = 0.1f; // finds point pairs of distance < 2*eps

// Boxes are just pointers to 3d points. The traits class adds the
// +- eps size to each interval around the point, effectively building
// on the fly a box of size 2*eps centered at the point.
struct Traits {
    typedef float          NT;
    typedef Point_3*       Box_parameter;
    typedef std::ptrdiff_t ID;

    static int   dimension() { return 3; }
    static float coord( Box_parameter b, int d) {
        return (d == 0) ? b->x() : ((d == 1) ? b->y() : b->z());
    }
    static float min_coord( Box_parameter b, int d) { return coord(b,d)-eps;}
    static float max_coord( Box_parameter b, int d) { return coord(b,d)+eps;}
    // id-function using address of current box,
    // requires to work with pointers to boxes later
    static std::ptrdiff_t id(Box_parameter b) { return (std::ptrdiff_t)(b); }
};

// callback function reports pairs in close proximity
void report( const Point_3* a, const Point_3* b) {
    float dist = std::sqrt( CGAL::squared_distance( *a, *b));
    if ( dist < 2*eps) {
        std::cout << "Point " << (a - &(points.front())) << " and Point "
                  << (b - &(points.front())) << " have distance " << dist
                  << "." << std::endl;
    }
}

int main(int argc, char*argv[]) {

  std::ifstream in((argc>1)?argv[1]:"data/points.xyz");
  Point_3 p;
  while(in >> p){
    points.push_back(p);
  }
  for(std::size_t i = 0; i< points.size();++i) {
    boxes.push_back( &points[i]);
  }
  // run the intersection algorithm and report proximity pairs
  CGAL::box_self_intersection_d( boxes.begin(), boxes.end(),
                                 report, Traits());
  return 0;
}
