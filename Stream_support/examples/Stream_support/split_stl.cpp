#include <CGAL/IO/STL.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>

#include <vector>
#include <array>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef K::Point_3                                            Point;


struct BB {
  CGAL::Bbox_3 bbox;
  BB()
    : bbox()
  {}

  void push_back(const Point& p)
  {
    bbox += p.bbox();
  }
};

template <typename T>
struct Empty_set {
  typedef typename std::vector<T>::const_iterator const_iterator;
  int m_size = 0;

  void push_back(const T& )
  {
    ++m_size;
  }

  int size() const
  {
    return m_size;
  }
};

namespace boost {
  template <>
  struct range_value<BB> {
    typedef Point type;
  };

  template <typename T>
  struct range_value<Empty_set<T>> {
    typedef T type;
  };
}

int main(int argc, char* argv[])
{
    CGAL::Timer t;
    t.start();
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/pig.stl");
  const std::string stem = (argc > 2) ? argv[2] : "split";
  int gridlength = (argc > 3) ? std::stoi(argv[3]) : 3;

  Empty_set<std::array<int,3>> count_faces;
  BB bbox_points;

  {
    std::ifstream in(filename, std::ios::binary);
    CGAL::IO::read_STL(in,
                       bbox_points,
                       count_faces,
                       CGAL::parameters::verbose(true).unique(false));

    std::cout << "#triangles = " << count_faces.size() << std::endl;
    std::cout << t.time() << " sec." << std::endl;
  }

  {
    std::ifstream in(filename, std::ios::binary);
    CGAL::IO::split_binary_STL(in, stem, bbox_points.bbox, gridlength, true);

  }
  return 0;
}
