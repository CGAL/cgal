#include <CGAL/Simple_cartesian.h>
#include <CGAL/Shape_regularization.h>

using Kernel = CGAL::Simple_cartesian<double>;
using FT     = typename Kernel::FT;

struct Custom_object {
  const std::string name;
  Custom_object(const std::string name_) :
  name(name_) { }
  // define your object here
};

struct Custom_neighbor_query_2 {
  void operator()(
    const std::size_t query_index, std::vector<std::size_t>& neighbors) {
    neighbors.clear();
    if (query_index == 0) { neighbors.push_back(1); } // first  object
    if (query_index == 1) { neighbors.push_back(0); } // second object
  }
};

struct Custom_regularization_2 {
  FT bound(const std::size_t ) const {
    return FT(5); // maximum angle change
  }
  FT target(const std::size_t , const std::size_t ) {
    return FT(0); // 0 angle change
  }
  void update(const std::vector<FT>& ) {
    // skip update
  }
};

template<typename NT>
class Custom_quadratic_program_traits  {
public:
  void  set_P(const std::size_t, const std::size_t, const FT) { }
  void  set_q(const std::size_t, const FT) { }
  void  set_r(const FT) { }
  void  set_A(const std::size_t, const std::size_t, const FT) { }
  void  set_l(const std::size_t, const FT) { }
  void  set_u(const std::size_t, const FT) { }
  void resize(const std::size_t, const std::size_t) { }

  template<typename OutputIterator>
  bool solve(OutputIterator solution) {

    // 3 = 2 objects + 1 edge between them
    for (std::size_t i = 0; i < 3; ++i) {
      *(++solution) = NT(0);
    }
    return true;
  }
};

using Objects = std::vector<Custom_object>;
using Neighbor_query = Custom_neighbor_query_2;
using Regularization_type = Custom_regularization_2;
using Quadratic_program = Custom_quadratic_program_traits<FT>;
using Regularizer = CGAL::Shape_regularization::QP_regularization<
  Kernel, Objects, Neighbor_query, Regularization_type, Quadratic_program>;

int main() {

  Neighbor_query neighbor_query;
  Regularization_type regularization_type;
  Quadratic_program quadratic_program;

  std::vector<Custom_object> objects = {
    Custom_object("first"), Custom_object("second")
  };

  Regularizer regularizer(
    objects, neighbor_query, regularization_type, quadratic_program);
  regularizer.regularize();

  std::cout << "* regularized 2 objects" << std::endl;
}
