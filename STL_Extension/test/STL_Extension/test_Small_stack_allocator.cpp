#include <CGAL/Small_stack_allocator.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Random.h>

#include <list>

constexpr std::size_t nb_tests = 10000000;

typedef std::list<int> Standard_list;

constexpr std::size_t max_size = 4;
typedef CGAL::Small_stack_allocator<int, max_size> SSA;
typedef std::list<int, SSA> List_with_small_stack;

template <typename ListType>
double speed_test (std::size_t nb_min, std::size_t nb_max)
{
  CGAL::Real_timer t;
  t.start();
  CGAL::Random rand(0);
  for (std::size_t i = 0; i < nb_tests; ++ i)
  {
    ListType list;
    std::size_t nb = rand.get_int(nb_min, nb_max);
    for (std::size_t j = 0; j < nb; ++ j)
      list.push_back (rand.get_int(0, 1000000));
  }
  t.stop();

  return t.time();
}

int main()
{
  std::cerr << "[Speed test with all insertions smaller than max_size]" << std::endl
            << " * Standard list:   "
            << speed_test<Standard_list>(max_size - 2, max_size)
            << "s" << std::endl
            << " * List with stack: "
            << speed_test<List_with_small_stack>(max_size - 2, max_size)
            << "s" << std::endl
            << "   (should be much faster)" << std::endl;

  std::cerr << "[Speed test with ~50% insertions smaller than max_size]" << std::endl
            << " * Standard list:   "
            << speed_test<Standard_list>(max_size - 2, max_size + 2)
            << "s" << std::endl
            << " * List with stack: "
            << speed_test<List_with_small_stack>(max_size - 2, max_size + 2)
            << "s" << std::endl
            << "   (should still be slightly faster)" << std::endl;

  std::cerr << "[Speed test with all insertions greater than max_size]" << std::endl
            << " * Standard list:   "
            << speed_test<Standard_list>(max_size + 2, max_size + 4)
            << "s" << std::endl
            << " * List with stack: "
            << speed_test<List_with_small_stack>(max_size + 2, max_size + 4)
            << "s" << std::endl
            << "   (shouldn't be faster)" << std::endl;

  return EXIT_SUCCESS;
}
