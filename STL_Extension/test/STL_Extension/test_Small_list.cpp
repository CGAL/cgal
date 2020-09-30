#include <CGAL/Small_list.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Random.h>

#include <list>

constexpr std::size_t nb_tests = 1000000;
constexpr std::size_t max_size = 10;

typedef std::list<int> Standard_list;
typedef CGAL::Small_list<int, max_size> Small_list;

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
  std::cerr << "[Small list sanity check]" << std::endl;
  {
    Small_list list;
    CGAL_assertion (list.size() == 0);
    CGAL_assertion (list.empty());

    std::cerr << "push_back() 0 1 2 3 4" << std::endl;
    list.push_back (0);
    list.push_back (1);
    list.push_back (2);
    list.push_back (3);
    list.push_back (4);
    CGAL_assertion (list.size() == 5);

    std::cerr << "begin() / end()" << std::endl;
    for (const int& i : list)
      std::cerr << i << " ";
    std::cerr << std::endl;
    std::cerr << "rbegin() / rend()" << std::endl;
    for (auto it = list.rbegin(); it != list.rend(); ++ it)
      std::cerr << *it << " ";
    std::cerr << std::endl;

    auto it1 = list.begin();
    auto it2 = list.end();
    ++ it1;
    -- it2;
    std::cerr << "erase(1, 4)" << std::endl;
    list.erase(it1, it2);
    CGAL_assertion (list.size() == 2);

    std::cerr << "begin() / end()" << std::endl;
    for (const int& i : list)
      std::cerr << i << " ";
    std::cerr << std::endl;
    std::cerr << "rbegin() / rend()" << std::endl;
    for (auto it = list.rbegin(); it != list.rend(); ++ it)
      std::cerr << *it << " ";
    std::cerr << std::endl;

    std::cerr << "insert(begin(), 6)" << std::endl;
    list.insert(list.begin(), 6);
    std::cerr << "begin() / end()" << std::endl;
    for (const int& i : list)
      std::cerr << i << " ";
    std::cerr << std::endl;
    std::cerr << "rbegin() / rend()" << std::endl;
    for (auto it = list.rbegin(); it != list.rend(); ++ it)
      std::cerr << *it << " ";
    std::cerr << std::endl;

    std::cerr << "erase(begin(), end())" << std::endl;
    list.erase(list.begin(), list.end());
    CGAL_assertion (list.size() == 0);
    CGAL_assertion (list.empty());
  }

  std::cerr << "[Speed test with all insertions smaller than max_size]" << std::endl
            << " * Standard list:   "
            << speed_test<Standard_list>(max_size - 2, max_size)
            << "s" << std::endl
            << " * List with stack: "
            << speed_test<Small_list>(max_size - 2, max_size)
            << "s" << std::endl
            << "   (should be much faster)" << std::endl;

  std::cerr << "[Speed test with ~50% insertions smaller than max_size]" << std::endl
            << " * Standard list:   "
            << speed_test<Standard_list>(max_size - 2, max_size + 2)
            << "s" << std::endl
            << " * List with stack: "
            << speed_test<Small_list>(max_size - 2, max_size + 2)
            << "s" << std::endl
            << "   (should still be slightly faster)" << std::endl;

  std::cerr << "[Speed test with all insertions greater than max_size]" << std::endl
            << " * Standard list:   "
            << speed_test<Standard_list>(max_size + 2, max_size + 4)
            << "s" << std::endl
            << " * List with stack: "
            << speed_test<Small_list>(max_size + 2, max_size + 4)
            << "s" << std::endl
            << "   (shouldn't be faster)" << std::endl;

  return EXIT_SUCCESS;
}
