#include <CGAL/Handle_for_virtual.h>

#ifndef CGAL_LINKED_WITH_TBB

#include <iostream>

int main()
{
  std::cout << "NOTICE: this test needs CGAL_LINKED_WITH_TBB, and will not be tested."
            << std::endl;
  return 0;
}

#else

#include <atomic>
#include <cstddef>
#include <tbb/parallel_for.h>

class Test_rep : public CGAL::Ref_counted_virtual
{
public:
  ~Test_rep() override { ++destructions; }

  static std::atomic<unsigned int> destructions;
};

std::atomic<unsigned int> Test_rep::destructions{0};

int main()
{
  using Handle = CGAL::Handle_for_virtual<Test_rep>;

  {
    Handle handle(Test_rep{});
    Test_rep::destructions.store(0, std::memory_order_relaxed);

    constexpr std::size_t copies = 800000;
    tbb::parallel_for(std::size_t(0), copies, [&handle](std::size_t)
    {
      Handle copy(handle);
    });
  }

  return Test_rep::destructions.load(std::memory_order_relaxed) == 1 ? 0 : 1;
}

#endif
