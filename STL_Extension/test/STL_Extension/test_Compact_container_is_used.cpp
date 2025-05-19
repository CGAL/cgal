// test program for Compact_container.

#include <CGAL/Compact_container.h>

class Node_1
{
  union {
    Node_1 * p;
    void * p_cc;
  };

public:

  Node_1() : p(nullptr)
  {}

  void *   for_compact_container() const { return p_cc; }
  void for_compact_container(void *p) { p_cc = p; }
};

int main()
{
  typedef CGAL::Compact_container<Node_1> C1;

  C1 c1;
  if (!c1.empty())
  {
    std::cout<<"PB new container is not empty."<<std::endl;
    return EXIT_FAILURE;
  }

  std::size_t nb=0;
  for (nb = 0 ; nb < 10000 ; ++nb)
  {
    C1::iterator it=c1.emplace();
    if ( !c1.is_used(it) )
    {
      std::cout<<"PB new emplace element is not used."<<std::endl;
      return EXIT_FAILURE;
    }
    if ( !c1.is_used(nb) )
    {
      std::cout<<"PB new emplace element is not used (2)."<<std::endl;
      return EXIT_FAILURE;
    }
    if ( c1.index(it)!=nb )
    {
      std::cout<<"PB index of new emplace element is not correct."<<std::endl;
      return EXIT_FAILURE;
    }
    if ( &(c1[nb])!=&*it )
    {
      std::cout<<"PB operator[] gives a wrong result."<<std::endl;
      return EXIT_FAILURE;
    }
  }

  nb=0;
  for (C1::iterator it = c1.begin(), itend=c1.end(); it!=itend; ++it, ++nb)
  {
    c1.erase(it);
    if ( c1.is_used(it) )
    {
      std::cout<<"PB erase element is used."<<std::endl;
      return EXIT_FAILURE;
    }
    if ( c1.is_used(nb) )
    {
      std::cout<<"PB erase element is used (2)."<<std::endl;
      return EXIT_FAILURE;
    }
  }

  if (!c1.empty())
  {
    std::cout<<"PB container at the end is not empty."<<std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
