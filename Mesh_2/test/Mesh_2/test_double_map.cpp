#include <CGAL/basic.h>
#include <CGAL/Double_map.h>
#include <utility>
#include <CGAL/Random.h>

#include <list>
#include <iostream>

#include <cstdlib>
#include <cassert>

typedef CGAL::Double_map<int, int> Map;

int main(int argc, char** argv)
{
  unsigned int number_of_elements = 500;

#ifdef CGAL_USE_BOOST_BIMAP
  std::cerr << "(Using the \"Boost.Bimap implementation\" of <CGAL/Double_map.h>...)\n\n";
#endif

  if(argc > 1)
    number_of_elements = std::atoi(argv[1]);

  /* FOR VISUAL DEBUGGING */
  std::cerr << "Creating empty maps f and f2...\n";
  Map f;
  Map f2;

  std::cerr << "Filling f with f(i)=2^i-1, i=8..3 ...\n";
  int power=256;
  for(int n=9; n>2; n--)
  {
    std::cerr << "insertion of (" << n << ", " << power-1 << ")...\n";
    f.insert(n, power-1);
    power /= 2;
  }
  std::cerr << "Display of f:\n";
  f.dump_direct_func(std::cerr);
  std::cerr << "Display of f^(-1):\n";
  f.dump_reverse_func(std::cerr);
  std::cerr << "erase(5)...\n";
  f.erase(5);
  std::cerr << "Display of f:\n";
  f.dump_direct_func(std::cerr);
  std::cerr << "Display of f^(-1):\n";
  f.dump_reverse_func(std::cerr);
  std::cerr << "Copy f into f2...\n";
  f2 = f;
  std::cerr << "Display of f2:\n";
  f2.dump_direct_func(std::cerr);
  std::cerr << "Display of f2^(-1):\n";
  f2.dump_reverse_func(std::cerr);

  std::cerr << "\nEmptying f and f2...\n\n";
  f.clear();
  f2.clear();

  assert(f.empty() && f2.empty());
  /* AUTOMATIC CHECKS */
  std::cerr << "Filling f with " << number_of_elements << " random integers...\n";
  for(unsigned int n=0; n<number_of_elements; n++)
  {
    int i=CGAL::get_default_random().get_int(0,1000);
    f.insert(i, i*i);
  }

  std::cerr << "Check f.size<=" << number_of_elements << ".\n";
  assert(f.size()<=number_of_elements);

  std::cerr << "Assignment f2=f...\n";
  f2 = f; // check the assignment 
  std::cerr << "Auto-assignment f=f...\n";
  f2 = f2; // check the auto-assignment 
  std::cerr << "Copy-construction...\n";
  Map f3(f); // check the copy constructor

  std::cerr << "Check sizes: "
	    << "f.size()=" << f.size()
	    << " ,f2.size()=" << f2.size()
	    << " ,f3.size()=" << f3.size()
	    << "\n";
  assert(f.size() == f2.size());
  assert(f.size() == f3.size());

  std::cerr << "Emptying f, f2, and f3, with different methods...\n";
  int i2 = 0;
  while(!f.empty())
  {
    if(argc > 2)
    {
      std::cerr << "Display of f:\n";
      f.dump_direct_func(std::cerr);
    }
    int i = f.front()->second;
    if(argc > 2)
      std::cerr << "f.front()->second = " << i << "\n";
#ifndef NDEBUG
    bool i_is_in_f2 =
#endif
      f2.erase(i);
#ifndef NDEBUG
    bool i_is_in_f3 = 
#endif
      f3.erase(i);
    assert( i_is_in_f2 );
    assert( i_is_in_f3 );
    int new_i2 = f.front()->first;
    if(argc > 2)
      std::cerr << "f.front()->first = " << new_i2 << "\n";
    assert(new_i2>=i2);
    i2=new_i2;
    assert(i2 == i*i);
    f.pop_front();
  }

  std::cerr << "Check sizes: "
	    << "f.size()=" << f.size()
	    << " ,f2.size()=" << f2.size()
	    << " ,f3.size()=" << f3.size()
	    << "\n";
  assert(f.size()==0);
  assert(f2.empty());
  assert(f3.empty());

  std::cerr << "Filling f with f(i)=i*i, i=0.." << number_of_elements-1 << "...\n";
  for(unsigned int n=0; n<number_of_elements; n++)
  {
    f.insert(n, n*n);
  }

  std::cerr << "Check size: "
	    << "f.size()=" << f.size() << "\n";
  assert(f.size()==number_of_elements);

  std::cerr << "Emptying f...\n";
  i2=0;
  unsigned int counter = 0;
  while(!f.empty())
  {
    int i = f.front()->second;
    int new_i2 = f.front()->first;
    assert(new_i2>=i2);
    i2=new_i2;
    assert(i2 == i*i);
    f.pop_front();
    ++counter;
  }

  std::cerr << "Check size: "
	    << "f.size()=" << f.size() << "\n";
  assert(f.size()==0);
  assert(counter==number_of_elements);

  std::cerr << "Filling f with f(i)=i*i, i=0.." << number_of_elements -1 << "...\n";
  for(unsigned int n=0; n<number_of_elements; n++)
  {
    f.insert(n, n*n);
  }
  std::cerr << "f.clear()...\n";
  f.clear();
  std::cerr << "Check size: "
	    << "f.size()=" << f.size() << "\n";
  assert(f.size()==0);

  return 0;
}
