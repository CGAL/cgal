#include <CGAL/basic.h>
#include <CGAL/Double_map.h>
#include <utility>
#include <CGAL/Random.h>

#include <list>
#include <iostream>

typedef CGAL::Double_map<int, int> Map;

int main(int, char**)
{
  std::cerr << "Creating empty maps f and f2...\n";
  Map f;
  Map f2;

  std::cerr << "Filling f with 500 random integers...\n";
  for(int n=0; n<500; n++)
    {
      int i=CGAL::default_random.get_int(0,1000);
      f.insert(i, i*i);
    }

  std::cerr << "Check f.size<=500.\n";
  CGAL_assertion(f.size()<=500);

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
  CGAL_assertion(f.size() == f2.size());
  CGAL_assertion(f.size() == f3.size());

  std::cerr << "Emptying f, f2, and f3, with different methods...\n";
  int i, i2;
  i2 = 0;
  while(!f.empty())
    {
      i = f.front()->second;
      f2.erase(i);
      f3.erase(i);
      int new_i2 = f.front()->first;
      CGAL_assertion(new_i2>=i2);
      i2=new_i2;
      CGAL_assertion(i2 == i*i);
      f.pop_front();
    }

  std::cerr << "Check sizes: "
	    << "f.size()=" << f.size()
	    << " ,f2.size()=" << f2.size()
	    << " ,f3.size()=" << f3.size()
	    << "\n";
  CGAL_assertion(f.size()==0);
  CGAL_assertion(f2.empty());
  CGAL_assertion(f3.empty());

  std::cerr << "Filling f with f(i)=i*i, i=0..499 ...\n";
  for(int n=0; n<500; n++)
    {
      f.insert(n, n*n);
    }

  std::cerr << "Check size: "
	    << "f.size()=" << f.size() << "\n";
  CGAL_assertion(f.size()==500);

  std::cerr << "Emptying f...\n";
  i2=0;
  int counter = 0;
  while(!f.empty())
    {
      i = f.front()->second;
      int new_i2 = f.front()->first;
      CGAL_assertion(new_i2>=i2);
      i2=new_i2;
      CGAL_assertion(i2 == i*i);
      f.pop_front();
      ++counter;
    }

  std::cerr << "Check size: "
	    << "f.size()=" << f.size() << "\n";
  CGAL_assertion(f.size()==0);
  CGAL_assertion(counter==500);

  std::cerr << "Filling f with f(i)=i*i, i=0..499 ...\n";
  for(int n=0; n<500; n++)
  {
    f.insert(n, n*n);
  }
  std::cerr << "f.clear()...\n";
  f.clear();
  std::cerr << "Check size: "
	    << "f.size()=" << f.size() << "\n";
  CGAL_assertion(f.size()==0);

  return 0;
}
