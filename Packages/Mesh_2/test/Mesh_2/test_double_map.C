#include <CGAL/basic.h>
#include <CGAL/Double_map.h>
#include <utility>
#include <CGAL/Random.h>

#include <list>
#include <iostream>

typedef CGAL::Double_map<int, int> Map;

int main(int, char**)
{
  Map f;
  Map f2;

  for(int n=0; n<500; n++)
    {
      int i=CGAL::default_random.get_int(0,1000);
      f.insert(i, i*i);
    }

  CGAL_assertion(f.size()<=500);

  f2 = f; // check the assignment 
  f2 = f2; // check the auto-assignment 
  Map f3(f); // check the copy constructor

  CGAL_assertion(f.size() == f2.size());
  CGAL_assertion(f.size() == f3.size());

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

  CGAL_assertion(f.size()==0);
  CGAL_assertion(f2.empty());
  CGAL_assertion(f3.empty());

  for(int n=0; n<500; n++)
    {
      f.insert(n, n*n);
    }

  CGAL_assertion(f.size()==500);
  
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

   CGAL_assertion(f.size()==0);
   CGAL_assertion(counter==500);

   for(int n=0; n<500; n++)
    {
      f.insert(n, n*n);
    }
   f.clear();
   CGAL_assertion(f.size()==0);

  return 0;
}
