#include <CGAL/event.h>
#include <iostream>

using namespace CGAL;
using std::cout;

void by_reference(int& val)
{ cout << "Call by reference:  " << val <<  " Change value! \n";
  val = 3;
}

void by_value(int val)
{ cout << "Call by value:      " << val << "\n"; }

int main()
{
    event e;
  
    event_item IT[2];
    IT[0] = attach(e,by_value);
    IT[1] = attach(e,by_reference);
  
    int i = 10;  
    occur(e,i);
    occur<int&>(e,i);  
    occur(e,i);
    
    event e2;
    occur(e2);
    
    detach(IT, 1);
    cout << "after detach ...\n"; cout.flush();
    
    occur(e,i);
    occur<int&>(e,i);  
    occur(e,i);    
  
    return 0;
}
