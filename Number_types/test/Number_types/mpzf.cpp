#include <CGAL/basic.h>
#include <CGAL/mpzf.h>

using CGAL::mpzf;
int main(){
  mpzf a=3;
  mpzf b=4.5;
  mpzf c=2*(a+b)-5*a;
  assert(CGAL::sign(c)==0);
}
