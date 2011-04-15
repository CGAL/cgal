// Use something defined not in headers but in the CGAL library to test that is was indeed properly built and linked to,

#include <CGAL/IO/pixmaps/alpha_shape.xpm>
#include <CGAL/auto_link/CGALQt3.h>

int main()
{
  volatile const char * xpm = alpha_shape_xpm[0] ;
  
  return (xpm != 0) ? 0 : 1 ;
}
