#ifndef CGAL_CONSTRAINED_TRIANGULATION_DEMO_2_H
#define CGAL_CONSTRAINED_TRIANGULATION_DEMO_2_H

#include <CGAL/Constrained_triangulation_2.h>
#include <CGAL/Constrained_triangulation_sweep_demo_2.h>
#include <CGAL/IO/Window_stream.h>

// template < class Gt,class Tds>
// class CGAL_Constrained_triangulation_sweep_demo_2;

template < class Gt,class Tds>
class CGAL_Constrained_triangulation_demo_2
  : public CGAL_Constrained_triangulation_2<Gt,Tds>
{
public:
  typedef CGAL_Constrained_triangulation_2<Gt,Tds> Constrained_triangulation;
  typedef CGAL_Constrained_triangulation_sweep_demo_2<Gt,Tds>  Sweep_demo;
  typedef typename Gt::Segment Segment;
  typedef CGAL_Window_stream Window_stream;

CGAL_Constrained_triangulation_demo_2() :  
  CGAL_Constrained_triangulation_2<Gt,Tds>() {} 
  
CGAL_Constrained_triangulation_demo_2(const Gt& gt=Gt()) 
  : CGAL_Constrained_triangulation_2<Gt,Tds>(gt) {}
  

CGAL_Constrained_triangulation_demo_2(Window_stream& W,
				     list<Constraint>& lc, const Gt& gt=Gt()) 
  : CGAL_Constrained_triangulation_2<Gt,Tds>(gt)
{
  Sweep_demo  sweep(W,lc, gt);
  CGAL_Constrained_triangulation_2<Gt,Tds> Tr( sweep.vertex(), gt);
  swap(Tr);
  CGAL_triangulation_postcondition( is_valid() );
}

};


#endif //CGAL_CONSTRAINED_TRIANGULATION_DEMO_2_H
