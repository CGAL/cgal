#ifndef CGAL_CONSTRAINED_TRIANGULATION_2_H
#define CGAL_CONSTRAINED_TRIANGULATION_2_H

#include <pair.h>
#include <list.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Constrained_triangulation_sweep.h>


template < class Gt, class Tds>
class CGAL_Constrained_triangulation_2
  : public CGAL_Triangulation_2<Gt,Tds>
{

public:
  typedef CGAL_Triangulation_2<Gt,Tds> Triangulation;
  typedef CGAL_Constrained_triangulation_2<Gt,Tds> Constrained_triangulation;
  typedef pair<Point,Point> Constraint;

  CGAL_Constrained_triangulation_2() : Triangulation() { }

  CGAL_Constrained_triangulation_2(const Gt& gt) : Triangulation(gt) { }






};

#endif CGAL_CONSTRAINED_TRIANGULATION_2_H
