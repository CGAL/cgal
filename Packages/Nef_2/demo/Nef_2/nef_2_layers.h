#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_Nef_2.h>

template <class Nef_polyhedron>
class Qt_layer_nef_blue : public CGAL::Qt_widget_layer
{
public:	
  Qt_layer_nef_blue(Nef_polyhedron &n): Nef(n){}
  void draw()
  {    
    *widget << CGAL::FillColor(CGAL::BLUE) << CGAL::GREEN;
    widget->setRasterOp(XorROP);
    *widget << Nef;
  };	
private:
  Nef_polyhedron &Nef;
};//end class 

template <class Nef_polyhedron>
class Qt_layer_nef_gray : public CGAL::Qt_widget_layer
{
public:
	
  Qt_layer_nef_gray(Nef_polyhedron &n): Nef(n){}
  void draw()
  {
    *widget << CGAL::FillColor(CGAL::GRAY) << CGAL::WHITE;
    widget->setRasterOp(XorROP);
    *widget << Nef;
  }
private:
  Nef_polyhedron &Nef;
};//end class 
