// MSVC workaround for non-template version of geowin_support.h

#include <CGAL/leda_integer.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Ray_2.h>

leda_point convert_to_leda(const CGAL::Point_2<CGAL::Cartesian<leda_integer> >& obj)
{
  double x = CGAL::to_double(obj.x());
  double y = CGAL::to_double(obj.y());
  leda_point p(x,y);
  return p;
}

ps_file& operator<<(ps_file& F,const CGAL::Point_2<CGAL::Cartesian<leda_integer> >& o) { F << convert_to_leda(o); return F; }

ps_file& operator<<(ps_file& F,const CGAL::Segment_2<CGAL::Cartesian<leda_integer> >& o) { return F; }
ps_file& operator<<(ps_file& F,const CGAL::Ray_2<CGAL::Cartesian<leda_integer> >& o) { return F; }

const char* leda_tname(CGAL::Point_2<CGAL::Cartesian<leda_integer> >* p) {  return "CGALPoint"; }

bool geowin_IntersectsBox(const CGAL::Point_2<CGAL::Cartesian<leda_integer> >& obj, double x1,double y1,double x2, double y2,bool f)
{
  double xw= CGAL::to_double(obj.x());
  double yw= CGAL::to_double(obj.y());  
  
  if (x1<=xw && x2>=xw && y1<=yw && y2>=yw) return true;
  return false;
}

void geowin_BoundingBox(const CGAL::Point_2<CGAL::Cartesian<leda_integer> >& obj,double& x1, double& x2,
	 double& y1, double& y2)
{
  CGAL::Bbox_2 bb= obj.bbox();
  x1= bb.xmin(); x2=bb.xmax(); y1=bb.ymin(); y2=bb.ymax();
}

void geowin_Translate(CGAL::Point_2<CGAL::Cartesian<leda_integer> >& obj, double dx, double dy)
{ 
  CGAL::Vector_2<CGAL::Cartesian<leda_integer> > vec;
  vec = CGAL::Vector_2<CGAL::Cartesian<leda_integer> >(leda_integer(dx),leda_integer(dy));
  obj = obj + vec;
}

void geowin_Rotate(CGAL::Point_2<CGAL::Cartesian<leda_integer> >& obj, double x, double y, double a)
{
  leda_point p2(CGAL::to_double(obj.x()), CGAL::to_double(obj.y()));
  p2 = p2.rotate(leda_point(x,y), a);
  obj = CGAL::Point_2<CGAL::Cartesian<leda_integer> >(leda_integer(p2.xcoord()), leda_integer(p2.ycoord())); 
}

// functions for the container

leda_string geowin_info_fcn(const std::list<CGAL::Point_2<CGAL::Cartesian<leda_integer> > >& L)
{
  leda_string str("~~~\\black \\tt STL-list of %d %ss", L.size(), " CGAL-point");  return str;
}

void geowin_generate_objects(GeoWin& gw, std::list<CGAL::Point_2<CGAL::Cartesian<leda_integer> > >& L)
{ 
  leda_list<leda_point> H;
  geowin_generate_objects(gw,H);

  //convert the contents
  CGAL::Point_2<CGAL::Cartesian<leda_integer> > p;
  leda_point mp;

  forall(mp,H){
   p= CGAL::Point_2<CGAL::Cartesian<leda_integer> >(leda_integer(mp.xcoord()), leda_integer(mp.ycoord()));
   L.push_front(p);
  }
}
