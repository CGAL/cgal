// argv[1] = File to store the disks. Defaults to stdout, in which case no 
// constraints are input
// argv[2] = File to store the constraints. If not given, no constraints
// are input.

// First input the disks. Then press enter. Then input the constraints.
// To begin a left_xx (respectively right_xx) constraint, click the source
// disk with the left (respectively right button). Then enter the target
// disk similarily.

// It is possible to pop the disks or constraints previously input
// by pressing backspace.

#include<CGAL/basic.h>

#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*){
  std::cout << "Sorry, this demo needs QT..." << std::endl; return 0;}
#else

#include<CGAL/Point_2.h>
#include<CGAL/Cartesian.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <iostream>
#include <fstream>
#include <string>


#include<CGAL/IO/Qt_widget.h>
#include<CGAL/Visibility_complex_2.h>



#include<algorithm>
#include<cmath>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Gmpz.h>


typedef CGAL::Gmpz N;
typedef CGAL::Cartesian<N> K;

#include <CGAL/Point_2.h>
typedef CGAL::Point_2<K> Point_2;

template<class Disk> CGAL::Bbox_2 bbox(const Disk&d) {
  return d.bbox();
}

#ifdef VC_SCENE_COMPUTE_BIT_CIRCLE

#include<CGAL/Visibility_complex_2/Circle_traits.h>
typedef CGAL::Visibility_complex_2_circle_traits<K> Gt;

#elif defined(VC_SCENE_COMPUTE_BIT_POLYGON)

#include<CGAL/Visibility_complex_2/Polygon_traits.h>
typedef CGAL::Visibility_complex_2_polygon_traits<K> Gt;


#include <CGAL/IO/Qt_widget_Polygon_2.h>


#elif defined(VC_SCENE_COMPUTE_BIT_SEGMENT)
#include<CGAL/Visibility_complex_2/Segment_traits.h>
typedef CGAL::Visibility_complex_2_segment_traits<K> Gt;


#elif defined(VC_SCENE_COMPUTE_BIT_POINT)
#include<CGAL/Visibility_complex_2/Point_traits.h>
typedef CGAL::Visibility_complex_2_point_traits<K> Gt;

#elif defined(VC_SCENE_COMPUTE_BIT_ELLIPSE)
#include<CGAL/Visibility_complex_2/Ellipse_traits.h>
#include<CGAL/IO/Qt_widget_Conic_2.h>
typedef CGAL::Visibility_complex_2_ellipse_traits<K> Gt;

void compute_x_span(double r,double s,double t,double u,double v,double w,
           double&xc,double&width) {
  double a=-4*s*r+t*t;
  double b=v*t-2*s*u;
  double c=v*v-4*s*w;
  xc=-b/a;
  width=CGAL_NTS sqrt(b*b-a*c)/CGAL_NTS abs(a);
}

CGAL::Bbox_2 bbox(Gt::Disk&d) {
  double r=CGAL_NTS to_double(d.r());
  double s=CGAL_NTS to_double(d.s());
  double t=CGAL_NTS to_double(d.t());
  double u=CGAL_NTS to_double(d.u());
  double v=CGAL_NTS to_double(d.v());
  double w=CGAL_NTS to_double(d.w());

  double xc,yc,width,height;
  compute_x_span(r,s,t,u,v,w,xc,width);
  compute_x_span(s,r,t,v,u,w,yc,height);
  return CGAL::Bbox_2(xc-width,yc-height,xc+width,yc+height);
}
#endif


typedef Gt::Bitangent_2 Bitangent_2;
typedef CGAL::Compute_free_bitangents_2<Gt> VC;
typedef VC::Constraint_input Constraint_input;


#include <CGAL/IO/Ostream_iterator.h>




typedef CGAL::Point_2<K> Point_2;
typedef Gt::Segment_2 Segment_2;

int main(int argc,char ** argv) {
  std::istream * di=&std::cin;
  if (argc>1) di=new std::ifstream(argv[1]);
  std::istream_iterator<Gt::Disk> disk_it(*di),disk_end;
  std::vector<Gt::Disk> disks(disk_it,disk_end);
  std::istream_iterator<Constraint_input> 
    constraint_it,constraint_end;
  std::istream * ci=0;
  if (argc>2) {
    ci=new std::ifstream(argv[2]);
    constraint_it=
      std::istream_iterator<Constraint_input>(*ci);
  }
  std::vector<Bitangent_2> vb;
  VC()(disks.begin(),disks.end(),
       constraint_it,constraint_end,std::back_inserter(vb));
  double xmin=999999999,xmax=-999999999,ymin=999999999,ymax=-999999999;
  for (std::vector<Gt::Disk>::iterator i=disks.begin();i!=disks.end();++i) {
    CGAL::Bbox_2 box=bbox(*i);
    xmin=std::min(box.xmin(),xmin);
    xmax=std::max(box.xmax(),xmax);
    ymin=std::min(box.ymin(),ymin);
    ymax=std::max(box.ymax(),ymax);
  }
  double dx=xmax-xmin;
  double dy=ymax-ymin;
  int width=600;
  int height=(int)(width*dy/dx);

  int ac=1;
  const char * ctitle="Display free bitangents";
  char title[24];
  std::copy(ctitle,ctitle+24,title);
  char * av[1]={title};
  QApplication app(ac,av);
  CGAL::Qt_widget* w;
  w = new CGAL::Qt_widget();
  app.setMainWidget( w );
  w->resize(width,height);
  w->set_window(xmin-0.1*dx,xmax+0.1*dx,ymin-0.1*dy,ymax+0.1*dy);
  w->show();
  w->lock();
  *w << CGAL::BLACK;
  std::copy(disks.begin(),disks.end(),
            CGAL::Ostream_iterator<Gt::Disk,CGAL::Qt_widget>(*w));
  *w << CGAL::RED;
  std::copy(vb.begin(),vb.end(),
            CGAL::Ostream_iterator<Segment_2,CGAL::Qt_widget>(*w));
  w->unlock();
  std::cout<<vb.size()<<"\n";
  return app.exec();
}
#endif
