#include "Viewer.h"
#include <CGAL/point_generators_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Exact_spherical_kernel_3.h>
#include <vector>


template <class Kernel>
void Viewer::draw_circle_on_unit_sphere(const typename CGAL::Point_3<Kernel>& center) const {
  const typename Kernel::Point_3 origin(0,0,0);
  const typename Kernel::Plane_3 plane(center, center-origin);
  typename Kernel::Vector_3 base1=plane.base1();
  typename Kernel::Vector_3 base2=plane.base2();
  base1=base1/CGAL::sqrt(base1.squared_length());
  base2=base2/CGAL::sqrt(base2.squared_length());
  const double radius=CGAL::sqrt( CGAL::to_double( 1 - CGAL::squared_distance(origin,center) ) );
  const double nb_pt_per_circle=100;
  const double step=2 * CGAL_PI / nb_pt_per_circle;
  ::glDisable(GL_LIGHTING);
  ::glBegin(GL_LINE_LOOP);
  for (double theta = 0; theta < 2 * CGAL_PI ; theta += step) {
    const typename Kernel::Point_3 point=center + ( radius*cos(theta)*base1 + radius*sin(theta)*base2 );
    ::glVertex3f( point.x(),point.y(),point.z() );
  }
  ::glEnd();
  ::glEnable(GL_LIGHTING);
}

void Viewer::draw()
{
  const int sphere_res=30;
  
  //draw central sphere
  ::glPushMatrix();
  ::glColor3f(1,1,1);
  ::gluSphere(qsphere,0.999,sphere_res,sphere_res);
  ::glCallList(dl_nb);
  ::glCallList(dl_nb+1);
  ::glPopMatrix();
}

void Viewer::init()
{
  // Restore previous viewer state.
  restoreStateFromFile();
  
  //init quadric to store sphere.
  qsphere=gluNewQuadric();
  gluQuadricOrientation(qsphere,GLU_OUTSIDE);
  
  //code for antialiasing
  ::glEnable(GL_BLEND);
  ::glEnable(GL_LINE_SMOOTH);
  ::glEnable(GL_POINT_SMOOTH);
  ::glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  ::glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  

  //set intersection point size
  ::glPointSize(5);
  
  //random generator of points within a sphere
  typedef CGAL::Creator_uniform_3<EPIC::FT,EPIC::Point_3>   Creator;
  CGAL::Random_points_in_sphere_3<EPIC::Point_3, Creator>   gen;
  
  const unsigned nb_circles=20;
  
  //vector to store input points
  std::vector<EPIC::Point_3> points;
  points.reserve(nb_circles);
  
  //disable lighting before drawing segments and points
  ::glEnable(GL_LIGHTING);
  
  //init call lists
  dl_nb=::glGenLists(2);
  
  //store circles in a display list
  ::glNewList(dl_nb,GL_COMPILE);
  ::glColor3f(1,0,0);
  
  for (unsigned i=0;i<nb_circles;++i){
    EPIC::Point_3 p=*++gen;
    //prevent great circles
    while (p.x()==0 && p.y()==0 && p.z()==0) {  p=*++gen; }
    draw_circle_on_unit_sphere(p);
    points.push_back(p);
  }
  ::glEndList();
  
  std::vector<EPIC::Point_3> intersections;
  naive_compute_intersection_points(points,std::back_inserter(intersections));

  ::glNewList(dl_nb+1,GL_COMPILE);
  ::glColor3f(0,1,0);
  //draw points as small spheres
  for (std::vector<EPIC::Point_3>::const_iterator it=intersections.begin();it!=intersections.end();++it){
    glPushMatrix();
    glTranslatef(it->x(),it->y(),it->z());
    gluSphere(qsphere,0.005,10,10);
    glPopMatrix();    
  }
  ::glEndList();
  
  //lighting
  ::glEnable(GL_LIGHT0);
  ::glEnable(GL_LIGHTING);
  float lpos[4] = { -.2f, .2f, .9797958971f, 0.0f };
  ::glLightfv(GL_LIGHT0,GL_POSITION,lpos);
  ::glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0f);

  ::glEnable(GL_NORMALIZE);
  ::glShadeModel(GL_SMOOTH);
}


template<class Output_iterator>
void Viewer::naive_compute_intersection_points(const std::vector<EPIC::Point_3>& points,Output_iterator out) const {
  typedef CGAL::Exact_spherical_kernel_3 SK;
  SK::Sphere_3 sphere(SK::Point_3(0,0,0),1);
  
  //converter point to exact SK point type
  CGAL::Cartesian_converter<EPIC,SK> to_exact;
  std::vector<SK::Circle_3> circles;
  
  //create circles from points: need to use a converter to change floating points coordinates into an exact NT.
  for (std::vector<EPIC::Point_3>::const_iterator it=points.begin();it!=points.end();++it){
    const SK::Point_3 center=to_exact(*it);
    circles.push_back( SK::Circle_3(sphere,SK::Plane_3(center,center-CGAL::ORIGIN) ) );
  }
  
  
  //Look for intersection points among pair of circles: use a naive and quadratic way
  for (std::vector<SK::Circle_3>::const_iterator it_f=circles.begin();it_f!=--circles.end();++it_f){
    std::vector<SK::Circle_3>::const_iterator it_s=it_f;
    ++it_s;
    for (;it_s!=circles.end();++it_s){
      std::vector <CGAL::Object> intersections;
      //ensure_circles are different
      CGAL_precondition(*it_s!=*it_f);
      CGAL::intersection(*it_f,*it_s,std::back_inserter(intersections));
      if (!intersections.empty()){
        for (std::vector <CGAL::Object>::const_iterator it_pt=intersections.begin();it_pt!=intersections.end();++it_pt){
          const std::pair<SK::Circular_arc_point_3,unsigned>* pt=
            CGAL::object_cast< std::pair<SK::Circular_arc_point_3,unsigned> > (&(*it_pt));
          assert(pt!=NULL);
          *out++=EPIC::Point_3( CGAL::to_double(pt->first.x()),
                                CGAL::to_double(pt->first.y()),
                                CGAL::to_double(pt->first.z()) 
          );
        }
      }
    }
  }
}
