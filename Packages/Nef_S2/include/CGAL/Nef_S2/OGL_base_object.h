#ifndef CGAL_OGL_BASE_OBJECT_H
#define CGAL_OGL_BASE_OBJECT_H

#include <CGAL/basic.h>
#include <GL/glut.h>

CGAL_BEGIN_NAMESPACE

namespace OGL {

  class OGL_base_object {    
  public:

    typedef CGAL::Simple_cartesian<double> DKernel;  
    typedef DKernel::Point_3               Double_point;
    typedef DKernel::Vector_3              Double_vector;
    typedef DKernel::Segment_3             Double_segment;
    typedef DKernel::Aff_transformation_3  Affine_3;

    virtual void draw() const  = 0;
    virtual void init() = 0;
    virtual void toggle(int) = 0;
    virtual void set_style(int) = 0;
    virtual ~OGL_base_object() {}
  };
}

CGAL_END_NAMESPACE
#endif // CGAL_OGL_BASE_OBJECT_H
