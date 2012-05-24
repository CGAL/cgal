#include "Volume_plane_intersection.h"
#include "Volume_plane_interface.h"

#include <CGAL/gl.h>


void Volume_plane_intersection::draw() const {
  glDisable(GL_LIGHTING);
  glLineWidth(4.0f);

  if(b && c) {
    glPushMatrix();
    glMultMatrixd(b->manipulatedFrame()->matrix());
    glMultMatrixd(c->manipulatedFrame()->matrix());

    glBegin(GL_LINES);
    {
      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f(x, 0.0f, 0.0f);
    }
    glEnd();

    glPopMatrix();
  }

  if(a && c) {
    glPushMatrix();
    glMultMatrixd(a->manipulatedFrame()->matrix());
    glMultMatrixd(c->manipulatedFrame()->matrix());

    glBegin(GL_LINES);
    {
      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f(0.0f, y, 0.0f);
    }
    glEnd();

    glPopMatrix();
  }

  if(a && b) {  
    glPushMatrix();
    glMultMatrixd(a->manipulatedFrame()->matrix());
    glMultMatrixd(b->manipulatedFrame()->matrix());

    glBegin(GL_LINES);
    {
      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f(0.0f, 0.0f, z);
    }
    glEnd();

    glPopMatrix();
  }

  glLineWidth(1.0f);
  glEnable(GL_LIGHTING);
}

#include "Volume_plane_intersection.moc"
