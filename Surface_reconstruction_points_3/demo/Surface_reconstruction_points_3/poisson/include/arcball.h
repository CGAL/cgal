#pragma once
#include "Director.h"

class Arcball 
{
private:
  // xside = half of the width of the window
  // yside = half of the height of the window
  // side = min(xside, yside)
  int side;
  int xside;
  int yside;
  float size;
  Director* b;

public:

  const float* get_camera_position()
  {
    return b->center();
  }

  Director* setDirector(float max_size,
    float cx, 
    float cy, 
    float cz)
  {
    Director* b = new Director();
    b->radius() = 1.0f/max_size;
    b->setOffset(cx, cy, cz);
    return b;
  }

  Arcball(int width,
    int height,
    float max_size, 
    float cx,
    float cy,
    float cz )
  {
    size = max_size;
    b = setDirector(max_size, cx, cy, cz);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.,1.,-1.,1.,-1000.,1000.); 
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    resize(width, height);
  }

  ~Arcball()
  {
    delete b;
  }

  void resize(int x,int y) {
    float r = (float)y/(float)x;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glViewport(0, 0, x, y);
    if (x<y){
      side = x/2;
      glOrtho(-1., 1., -r, r, -1000., 1000.); 
    }
    else {
      glOrtho(-1./r, 1./r, -1., 1., -1000., 1000.); 
      side = y/2;
    }
    xside = x/2;
    yside = y/2;

    glMatrixMode(GL_MODELVIEW);
  }

  void draw_sphere()
  {
    const float* offset = b->getOffset();
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

    glPushMatrix();

    glDisable(GL_LIGHTING);

    glTranslated(-offset[0], -offset[1], -offset[2]);

    const float scale = 1.2f;
    glScaled(size * scale, size * scale, size * scale);

    glColor3d(0., 0., 1.);
    unitCircle(1.);

    glRotated(90., 1., 0., 0.);
    glColor3d(0., 1., 0.);
    unitCircle(1.);

    glRotated(90., 0., 1., 0.);
    glColor3d(1., 0., 0.);
    unitCircle(1.);

    glPopMatrix();
  }

  void draw_sphere_scene()
  {
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

    glPushMatrix();

    glDisable(GL_LIGHTING);

    const float scale = 1.f;
    glScaled(scale, scale, scale);

    glColor3d(0., 0., 1.);
    unitCircle(1.);

    glRotated(90., 1., 0., 0.);
    glColor3d(0., 1., 0.);
    unitCircle(1.);

    glRotated(90., 0., 1., 0.);
    glColor3d(1., 0., 0.);
    unitCircle(1.);

    glPopMatrix();
  }


  void unitCircle(const double r)
  {
    glBegin(GL_LINE_LOOP);
    for(double i = 0.; i < 320.; ++i)
      glVertex3d(r*cos(i/160.*3.1416), r*sin(i/160.*3.1416), 0.);
    glEnd();
  }

  void setCamera() {
    const float* offset;
    offset = b->getOffset();

    glMatrixMode(GL_MODELVIEW);

    glLoadMatrixf(b->mat());
    glTranslated(offset[0], offset[1], offset[2]);
  }

  void setCamera_no_offset_no_scale() {
    glMatrixMode(GL_MODELVIEW);

    glLoadMatrixf(b->matSansTransZoom());
  }


  void zoom(float factor) {
    b->radius() *= factor;
  }

  void recenter() {
    b->recenter(0., 0., 0.);
  }

  void reset()
  {
    b->radius() = 1.0f/size;
    recenter();
  }

  void translate(int xold, int yold, int x, int y) {
    float Xold = ((float) xold - xside)/side;
    float Yold = ((float) yold - yside)/side;
    float X = ((float) x - xside)/side;
    float Y = ((float) y - yside)/side;

    float vx = X - Xold;
    float vy = Y - Yold;
    b->translation(vx, -vy, 0.);
  }

  void rotate_obj(int xold, int yold, int x, int y) 
  {
    float Xold = ((float) xold - xside)/side;
    float Yold = ((float) yold - yside)/side;
    float X = ((float) x - xside)/side;
    float Y = ((float) y - yside)/side;
    float imCenter[3];
    const float* offset = b->getOffset();

    imCenter[0] = -offset[0];
    imCenter[1] = -offset[1];
    imCenter[2] = -offset[2];
    b->image(imCenter);

    const float scale = 1.2f;

    b->rotationObj((Xold - imCenter[0])/(b->radius()*size*scale), 
      (-Yold - imCenter[1])/(b->radius()*size*scale), 
      (X - imCenter[0])/(b->radius()*size*scale), 
      (-Y - imCenter[1])/(b->radius()*size*scale));
  }

  void rotate_scene(int xold, int yold, int x, int y) {
    float Xold = ((float) xold - xside)/side;
    float Yold = ((float) yold - yside)/side;
    float X = ((float) x - xside)/side;
    float Y = ((float) y - yside)/side;

    b->rotation(Xold, -Yold, X, -Y);

  }
};

