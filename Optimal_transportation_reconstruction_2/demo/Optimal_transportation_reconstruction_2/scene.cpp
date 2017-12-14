#include "scene.h"
#include "glviewer.h"
void Scene::draw_samples(const float point_size, GlViewer* viewer) {

  viewer->glPointSize(point_size);
  viewer->glBegin(GL_POINTS);

  std::vector<Sample_>::const_iterator it;
  for (it = m_samples.begin(); it != m_samples.end(); it++) {
    double mass = it->mass();

    float value = mass;
    float grey = 0.9 * (1.0f - value);
    viewer->glColor3f(grey, grey, grey);
    const Point& p = it->point();
    viewer->glVertex2d(p.x(), p.y());
  }
  viewer->glEnd();
}
