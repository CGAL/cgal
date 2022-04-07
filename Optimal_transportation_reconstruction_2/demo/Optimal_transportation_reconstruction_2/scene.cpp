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

void Scene::draw_tolerance(GlViewer* viewer)
{
  double tol = m_pwsrec->tolerance();
  if (tol < 0.)
    return;

  const std::size_t resolution = 16;
  viewer->glColor3f (0.6f, 1.0f, 0.8f);
  std::vector<Sample_>::const_iterator it;
  for (it = m_samples.begin(); it != m_samples.end(); it++)
    {
      viewer->glBegin (GL_TRIANGLE_FAN);
      const Point& center = it->point();
      for (std::size_t i = 0; i < resolution; ++ i)
        {
          Point p = center + Vector (tol * std::cos (2. * M_PI * (i / (double)resolution)),
                                     tol * std::sin (2. * M_PI * (i / (double)resolution)));
          viewer->glVertex2d(p.x(), p.y());
        }
      viewer->glEnd();
    }
}
