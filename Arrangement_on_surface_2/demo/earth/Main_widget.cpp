
#include "Main_widget.h"

#include <cmath>
#include <iostream>
#include <string>

#include <QMouseEvent>

#include "Geodesic_arcs.h"
#include "Tools.h"


Main_widget::~Main_widget()
{
  // Make sure the context is current when deleting the texture and the buffers.
  makeCurrent();
  doneCurrent();
}

void Main_widget::set_mouse_button_pressed_flag(QMouseEvent* e, bool flag)
{
  switch (e->button())
  {
  case Qt::LeftButton:
    m_left_mouse_button_down = flag;
    break;

  case Qt::MiddleButton:
    m_middle_mouse_button_down = flag;
    break;
  }
}
void Main_widget::mousePressEvent(QMouseEvent* e)
{
  set_mouse_button_pressed_flag(e, true);
  m_mouse_press_pos = m_last_mouse_pos = QVector2D(e->position());

  // for the backprojected diff-vector method:
  if (m_left_mouse_button_down)
  {
    m_camera.save_config();
  }
}
void Main_widget::mouseMoveEvent(QMouseEvent* e)
{
  auto current_mouse_pos = QVector2D(e->position());
  const auto diff = current_mouse_pos - m_last_mouse_pos;

  if (m_left_mouse_button_down)
  {
    const float rotation_scale_factor = 0.1f;

    if(1)
    {
      // OUR CUSTOM AD-HOC CAMERA ROTATION
      m_theta += rotation_scale_factor * diff.x();
      m_phi += rotation_scale_factor * diff.y();
      m_camera.rotate_from_init_config(-m_theta, -m_phi);
    }
    else
    {
      // ROTATION AROUND AN AXIS ORTHOGONAL TO THE BACKPROJECTED DIF-VECTOR
      //QVector3D p0(m_last_mouse_pos.x(), m_vp_height - m_last_mouse_pos.y(), 0);
      QVector3D p0(m_mouse_press_pos.x(), m_vp_height - m_mouse_press_pos.y(), 0);
      QVector3D p1(current_mouse_pos.x(), m_vp_height - current_mouse_pos.y(), 0);
      auto dp = p1 - p0; // difference vector in OpenGL window coords.
      QVector3D rdp(-dp.y(), dp.x(), 0); // rotate diff-vector CCW by 90-deg
      QVector3D rp = p0 + rdp; // r1 rotated CCW by 90 deg
     
      QMatrix4x4 model; // this is different from Sphere's model matrix!!!
      auto proj = m_camera.get_projection_matrix();
      auto view = m_camera.get_view_matrix();
      auto model_view = view * model;
      QRect viewport(0, 0, m_vp_width, m_vp_height);
      auto wp0 = p0.unproject(model_view, proj, viewport);
      auto wrp = rp.unproject(model_view, proj, viewport);

      // rotation axis & angle
      auto rot_axis = wrp - wp0;
      rot_axis.normalize();
      const auto rot_angle = rotation_scale_factor * dp.length();
    
      QMatrix4x4 rot_matrix;
      rot_matrix.rotate(-rot_angle, rot_axis);
      
      m_camera.rotate_from_saved_config(rot_matrix);
    }
  }
  else if(m_middle_mouse_button_down)
  {
    const float zoom_scale_factor = 0.01f;
    const auto distance = zoom_scale_factor * diff.y();
    m_camera.move_forward(distance);
  }

  m_last_mouse_pos = current_mouse_pos;
}
void Main_widget::mouseReleaseEvent(QMouseEvent* e)
{
  set_mouse_button_pressed_flag(e, false);
}
void Main_widget::timerEvent(QTimerEvent*)
{
  update();
}

#include <qfile.h>
#include <qiodevice.h>
#include <QJsonArray.h>
#include <qjsonobject.h>
#include <QJsonParseError>
#include <qstring.h>


#include <QtXml/qdom.h>
#include <qxmlstream.h>

namespace {
  struct Node
  {
    double lon, lat;
  };

  struct LinearRing
  {
    std::vector<Node> nodes;
    QString str;
  };

  struct MultiGeometry
  {
    std::vector<LinearRing> polygons;
  };
  
  struct Placemark
  {
    MultiGeometry geometry;
    QString name;
  };
}


void Main_widget::initializeGL()
{
  QString file_name("C:/work/gsoc2023/data/world_countries.kml");
  
  if(0)
  {
    QDomDocument doc("mydoc");
    QFile file(file_name);
    if (!file.open(QFile::ReadOnly | QFile::Text)) 
    {
      qDebug() << "could not open file!";
      return;
    }
    if (!doc.setContent(&file)) {
      file.close();
      return;
    }
    file.close();

    //QDomElement docElem = doc.documentElement();
    //QDomNode n = docElem.firstChild();
    //qDebug() << n.nodeName();
    //qDebug() << n.isElement();
    //auto doc_elem = n.toElement();

    auto placemarks = doc.elementsByTagName("Placemark");
    qDebug() << placemarks.count();

  }
  else
  {
    std::vector<Placemark> placemarks;

    Placemark placemark;
    MultiGeometry mgeometry;
    LinearRing   lring;


    QFile file(file_name);
    if (file.open(QIODevice::ReadOnly)) 
    {
      QXmlStreamReader xmlReader;
      xmlReader.setDevice(&file);
      
      xmlReader.readNext();
      
      // Reading from the file
      while (!xmlReader.isEndDocument())
      {
        QString name = xmlReader.name().toString();
       // qDebug() << "----------------------";
       // qDebug() << name;
       // qDebug() << xmlReader.text();


        if (xmlReader.isStartElement())
        {
         // qDebug() << "START ELEMENT";
          if (name == "Placemark")
          {
           // qDebug() << "Placemark - Start";
            placemark = Placemark{};
          }
          else if (name == "MultiGeometry")
          {
           // qDebug() << "MultiGeometry - Start";
            mgeometry = MultiGeometry{};
          }
          else if (name == "LinearRing")
          {
           // qDebug() << "LinearRing - Start";
            lring = LinearRing{};
          }
          else if (name == "coordinates")
          {
           // qDebug() << "coordinates - Start";
            xmlReader.readNext();
            lring.str = xmlReader.text().toString();
           // qDebug() << lring.str;
          }
          else if (name == "SimpleData")
          {
           // qDebug() << "SimpleData - Start";
            auto attributes = xmlReader.attributes();
            auto attr_name = attributes[0].name().toString();
            auto attr_value = attributes[0].value().toString();
            if ((attr_name == "name") && (attr_value == "name"))
            {
              xmlReader.readNext();
              placemark.name = xmlReader.text().toString();
             // qDebug() << "country name = " << placemark.name;
            }
            //qDebug() << "num attribues = " << attributes.size();
            //qDebug() << "attribute[0].name = " << attributes[0].name();
            //qDebug() << "attribute[0].value = " << attributes[0].value();
          }
        }
        else if (xmlReader.isEndElement())
        {
         // qDebug() << "END ELEMENT";
          if (name == "Placemark")
          {
           // qDebug() << "Placemark - End";
            placemarks.push_back(placemark); // move?
          }
          else if (name == "MultiGeometry")
          {
           // qDebug() << "MultiGeometry - End";
            placemark.geometry = mgeometry; // move?
          }
          else if (name == "LinearRing")
          {
           // qDebug() << "LinearRing - End";
            mgeometry.polygons.push_back(lring); // move?
          }
          else if (name == "coordinates")
          {
           // qDebug() << "coordinates - End";
            // no need to do anything here: the coordinates are read above!
          }
        }

        xmlReader.readNext();
      }

      if (xmlReader.hasError())
      {
        std::cout << "XML error: " << xmlReader.errorString().data() << std::endl;
      }
    }
  }

  initializeOpenGLFunctions();

  init_camera();
  init_geometry();
  init_shader_programs();

  {
    // TO-DO: move this code to resizeGL (when viewport is initialized)
    // has to be defined after camera has been defined:
    // because we want to compute the error based on camera parameters!
    Geodesic_arcs ga;
    const double error = 0.001; // calculate this from cam parameters!
    auto lsa = ga.get_approximate_arcs(error);
    m_geodesic_arcs = std::make_unique<Line_strips>(lsa);
  }

  glClearColor(0, 0, 0, 1);
  glEnable(GL_DEPTH_TEST);  // Enable depth buffer
  //glEnable(GL_CULL_FACE); // Enable back face culling

  // Use QBasicTimer because its faster than QTimer
  m_timer.start(12, this);
}



void Main_widget::init_camera()
{
  m_camera.set_pos(0, 0, 3);
  //m_camera.rotate_around_x(-90);
}
void Main_widget::init_geometry()
{
  int num_slices, num_stacks;
  num_slices = num_stacks = 64;
  float r = 1;
  m_sphere = std::make_unique<Sphere>(num_slices, num_stacks, r);
  const float c = 0.8;
  m_sphere->set_color(c, c, c, 1);


  const float axes_length = 2;
  m_world_coord_axes = std::make_unique<World_coord_axes>(axes_length);
}
void Main_widget::init_shader_programs()
{
  Shader_program::set_shader_path("shaders/");
  m_sp_smooth.init_with_vs_fs("smooth");;
  m_sp_per_vertex_color.init_with_vs_fs("per_vertex_color");
  m_sp_arc.init_with_vs_fs("arc");
}

void Main_widget::find_minimum_projected_error_on_sphere(float we)
{
  QRect vp(0, 0, m_vp_width, m_vp_height);
  auto proj = m_camera.get_projection_matrix();
  auto view = m_camera.get_view_matrix();
  QMatrix4x4 model;
  auto model_view = view * model;

  float max_err = 0;
  float max_theta = -1;
  float max_phi = -1;

  int num_divs = 200;
  const float dtheta = M_PI_2 / num_divs;
  const float dphi = M_PI_2 / num_divs;

  const float r1 = 1.f;
  const float r2 = r1 - we;
  for (int i = 0; i <= num_divs; i++)
  {
    const float theta = dtheta * i;
    const float cos_theta = std::cos(theta);
    const float sin_theta = std::sin(theta);

    for (int j = 0; j <= num_divs; j++)
    {
      QVector3D p1, p2;
      const float phi = dphi * j;
      const float cos_phi = std::cos(phi);
      const float sin_phi = std::sin(phi);

      // p1
      const float r1xz = r1 * sin_phi;
      p1.setY(r1 * cos_phi);
      p1.setX(r1xz * cos_theta);
      p1.setZ(r1xz * sin_theta);

      // p2
      const float r2xz = r2 * sin_phi;
      p2.setY(r2 * cos_phi);
      p2.setX(r2xz * cos_theta);
      p2.setZ(r2xz * sin_theta);

      auto wp1 = p1.project(model_view, proj, vp);
      auto wp2 = p2.project(model_view, proj, vp);

      const auto pe = wp1.distanceToPoint(wp2);
      if (max_err < pe)
      {
        max_err = pe;
        max_theta = theta;
        max_phi = phi;
      }
    }
  }

  std::cout << "max err = " << max_err << std::endl;
  std::cout << "max phi = " << max_phi * 180 / M_PI << std::endl;
  std::cout << "max theta = " << max_theta * 180 / M_PI << std::endl;

  auto wp1 = QVector3D(0, r1, 0).project(model_view, proj, vp);
  auto wp2 = QVector3D(0, r2, 0).project(model_view, proj, vp);
  auto pe = wp1.distanceToPoint(wp2);
  std::cout << "polar err = " << pe << std::endl;

  wp1 = QVector3D(r1, 0, 0).project(model_view, proj, vp);
  wp2 = QVector3D(r2, 0, 0).project(model_view, proj, vp);
  pe = wp1.distanceToPoint(wp2);
  std::cout << "x-axis err = " << pe << std::endl;

  wp1 = QVector3D(0, 0, 1).project(model_view, proj, vp);
  wp2 = QVector3D(we, 0, 1).project(model_view, proj, vp);
  pe = wp1.distanceToPoint(wp2);
  std::cout << "nearest proj err = " << pe << std::endl;

  wp1 = QVector3D(0, 0, -1).project(model_view, proj, vp);
  wp2 = QVector3D(we, 0, -1).project(model_view, proj, vp);
  pe = wp1.distanceToPoint(wp2);
  std::cout << "farthest proj err = " << pe << std::endl;

  // project the origin on the screen (to check if it projects to the mid-vp)
  //std::cout << QVector3D(0, 0, 0).project(model_view, proj, vp) << std::endl;
}

void Main_widget::resizeGL(int w, int h)
{
  m_vp_width = w;
  m_vp_height = h;

  // Reset projection
  qreal aspect = qreal(w) / qreal(h ? h : 1);
  const qreal z_near = 1.0, z_far = 100.0, fov = 45.0;
  m_camera.perspective(fov, aspect, z_near, z_far);

  {
    // compute the back-projected error
    QRect vp(0, 0, m_vp_width, m_vp_height);
    auto proj = m_camera.get_projection_matrix();
    auto view = m_camera.get_view_matrix();
    QMatrix4x4 model;
    auto model_view = view * model;

    QVector3D p0(m_vp_width / 2, m_vp_height / 2, 0);
    QVector3D p1(p0.x() + 1, p0.y(), 0);
    auto wp0 = p0.unproject(model_view, proj, vp);
    auto wp1 = p1.unproject(model_view, proj, vp);
    const float z_near = m_camera.get_z_near();
    const float r = 1.f; // sphere radius
    const QVector3D origin(0, 0, 0);
    const float dist_to_cam = m_camera.get_pos().distanceToPoint(origin);
    
    float d = dist_to_cam - r;
    float err = wp0.distanceToPoint(wp1) * (d / z_near);
    std::cout << "error = " << err << std::endl;


    // find the minimum error over the sphere
    //find_minimum_projected_error_on_sphere(err);
  }

}
void Main_widget::paintGL()
{
  QMatrix4x4 model;
  model.rotate(-90, 1,0,0); // this makes z-axes point upwards!
  const auto view = m_camera.get_view_matrix();
  const auto projection = m_camera.get_projection_matrix();
  const auto mvp = projection * view * model;

  // Clear color and depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // SPHERE
  {
    glEnable(GL_DEPTH_TEST);

    auto& sp = m_sp_smooth;
    sp.use();
    sp.set_uniform("u_mvp", mvp);
    sp.set_uniform("u_color", m_sphere->get_color());
    
    m_sphere->draw();

    sp.unuse();
  }

  // WORLD COORDINATE AXES
  {
    auto& sp = m_sp_per_vertex_color;
    sp.use();
    sp.set_uniform("u_mvp", mvp);

    m_world_coord_axes->draw();

    sp.unuse();
  }

  // GEODESIC ARCS
  {
    glDisable(GL_DEPTH_TEST);

    auto& sp = m_sp_arc;
    sp.use();
    sp.set_uniform("u_mvp", mvp);

    // compute the cutting plane
    // remember that we are passing the local vertex positions of the sphere 
    // between the vertex and fragment shader stages, so we need to convert
    // the camera-pos in world coords to sphere's local coords!
    auto c =  model.inverted() * m_camera.get_pos();
    const auto d = c.length();
    const auto r = 1.0f;
    const auto sin_alpha = r / d;
    const auto n = (c / d); // plane unit normal vector
    const auto cos_beta = sin_alpha;
    const auto p = (r * cos_beta) * n;
    QVector4D plane(n.x(), n.y(), n.z(), -QVector3D::dotProduct(p, n));
    const QVector4D arc_color(0, 0.5, 1, 1);
    glLineWidth(5);
    sp.set_uniform("u_color", arc_color);
    sp.set_uniform("u_plane", plane);
    m_geodesic_arcs->draw();

    sp.unuse();
  }
}