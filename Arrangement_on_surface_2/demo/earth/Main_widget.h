// Copyright (C) 2016 The Qt Company Ltd.
// SPDX-License-Identifier: LicenseRef-Qt-Commercial OR BSD-3-Clause

#ifndef MAIN_WIDGET_H
#define MAIN_WIDGET_H

#include <QOpenGLWidget>
#include <QMatrix4x4>
#include <QQuaternion>
#include <QVector2D>
#include <QBasicTimer>

#include <functional>
#include <memory>

#include <qopenglwidget.h>

#include "Camera.h"
#include "Camera_manip.h"
#include "Common_defs.h"
#include "Kml_reader.h"
#include "Line_strips.h"
#include "Shader_program.h"
#include "Sphere.h"
#include "Vertices.h"
#include "World_coordinate_axes.h"


class Main_widget : public QOpenGLWidget, protected OpenGLFunctionsBase
{
  Q_OBJECT

public:
  using QOpenGLWidget::QOpenGLWidget;
  ~Main_widget();

protected:
  void mousePressEvent(QMouseEvent* e) override;
  void mouseMoveEvent(QMouseEvent* e) override;
  void mouseReleaseEvent(QMouseEvent* e) override;
  void timerEvent(QTimerEvent* e) override;
  void keyPressEvent(QKeyEvent* event) override;

  void initializeGL() override;
  void resizeGL(int w, int h) override;
  void paintGL() override;

    
  void init_camera();
  void init_geometry();  
  void init_shader_programs();

  void init_country_borders(float error);
  void init_country_selection();
  

  float compute_backprojected_error(float pixel_error);
  // Use this to find the approximate of the true minimum projected error.
  // we are ot using this complicated method, but provide it for completeness.
  void find_minimum_projected_error_on_sphere(float we);

  // verify that the node (180.0, -84.71338) in Antarctica is redundant
  void verify_antarctica_node_is_redundant();

  // init problematic vertices: these are the vertices incident to deg-4 vertex
  void init_problematic_nodes();

private:
  // Objects in the scene
  std::unique_ptr<Sphere>           m_sphere;
  std::unique_ptr<World_coord_axes> m_world_coord_axes;
  std::unique_ptr<Line_strips>      m_geodesic_arcs;
  std::unique_ptr<Vertices>         m_vertices, m_problematic_vertices;
  std::unique_ptr<Line_strips>      m_identification_curve;

  // COUNTRY DATA
  Kml::Placemarks                             m_countries;
  std::vector<std::string>                    m_country_names;
  std::vector<std::unique_ptr<Line_strips>>   m_country_borders;

  // now we draw boundary-arcs by country
  int             m_selected_country_index, m_selected_arc_index;
  Kml::Nodes      m_selected_country_nodes;
  Kml::Arcs       m_selected_country_arcs;
  Kml::Placemark* m_selected_country;
 

  // Shaders
  Shader_program  m_sp_smooth;
  Shader_program  m_sp_per_vertex_color;
  Shader_program  m_sp_arc;
  
  // Camera & controls
  Camera  m_camera;
  std::unique_ptr<Camera_manip>  m_camera_manip_rot;
  std::unique_ptr<Camera_manip>  m_camera_manip_zoom;

  // view-port 
  int m_vp_width = 0, m_vp_height = 0;

  // Timer for continuous screen-updates
  QBasicTimer m_timer;
};

#endif
