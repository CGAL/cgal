// Copyright (C) 2016 The Qt Company Ltd.
// SPDX-License-Identifier: LicenseRef-Qt-Commercial OR BSD-3-Clause

#include "mainwidget.h"

#include <QMouseEvent>

#include <cmath>
#include <iostream>
#include <string>


namespace {
  // vertex shader
  const char* vertex_shader_code = R"vs(
#version 330

layout (location = 0) in vec3 pos;
layout (location = 1) in vec3 normal;

//out vec4 vCol;
//out vec3 vpos;
flat out vec3 vNormal;

uniform mat4 MVP; 

void main()
{
	//vpos = pos;
  vNormal = normal;
	gl_Position = MVP * vec4(pos.xyz, 1);
  //gl_Position = vec4(pos.xyz, 1);
}
)vs";


  // GEOMETRY SHADER
  // * I am using the geometry shader to compute the face-normals in the GPU on the fly
  const char* geometry_shader_code = R"gs(
#version 330

in vec3 vpos[];
out vec4 vCol;

layout (triangles) in;
layout (triangle_strip, max_vertices = 3) out;

void main()
{ 
	const vec3 lightDir = normalize(vec3(1,.5,.5));

	// compute the normal for the current triangle
	vec3 triNormal = normalize(cross(vpos[1]-vpos[0], vpos[2]-vpos[0]));
	//float c = clamp(dot(lightDir,triNormal), 0, 1);
	float c = abs(dot(lightDir,triNormal));
	vCol = vec4(.2, .2,0,1) + vec4(c,c,0,1);

	gl_Position = gl_in[0].gl_Position; EmitVertex();
	gl_Position = gl_in[1].gl_Position; EmitVertex();
	gl_Position = gl_in[2].gl_Position; EmitVertex();
	EndPrimitive();
}
)gs";


  // FRAGMENT SHADER
  static const char* fragment_shader_code = R"fs(
#version 330

//in vec4 vCol;
flat in vec3 vNormal;

out vec4 color;

void main()
{
	const vec3 lightDir = normalize(vec3(1,.5,.5));

	//float c = clamp(dot(lightDir,triNormal), 0, 1);
	vec3 n = normalize(vNormal);
  float c = abs( dot(lightDir, n) );
	color = vec4(.2, .2,0,1) + 0.8*vec4(c,c,0,1);

	//color = vec4(1,1,0,1);
  //color = vCol;
}
)fs";
}


MainWidget::~MainWidget()
{
  // Make sure the context is current when deleting the texture and the buffers.
  makeCurrent();
  doneCurrent();
}


void MainWidget::mousePressEvent(QMouseEvent *e)
{
  m_mouse_pressed = true;
  m_last_mouse_pos = QVector2D(e->position());
}




void MainWidget::mouseMoveEvent(QMouseEvent* e)
{
  auto current_mouse_pos = QVector2D(e->position());

  if (m_mouse_pressed)
  {  
    const auto diff = current_mouse_pos - m_last_mouse_pos;
    const float scale_factor = 0.1f;
    const float theta_around_x = scale_factor * diff.y();
    const float theta_around_y = scale_factor * diff.x();
  
    m_camera.rotate(theta_around_x, theta_around_y);
  }

  m_last_mouse_pos = current_mouse_pos;
}
void MainWidget::mouseReleaseEvent(QMouseEvent *e)
{
  m_mouse_pressed = false;

}
void MainWidget::timerEvent(QTimerEvent *)
{
  update();
}



void MainWidget::initializeGL()
{
    initializeOpenGLFunctions();

    glClearColor(0, 0, 0, 1);

    m_camera.set_pos(0,0,10);
    init_geometry();
    init_shader_program();

    // Enable depth buffer
    glEnable(GL_DEPTH_TEST);

    // Enable back face culling
    //glEnable(GL_CULL_FACE);

    // Use QBasicTimer because its faster than QTimer
    m_timer.start(12, this);
}

void MainWidget::add_shader(GLuint the_program, const char* shader_code, 
                                                            GLenum shader_type)
{
  GLuint the_shader = glCreateShader(shader_type);

  const GLchar* the_code[] = { shader_code };
  GLint code_length[] = { strlen(shader_code) };

  glShaderSource(the_shader, 1, the_code, code_length);
  glCompileShader(the_shader);


  GLint result = 0;
  GLchar elog[1024] = { 0 };
  glGetShaderiv(the_shader, GL_COMPILE_STATUS, &result);
  if (!result)
  {
    std::string shader_type_name;
    switch (shader_type)
    {
    case GL_VERTEX_SHADER:   shader_type_name = "VERTEX"; break;
    case GL_GEOMETRY_SHADER: shader_type_name = "GEOMETRY"; break;
    case GL_FRAGMENT_SHADER: shader_type_name = "FRAGMENT"; break;
    }
    glGetShaderInfoLog(the_shader, sizeof(elog), NULL, elog);
    std::cout << "! error compiling the " << shader_type_name << 
                 " shader:\n" << elog << std::endl;
    return;
  }

  glAttachShader(the_program, the_shader);
}
void MainWidget::init_shader_program()
{
  m_shader = glCreateProgram();
  if (!m_shader)
  {
    std::cout << "error creating shader program!\n";
    return;
  }

  add_shader(m_shader, vertex_shader_code, GL_VERTEX_SHADER);
  //add_shader(shader, geometry_shader_code, GL_GEOMETRY_SHADER);
  add_shader(m_shader, fragment_shader_code, GL_FRAGMENT_SHADER);

  GLint result = 0;
  GLchar elog[1024] = { 0 };

  glLinkProgram(m_shader);
  glGetProgramiv(m_shader, GL_LINK_STATUS, &result);
  if (!result)
  {
    glGetProgramInfoLog(m_shader, sizeof(elog), NULL, elog);
    std::cout << "! error linking program:\n" << elog << std::endl;
    return;
  }

  glValidateProgram(m_shader);
  glGetProgramiv(m_shader, GL_VALIDATE_STATUS, &result);
  if (!result)
  {
    glGetProgramInfoLog(m_shader, sizeof(elog), NULL, elog);
    std::cout << "! error validating program:\n" << elog << std::endl;
    return;
  }

  m_uniform_mvp = glGetUniformLocation(m_shader, "MVP");
  std::cout << "uniform loc = " << m_uniform_mvp << std::endl;
}



void MainWidget::init_geometry()
{
  int num_slices, num_stacks;
  num_slices = num_stacks = 64;
  float r = 3;
  m_sphere = std::make_unique<Sphere>(num_slices, num_stacks, r);
}



void MainWidget::resizeGL(int w, int h)
{
  // Calculate aspect ratio
  qreal aspect = qreal(w) / qreal(h ? h : 1);

  // near and far plane locations and vertical field-of-view angle in degrees
  const qreal z_near = 1.0, z_far = 100.0, fov = 45.0;

  // Reset projection
  m_camera.perspective(fov, aspect, z_near, z_far);
}


void MainWidget::paintGL()
{
  QMatrix4x4 model;
  const auto view = m_camera.get_view_matrix();
  const auto projection = m_camera.get_projection_matrix();
  const auto mvp = projection * view * model;
  

  // Clear color and depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  {
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUseProgram(m_shader);
    glUniformMatrix4fv(m_uniform_mvp, 1, GL_FALSE, mvp.data());
    
    m_sphere->draw();

    glUseProgram(0);
  }
}


