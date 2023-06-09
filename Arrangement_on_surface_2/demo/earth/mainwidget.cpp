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
out vec3 vNormal;

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
in vec3 vNormal;

out vec4 color;

void main()
{
	const vec3 lightDir = normalize(vec3(1,.5,.5));

	//float c = clamp(dot(lightDir,triNormal), 0, 1);
	vec3 n = normalize(vNormal);
  float c = abs( dot(lightDir, n) );
	color = vec4(.2, .2,0,1) + vec4(c,c,0,1);

	//color = vec4(1,1,0,1);
  //color = vCol;
}
)fs";
}


MainWidget::~MainWidget()
{
    // Make sure the context is current when deleting the texture
    // and the buffers.
    makeCurrent();
    doneCurrent();
}


void MainWidget::mousePressEvent(QMouseEvent *e)
{
    // Save mouse press position
    m_mouse_press_position = QVector2D(e->position());
}
void MainWidget::mouseReleaseEvent(QMouseEvent *e)
{
    // Mouse release position - mouse press position
    QVector2D diff = QVector2D(e->position()) - m_mouse_press_position;

    // Rotation axis is perpendicular to the mouse position difference
    // vector
    QVector3D n = QVector3D(diff.y(), diff.x(), 0.0).normalized();

    // Accelerate angular speed relative to the length of the mouse sweep
    qreal acc = diff.length() / 100.0;

    // Calculate new rotation axis as weighted sum
    m_rotation_axis = (m_rotation_axis * m_angular_speed + n * acc).normalized();

    // Increase angular speed
    m_angular_speed += acc;
}
void MainWidget::timerEvent(QTimerEvent *)
{
    // Decrease angular speed (friction)
    m_angular_speed *= 0.99;

    // Stop rotation when speed goes below threshold
    if (m_angular_speed < 0.01) {
        m_angular_speed = 0.0;
    } else {
        // Update rotation
        m_rotation = QQuaternion::fromAxisAndAngle(m_rotation_axis, m_angular_speed) * 
                                                                       m_rotation;

        // Request an update
        update();
    }

    // redraw every 12 miliseconds
    update();
}



void MainWidget::initializeGL()
{
    initializeOpenGLFunctions();

    glClearColor(0, 0, 0, 1);

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
  shader = glCreateProgram();
  if (!shader)
  {
    std::cout << "error creating shader program!\n";
    return;
  }

  add_shader(shader, vertex_shader_code, GL_VERTEX_SHADER);
  //add_shader(shader, geometry_shader_code, GL_GEOMETRY_SHADER);
  add_shader(shader, fragment_shader_code, GL_FRAGMENT_SHADER);

  GLint result = 0;
  GLchar elog[1024] = { 0 };

  glLinkProgram(shader);
  glGetProgramiv(shader, GL_LINK_STATUS, &result);
  if (!result)
  {
    glGetProgramInfoLog(shader, sizeof(elog), NULL, elog);
    std::cout << "! error linking program:\n" << elog << std::endl;
    return;
  }

  glValidateProgram(shader);
  glGetProgramiv(shader, GL_VALIDATE_STATUS, &result);
  if (!result)
  {
    glGetProgramInfoLog(shader, sizeof(elog), NULL, elog);
    std::cout << "! error validating program:\n" << elog << std::endl;
    return;
  }

  m_uniform_mvp = glGetUniformLocation(shader, "MVP");
  std::cout << "uniform loc = " << m_uniform_mvp << std::endl;
}



void MainWidget::init_geometry()
{
  int num_slices, num_stacks;
  num_slices = num_stacks = 64;
  float r = 3;
  create_sphere(num_slices, num_stacks, r);
}

void MainWidget::create_sphere(int num_slices, int num_stacks, float r)
{
  num_stacks = std::max<int>(2, num_stacks);
  std::vector<QVector3D> vertices, normals;

  // NORTH POLE
  vertices.push_back(QVector3D(0, 0, r));
  normals.push_back(QVector3D(0, 0, 1));

  // SOUTH POLE
  vertices.push_back(QVector3D(0, 0, -r));
  normals.push_back(QVector3D(0, 0, -1));
  int starting_index_of_middle_vertices = vertices.size();

  for (int j = 1; j < num_stacks; ++j)
  {
    // Calculate the latitude (vertical angle) for the current stack
    float lat = M_PI * j / num_stacks;
    float rxy = r * std::sin(lat);
    float z = r * std::cos(lat);

    for (int i = 0; i < num_slices; ++i)
    {
      // Calculate the longitude (horizontal angle) for the current slice
      float lon = 2 * M_PI * i / num_slices;

      // Convert spherical coordinates to Cartesian coordinates
      float x = rxy * std::cos(lon);
      float y = rxy * std::sin(lon);

      auto p = QVector3D(x, y, z);
      auto n = p / p.length();
      vertices.push_back(p);
      normals.push_back(n);
    }
  }

  // strided vertex-data
  std::vector<QVector3D> vertex_data;
  for (int i = 0; i < vertices.size(); ++i)
  {
    vertex_data.push_back(vertices[i]);
    vertex_data.push_back(normals[i]);
  }


  // add the indices for all triangles
  std::vector<GLuint> indices;

  // NORTH CAP
  const int north_vertex_index = 0;
  const int north_cap_vertex_index_start = starting_index_of_middle_vertices;
  for (int i = 0; i < num_slices; i++)
  {
    indices.push_back(north_vertex_index);
    indices.push_back(north_cap_vertex_index_start + i);
    indices.push_back(north_cap_vertex_index_start + (i + 1) % num_slices);
  }

  // 0 = NORTH VERTEX
  // 1 = SOUTH VERTEX
  // [2, 2 + (numSlices-1)] = bottom vertices of the stack #1
  // [2+numSlices, 2 + (2*numSlices - 1)] = bottom vertices of the stack #2
  // ...
  // [2+(k-1)*numSlices, 2 + (k*numSlices -1)] = bottom vertices of the stack #k
  // ..
  // [2+(numStacks-1)*numSlices, 2+(numStacks*numSlices-1)] = bottom vertices of
  //                                                the last stack (# numStacks)

  // SOUTH CAP
  const int south_vertex_index = 1;
  const int south_cap_index_start = starting_index_of_middle_vertices + 
                                                  (num_stacks - 2) * num_slices;
  for (int i = 0; i < num_slices; ++i)
  {
    const auto vi0 = south_vertex_index;
    const auto vi1 = south_cap_index_start + i;
    const auto vi2 = south_cap_index_start + (i + 1) % num_slices;
    indices.push_back(vi2);
    indices.push_back(vi1);
    indices.push_back(vi0);
  }

  // MIDDLE TRIANGLES
  for (int k = 0; k < num_stacks - 2; ++k)
  {
    const int stack_start_index = starting_index_of_middle_vertices + 
                                                                k * num_slices;
    const int next_stack_start_index = stack_start_index + num_slices;
    for (int i = 0; i < num_slices; ++i)
    {
      // check why the following code snippet does not work (winding order?)
      //int vi0 = stackStartIndex + i;
      //int vi1 = nextStackStartIndex + i;
      //int vi2 = nextStackStartIndex + (i + 1) % numSlices;
      //int vi3 = stackStartIndex + (i + 1) % numSlices;
      int vi0 = stack_start_index + i;
      int vi1 = stack_start_index + (i + 1) % num_slices;
      int vi2 = next_stack_start_index + i;
      int vi3 = next_stack_start_index + (i + 1) % num_slices;

      indices.push_back(vi0);
      indices.push_back(vi2);
      indices.push_back(vi1);
      //
      indices.push_back(vi2);
      indices.push_back(vi3);
      indices.push_back(vi1);
    }
  }
  m_num_indices = indices.size();


  // DEFINE OPENGL BUFFERS
  glGenVertexArrays(1, &m_vao);
  glBindVertexArray(m_vao);

  // Index buffer
  glGenBuffers(1, &m_ibo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_ibo);
  auto indices_size = sizeof(GLuint) * indices.size();
  auto indices_data = reinterpret_cast<const void*>(indices.data());
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, 
               indices_size, 
               indices_data,
               GL_STATIC_DRAW);

  // Vertex Buffer
  glGenBuffers(1, &m_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
  auto vertex_buffer_size = sizeof(QVector3D) * vertex_data.size();
  auto vertex_buffer_data = reinterpret_cast<const void*>(vertex_data.data());
  glBufferData(GL_ARRAY_BUFFER, 
               vertex_buffer_size, 
               vertex_buffer_data, 
               GL_STATIC_DRAW);

    // Position Vertex-Attribute
    GLint position_attrib_index = 0;
    const void* position_offset = 0;
    GLsizei stride = 6 * sizeof(float);
    glVertexAttribPointer(position_attrib_index, 
                          3, 
                          GL_FLOAT, GL_FALSE, 
                          stride,
                          position_offset);
    glEnableVertexAttribArray(position_attrib_index);
    
    // Normal Vertex-Attribute
    GLint normal_attrib_index = 1;
    auto* normal_offset = reinterpret_cast<const void*>(3 * sizeof(float));
    glVertexAttribPointer(normal_attrib_index, 
                          3, 
                          GL_FLOAT, 
                          GL_FALSE, 
                          stride, 
                          normal_offset);
    glEnableVertexAttribArray(normal_attrib_index);

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  // Note: calling this before glBindVertexArray(0) results in no output!
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


void MainWidget::resizeGL(int w, int h)
{
  // Calculate aspect ratio
  qreal aspect = qreal(w) / qreal(h ? h : 1);

  // near and far plane locations and vertical field-of-view angle in degrees
  const qreal z_near = 1.0, z_far = 100.0, fov = 45.0;

  // Reset projection
  m_projection.setToIdentity();
  m_projection.perspective(fov, aspect, z_near, z_far);
}

void MainWidget::paintGL()
{
  QMatrix4x4 view;
  const QVector3D eye(0, 10, 10), center(0, 0, 0), up(0, 1, 0);
  view.lookAt(eye, center, up);

  QMatrix4x4 model;
  static float angle = 0;
  angle += 1;
  model.rotate(angle, up);
  
  // Clear color and depth buffer
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  {
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    glUseProgram(shader);
    auto mvp = m_projection * view * model;
    glUniformMatrix4fv(m_uniform_mvp, 1, GL_FALSE, mvp.data());

    {
      // DRAW TRIANGLE
      glBindVertexArray(m_vao);
      {
        //glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        glDrawElements(GL_TRIANGLES, m_num_indices, GL_UNSIGNED_INT, 0);
      }
      glBindVertexArray(0);
    }
    glUseProgram(0);
  }
}
