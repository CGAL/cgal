// Copyright (C) 2016 The Qt Company Ltd.
// SPDX-License-Identifier: LicenseRef-Qt-Commercial OR BSD-3-Clause

#include "mainwidget.h"

#include <QMouseEvent>

#include <cmath>
#include <iostream>
#include <string>
using namespace std;

namespace {
  // vertex shader
  const char* vShader = R"vs(
#version 330

layout (location = 0) in vec3 pos;
//out vec4 vCol;
out vec3 vpos;

uniform mat4 MVP; 

void main()
{
	vpos = pos;
	gl_Position = MVP * vec4(pos.xyz, 1);
  //gl_Position = vec4(pos.xyz, 1);
}
)vs";


  // GEOMETRY SHADER
  // * I am using the geometry shader to compute the face-normals in the GPU on the fly
  const char* gShader = R"gs(
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
  static const char* fShader = R"fs(
#version 330

in vec4 vCol;
out vec4 color;

void main()
{
	//color = vec4(1,1,0,1);
  color = vCol;
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
    mousePressPosition = QVector2D(e->position());
}
void MainWidget::mouseReleaseEvent(QMouseEvent *e)
{
    // Mouse release position - mouse press position
    QVector2D diff = QVector2D(e->position()) - mousePressPosition;

    // Rotation axis is perpendicular to the mouse position difference
    // vector
    QVector3D n = QVector3D(diff.y(), diff.x(), 0.0).normalized();

    // Accelerate angular speed relative to the length of the mouse sweep
    qreal acc = diff.length() / 100.0;

    // Calculate new rotation axis as weighted sum
    rotationAxis = (rotationAxis * angularSpeed + n * acc).normalized();

    // Increase angular speed
    angularSpeed += acc;
}
void MainWidget::timerEvent(QTimerEvent *)
{
    // Decrease angular speed (friction)
    angularSpeed *= 0.99;

    // Stop rotation when speed goes below threshold
    if (angularSpeed < 0.01) {
        angularSpeed = 0.0;
    } else {
        // Update rotation
        rotation = QQuaternion::fromAxisAndAngle(rotationAxis, angularSpeed) * rotation;

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

    //initGeometry();
    createSphere(20, 10, 3);
    initShaderProgram();

    // Enable depth buffer
    glEnable(GL_DEPTH_TEST);

    // Enable back face culling
    //glEnable(GL_CULL_FACE);

    // Use QBasicTimer because its faster than QTimer
    timer.start(12, this);
}

void MainWidget::addShader(GLuint theProgram, const char* shaderCode, GLenum shaderType)
{
  GLuint theShader = glCreateShader(shaderType);

  const GLchar* theCode[] = { shaderCode };
  GLint codeLength[] = { strlen(shaderCode) };

  glShaderSource(theShader, 1, theCode, codeLength);
  glCompileShader(theShader);


  GLint result = 0;
  GLchar elog[1024] = { 0 };
  glGetShaderiv(theShader, GL_COMPILE_STATUS, &result);
  if (!result)
  {
    string shaderTypeName;
    switch (shaderType)
    {
    case GL_VERTEX_SHADER:   shaderTypeName = "VERTEX"; break;
    case GL_GEOMETRY_SHADER: shaderTypeName = "GEOMETRY"; break;
    case GL_FRAGMENT_SHADER: shaderTypeName = "FRAGMENT"; break;
    }
    glGetShaderInfoLog(theShader, sizeof(elog), NULL, elog);
    cout << "! error compiling the " << shaderTypeName << " shader:\n" << elog << endl;
    return;
  }

  glAttachShader(theProgram, theShader);
}
void MainWidget::initShaderProgram()
{
  shader = glCreateProgram();
  if (!shader)
  {
    cout << "error creating shader program!\n";
    return;
  }

  addShader(shader, vShader, GL_VERTEX_SHADER);
  addShader(shader, gShader, GL_GEOMETRY_SHADER);
  addShader(shader, fShader, GL_FRAGMENT_SHADER);

  GLint result = 0;
  GLchar elog[1024] = { 0 };

  glLinkProgram(shader);
  glGetProgramiv(shader, GL_LINK_STATUS, &result);
  if (!result)
  {
    glGetProgramInfoLog(shader, sizeof(elog), NULL, elog);
    cout << "! error linking program:\n" << elog << endl;
    return;
  }

  glValidateProgram(shader);
  glGetProgramiv(shader, GL_VALIDATE_STATUS, &result);
  if (!result)
  {
    glGetProgramInfoLog(shader, sizeof(elog), NULL, elog);
    cout << "! error validating program:\n" << elog << endl;
    return;
  }

  uniformMVP = glGetUniformLocation(shader, "MVP");
  cout << "uniform loc = " << uniformMVP << endl;
}



void MainWidget::initGeometry()
{
  const float c = 0.5;
  GLfloat vertices[] = {
    -c, -c,  0,
     c, -c,  0,
     0,  c,  0
  };

  GLuint indices[] = { 0,1,2 };

  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  {
    glGenBuffers(1, &ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    {
      glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

      GLint index = 0;
      glVertexAttribPointer(index, 3, GL_FLOAT, GL_FALSE, 0, 0);
      glEnableVertexAttribArray(index);
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);

  }
  glBindVertexArray(0);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}
void MainWidget::createSphere(int numSlices, int numStacks, float r)
{
  numStacks = std::max<int>(2, numStacks);
  std::vector<QVector3D> vertices, normals;

  // NORTH POLE
  vertices.push_back(QVector3D(0, 0, r));
  normals.push_back(QVector3D(0, 0, 1));

  // SOUTH POLE
  vertices.push_back(QVector3D(0, 0, -r));
  normals.push_back(QVector3D(0, 0, -1));
  int startingIndexOfMiddleVertices = vertices.size();

  for (int j = 1; j < numStacks; ++j)
  {
    // Calculate the latitude (vertical angle) for the current stack
    float lat = M_PI * j / numStacks;
    float rxy = r * sin(lat);
    float z = r * cos(lat);

    for (int i = 0; i < numSlices; ++i)
    {
      // Calculate the longitude (horizontal angle) for the current slice
      float lon = 2 * M_PI * i / numSlices;

      // Convert spherical coordinates to Cartesian coordinates
      float x = rxy * cos(lon);
      float y = rxy * sin(lon);

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
  const int northVertexIndex = 0;
  const int northCapVertexIndexStart = startingIndexOfMiddleVertices;
  for (int i = 0; i < numSlices; i++)
  {
    indices.push_back(northVertexIndex);
    indices.push_back(northCapVertexIndexStart + i);
    indices.push_back(northCapVertexIndexStart + (i + 1) % numSlices);
  }

  // 0 = NORTH VERTEX
  // 1 = SOUTH VERTEX
  // [2, 2 + (numSlices-1)] = bottom vertices of the stack #1
  // [2+numSlices, 2 + (2*numSlices - 1)] = bottom vertices of the stack #2
  // ...
  // [2+(k-1)*numSlices, 2 + (k*numSlices -1) ] = bottom vertices of the stack #k
  // ..
  // [2+(numStacks-1)*numSlices, 2+(numStacks*numSlices-1)] = bottom vertices of the last stack (# numStacks)

  // SOUTH CAP
  const int southVertexIndex = 1;
  const int southCapIndexStart = startingIndexOfMiddleVertices + (numStacks - 2) * numSlices;
  for (int i = 0; i < numSlices; i++)
  {
    const auto vi0 = southVertexIndex;
    const auto vi1 = southCapIndexStart + i;
    const auto vi2 = southCapIndexStart + (i + 1) % numSlices;
    indices.push_back(vi2);
    indices.push_back(vi1);
    indices.push_back(vi0);
  }

  // MIDDLE TRIANGLES
  for (int k = 0; k < numStacks - 2; k++)
  {
    const int stackStartIndex = startingIndexOfMiddleVertices + k * numSlices;
    const int nextStackStartIndex = stackStartIndex + numSlices;
    for (int i = 0; i < numSlices; i++)
    {
      //int vi0 = stackStartIndex + i;
      //int vi1 = nextStackStartIndex + i;
      //int vi2 = nextStackStartIndex + (i + 1) % numSlices;
      //int vi3 = stackStartIndex + (i + 1) % numSlices;
      int vi0 = stackStartIndex + i;
      int vi1 = stackStartIndex + (i + 1) % numSlices;
      int vi2 = nextStackStartIndex + i;
      int vi3 = nextStackStartIndex + (i + 1) % numSlices;

      indices.push_back(vi0);
      indices.push_back(vi2);
      indices.push_back(vi1);
      //
      indices.push_back(vi2);
      indices.push_back(vi3);
      indices.push_back(vi1);
    }
  }

  numIndices = indices.size();


  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  {
    glGenBuffers(1, &ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * indices.size(), reinterpret_cast<const void*>(indices.data()), GL_STATIC_DRAW);

    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    {
      glBufferData(GL_ARRAY_BUFFER, sizeof(QVector3D) * vertex_data.size(), reinterpret_cast<const void*>(vertex_data.data()), GL_STATIC_DRAW);

      // positions
      GLint positionAttribIndex = 0;
      GLsizei stride = 6 * sizeof(float);
      glVertexAttribPointer(positionAttribIndex, 3, GL_FLOAT, GL_FALSE, stride, 0);
      glEnableVertexAttribArray(positionAttribIndex);
      //normals
      GLint normalAttribIndex = 1;
      auto* normal_offset = reinterpret_cast<const void*>(3 * sizeof(float));
      glVertexAttribPointer(normalAttribIndex, 3, GL_FLOAT, GL_FALSE, stride, normal_offset);
      glEnableVertexAttribArray(normalAttribIndex);
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);

  }
  glBindVertexArray(0);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}
//! [3]


void MainWidget::resizeGL(int w, int h)
{
    // Calculate aspect ratio
    qreal aspect = qreal(w) / qreal(h ? h : 1);

    // near and far plane locations and vertical field-of-view angle in degrees
    const qreal z_near = 1.0, z_far = 100.0, fov = 45.0;

    // Reset projection
    projection.setToIdentity();
    projection.perspective(fov, aspect, z_near, z_far);
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
    auto mvp = projection * view * model;
    glUniformMatrix4fv(uniformMVP, 1, GL_FALSE, mvp.data());

    {
      // DRAW TRIANGLE
      glBindVertexArray(vao);
      {
        //glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
        glDrawElements(GL_TRIANGLES, numIndices, GL_UNSIGNED_INT, 0);
      }
      glBindVertexArray(0);
    }
    glUseProgram(0);
  }
}
