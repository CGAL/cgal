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
//out vec3 vpos;

//uniform mat4 MVP; 

void main()
{
	//vpos = pos;
	//gl_Position = MVP * vec4(pos.xyz, 1);
  gl_Position = vec4(pos.xyz, 1);
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
	float c = clamp(dot(lightDir,triNormal), 0, 1);
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

//in vec4 vCol;
out vec4 color;

void main()
{
	color = vec4(1,1,0,1);
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

//! [0]
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
//! [0]

//! [1]
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
}
//! [1]



void MainWidget::initializeGL()
{
    initializeOpenGLFunctions();

    glClearColor(0, 0, 0, 1);

    initGeometry();
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
  //addShader(shader, gShader, GL_GEOMETRY_SHADER);
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
//! [3]


//! [5]
void MainWidget::resizeGL(int w, int h)
{
    // Calculate aspect ratio
    qreal aspect = qreal(w) / qreal(h ? h : 1);

    // Set near plane to 3.0, far plane to 7.0, field of view 45 degrees
    const qreal zNear = 3.0, zFar = 7.0, fov = 45.0;

    // Reset projection
    projection.setToIdentity();

    // Set perspective projection
    projection.perspective(fov, aspect, zNear, zFar);
}
//! [5]

void MainWidget::paintGL()
{
    // Clear color and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    {
      glClearColor(0, 0, 0, 1);
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


      glUseProgram(shader);
      {
        // DRAW TRIANGLE
        glBindVertexArray(vao);
        {
          //glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
          glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, 0);
        }
        glBindVertexArray(0);
      }
      glUseProgram(0);
    }
}
