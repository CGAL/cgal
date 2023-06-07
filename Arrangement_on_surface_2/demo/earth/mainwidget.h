// Copyright (C) 2016 The Qt Company Ltd.
// SPDX-License-Identifier: LicenseRef-Qt-Commercial OR BSD-3-Clause

#ifndef MAINWIDGET_H
#define MAINWIDGET_H


#include <QOpenGLWidget>
#include <QMatrix4x4>
#include <QQuaternion>
#include <QVector2D>
#include <QBasicTimer>

#include <qopenglwidget.h>
#include <qopenglfunctions_3_3_core.h>
#include <qopenglfunctions_4_5_core.h>


class GeometryEngine;

class MainWidget : public QOpenGLWidget, protected QOpenGLFunctions_3_3_Core
{
    Q_OBJECT

public:
    using QOpenGLWidget::QOpenGLWidget;
    ~MainWidget();

protected:
    void mousePressEvent(QMouseEvent *e) override;
    void mouseReleaseEvent(QMouseEvent *e) override;
    void timerEvent(QTimerEvent *e) override;

    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;


    GLuint vao, vbo, ibo, shader, numIndices;
    GLuint uniformMVP; // ModelViewProjection
    void addShader(GLuint theProgram, const char* shaderCode, GLenum shaderType);
    void initShaderProgram();
    
    void initGeometry();
    void createSphere(int numSlices, int numStacks, float r);


private:
    QBasicTimer timer;

    QMatrix4x4 projection;

    QVector2D mousePressPosition;
    QVector3D rotationAxis;
    qreal angularSpeed = 0;
    QQuaternion rotation;
};

#endif // MAINWIDGET_H
