// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008 Benoit Jacob <jacob.benoit.1@gmail.com>
//
// Eigen is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// Eigen is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#ifndef MANDELBROT_H
#define MANDELBROT_H

#include <Eigen/Core>
#include <QtGui/QApplication>
#include <QtGui/QWidget>
#include <QtCore/QThread>

class MandelbrotWidget;

class MandelbrotThread : public QThread
{
    friend class MandelbrotWidget;
    MandelbrotWidget *widget;
    long long total_iter;
    int id, max_iter;
    bool single_precision;

  public:
    MandelbrotThread(MandelbrotWidget *w, int i) : widget(w), id(i) {}
    void run();
    template<typename Real> void render(int img_width, int img_height);
};

class MandelbrotWidget : public QWidget
{
    Q_OBJECT

    friend class MandelbrotThread;
    Eigen::Vector2d center;
    double xradius;
    int size;
    unsigned char *buffer;
    QPoint lastpos;
    int draft;
    MandelbrotThread **threads;
    int threadcount;

  protected:
    void resizeEvent(QResizeEvent *);
    void paintEvent(QPaintEvent *);
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);

  public:
    MandelbrotWidget() : QWidget(), center(0,0), xradius(2),
                         size(0), buffer(0), draft(16)
    {
      setAutoFillBackground(false);
      threadcount = QThread::idealThreadCount();
      threads = new MandelbrotThread*[threadcount];
      for(int th = 0; th < threadcount; th++) threads[th] = new MandelbrotThread(this, th);
    }
    ~MandelbrotWidget()
    {
      if(buffer) delete[]buffer;
      for(int th = 0; th < threadcount; th++) delete threads[th];
      delete[] threads;
    }
};

#endif // MANDELBROT_H
