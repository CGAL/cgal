/***************************************************************************

    begin                : jan 02
    copyright            : (C) 2002 by Pierre Alliez
    email                : pierre.alliez@sophia.inria.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include <qgl.h>
#include "Camera.h"
#include "Viewport.h"
#include "Arcball.h"

#ifndef _GLVIEW_
#define _GLVIEW_

class QGLView : public QGLWidget
{
    Q_OBJECT
public:
    QGLView(const QGLFormat & format,
            QWidget *parent = 0,
            const char *name = 0,
            int wflags = 0);
    ~QGLView();

	CCamera   theCamera;
	CViewport theViewport;
	CArcball  theArcball;

  enum Rendering_mode{VERTEX, LINE, FILL};

	// rendering options
	bool m_Smooth;
	bool m_Culling;
	bool m_Antialiasing;
	bool m_Lighting;
	bool m_Superimpose_edges;
  bool m_Superimpose_vertices;
	void changeMode(int mode);
	void clearColor(float r,float g,float b);
	void toggleCulling();
  void toggleSmooth();
	void toggleLighting();
	void toggleAntialiasing();
	void toggleSuperimposing();
  void toggleSuperimposeV();
  void toggleMode(Rendering_mode m);
	void enableTexture(bool enable);
	bool m_Texture;
  float m_clear_color[3];
  float m_fore_color[3];
  Rendering_mode m_mode;

public:
  void set_clear_color(float r,float g,float b);
  void set_fore_color(float r,float g,float b);
protected:
	virtual void initializeGL();
	virtual void paintGL();
	virtual void resizeGL( int w, int h );
	virtual GLuint makeObject();
};

#endif // _GLVIEW_

