// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_QT_WIDGET_OPENGL_H
#define CGAL_QT_WIDGET_OPENGL_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Nef_S2/OGL_base_object.h>
#include <CGAL/Nef_S2/Sphere_geometry_OGL.h>

#include <qapplication.h>
#include <qgl.h>
#include <qwidget.h>
#include <qevent.h>
#include <qmenu.h>
#include <qpoint.h>

namespace CGAL {

    class Qt_widget_OpenGL : public QGLWidget {
        Q_OBJECT

            typedef CGAL::OGL::OGL_base_object::Affine_3       Affine_3;
        typedef CGAL::OGL::OGL_base_object::Double_vector  Double_vector;
        //    typedef CGAL::Simple_cartesian<double>::Aff_transformation_3  Affine_3;

        protected:
        enum MenuEntriesS2 { ROTATE, SCALE, TRANSLATE, TRANS_Z, RESET_CONTROL, 
            PERSP, FULLSCREEN, QUIT };

        size_t object_index_;
        std::vector<CGAL::OGL::OGL_base_object*> object_;
        std::vector<std::string> label_;

        int window_width;          // Breite und
        int window_height;         // Hoehe des Fensters
        int window_radius;           // min(width,height) / 2

        bool perspective;
        bool fullscreen;

        int mouse_x, mouse_y;      // Mauskoordinaten linker button
        int interaction;                   // type of interaction in motion fct.
        int motion_mode;           // Bewegen der Maus bei Mouse1 gedrueckt
        int submenu1, submenu2;
        long double dx;                // Translation
        long double dy;                // Translation
        long double dz;                     // Translation in Z
        long double s;                 // Skalierung
        long double init_s;

        Affine_3 rotation;   // Rotation

        long double factor_s;              // Umrechnungsfaktor fuer Skalierung
        long double factor_w;
        long double factor_d;

        QMenu* main;
        QMenu* sub1;
        QMenu* sub2;
        QMenu* sub3;

        protected:
        Affine_3 virtual_sphere_transformation( double old_x, double old_y, 
                double new_x, double new_y);
        void controlMenu(int index);
        void renderMenu(int index);
        void optionsMenu(int index);

        protected slots:
        void slotControlMenuReset() {controlMenu(RESET_CONTROL);}
        void slotControlMenuRotate() {controlMenu(ROTATE);}
        void slotControlMenuScale() {controlMenu(SCALE);}
        void slotControlMenuTranslate() {controlMenu(TRANSLATE);}
        void slotRenderMenuFaces() {renderMenu(CGAL::OGL::SM_FACES);}
        void slotRenderMenuSkeleton() {renderMenu(CGAL::OGL::SM_SKELETON);}
        void slotRenderMenuTriangulation() {renderMenu(CGAL::OGL::SM_TRIANGULATION);}
        void slotOptionsMenuAxes() {optionsMenu(CGAL::OGL::SM_AXES);}
        void slotOptionsMenuUnityCube() {optionsMenu(CGAL::OGL::SM_CUBE);}
        void slotFullscreen();
        void slotPerspective();
        void slotNextObject();
        void slotPrevObject();

        protected:
        virtual void mouseMoveEvent( QMouseEvent* event);
        virtual void mousePressEvent( QMouseEvent* event);
        virtual void mouseReleaseEvent( QMouseEvent* event);

        public:
        Qt_widget_OpenGL(int width, int height, double scale);
        virtual ~Qt_widget_OpenGL();
        virtual void paintGL();
        virtual void initializeGL();
        virtual void resizeGL(int w, int h);
    };

} // namespace CGAL
#endif // CGAL_QT_WIDGET_OPENGL_H
