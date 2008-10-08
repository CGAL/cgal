// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// $URL
// $Id 
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado, 
//             Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#include <qgl.h>
#include <CGAL/IO/Qt_widget_OpenGL.h>
#include <CGAL/IO/Qt_widget.h>

// INCLUDES DE QT
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qinputdialog.h>

// INCLUDES DES IMAGES

#include "images/button_solid.xpm"
#include "images/button_wire.xpm"
#include "images/arrow_01_down.xpm"
#include "images/arrow_01_up.xpm"
#include "images/arrow_01_left.xpm"
#include "images/arrow_01_right.xpm"
#include "images/zoom_in.xpm"
#include "images/zoom_out.xpm"

#include <cmath>
#include <cstdio>
#include <cstdlib>
//#include <unistd.h> // Header File For sleeping.
#include <vector>

// INCLUDES D'OPENGL
//#include <GL/glut.h> // Header File For The GLUT Library
#include <CGAL/gl.h> // Header File For The OpenGL32 Library
#include <CGAL/glu.h> // Header File For The GLu32 Library

#include <CGAL/Cartesian.h>
#include <CGAL/Algebraic_kernel_for_spheres_2_3.h>
#include <CGAL/Spherical_kernel_3.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/IO/Spherical_circle_gl.h>

// DEFINITION DES TYPES
typedef double                                           NT;
typedef CGAL::Cartesian<NT>                              Linear_k1;
typedef CGAL::Algebraic_kernel_for_spheres_2_3<NT>       Algebraic_k1;
typedef CGAL::Spherical_kernel_3<Linear_k1,Algebraic_k1> SK;
typedef SK::Plane_3                                      Plane_3;
typedef SK::Point_3                                      Point_3;
typedef SK::Sphere_3                                     Sphere_3;
typedef SK::Circle_3                                     Circle_3;

class MyGLDrawer : public QGLWidget {

	Q_OBJECT

	public:

		MyGLDrawer (QWidget* parent, const char* name)
			: QGLWidget(parent, name), showContour(false) {
				this->x = 0.f;
				this->y = 0.f;
				this->z = -6.f;
				this->rx = 0.f;
				this->ry = 0.f;
				this->rz = 0.f;
				this->list_cercle = 0;
			}

		void paintGL () {
			//std::cout << "MyGLDrawer.paintGL();" << std::endl;
			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear The Screen And The Depth Buffer
			glLoadIdentity(); // Reset The View
			glTranslatef(this->x, this->y, this->z);
			glRotatef(this->rx, 1.0f, 0.0f, 0.0f);
			glRotatef(this->ry, 0.0f, 1.0f, 0.0f);
			glRotatef(this->rz, 0.0f, 0.0f, 1.0f);
			glCallList(this->list_cercle);
			this->swapBuffers();
		}

		void initializeGL () {
			glClearColor(1.0f, 1.0f, 1.0f, 0.0f); // This Will Clear The Background Color To Black
			glClearDepth(1.0); // Enables Clearing Of The Depth Buffer
			glDepthFunc(GL_LESS); // The Type Of Depth Test To Do
			glEnable(GL_DEPTH_TEST); // Enables Depth Testing
			glShadeModel(GL_SMOOTH); // Enables Smooth Color Shading
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity(); // Reset The Projection Matrix
			glMatrixMode(GL_MODELVIEW);

			// creation de la liste openGL
			if (this->list_cercle != 0) glDeleteLists(this->list_cercle, 1);
			this->list_cercle = glGenLists(1);
			//On cree une liste dans laquelle on met l'objet 3ds
			//On cree une liste dans laquelle on met l'objet 3ds
			glNewList(this->list_cercle, GL_COMPILE);
			for (std::size_t i = 0; i < this->cercles.size(); i++) {
			  CGAL::dessiner_spherical_circle<SK>(this->cercles[i].first, this->cercles[i].second);
			}
			glEndList();

		}

		void resizeGL (int w, int h) {
			glViewport(0, 0, w, h); // Reset The Current Viewport And Perspective Transformation
			glMatrixMode(GL_PROJECTION);
			glLoadIdentity();
			gluPerspective(45.0f, (GLfloat) w / (GLfloat) w, 0.1f, 100.0f);
			glMatrixMode(GL_MODELVIEW);
		}

		inline void add_cercle (const Circle_3& c, int precision) {
			this->cercles.push_back(std::pair<Circle_3, int>(c, precision));
		}

	public slots:

		void wire_display () {
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			this->paintGL();
		}

		void surface_display () {
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			this->paintGL();
		}

		void translate(float x, float y, float z){
			this->x = x;
			this->y = y;
			this->z = z;
		}

		void rotate(float rx, float ry, float rz){
			this->rx = rx;
			this->ry = ry;
			this->rz = rz;
		}

		inline void show_contour () {this->showContour = true; paintGL();}
		inline void hide_contour () {this->showContour = false; paintGL();}

	private:
		float x, y, z, rx, ry, rz;
		int list_cercle;
		bool showContour;
		std::vector<std::pair<Circle_3, int> > cercles;

};


class MyWindow : public QMainWindow {
	Q_OBJECT

	public:

		MyWindow (int w, int h) {
			//initialization of coordonnate of camera
			this->x = 0.f;
			this->y = 0.f;
			this->z = -6.f;
			this->rx = 0.f;
			this->ry = 0.f;
			this->rz = 0.f;

			this->resize(w, h);
			wGL = new MyGLDrawer(this, "Demo Circle 3D");
			this->setCentralWidget(wGL);
			this->setCaption("Demo Circle 3D");
			//File Menu
			QPopupMenu* file = new QPopupMenu(this);
			file->insertItem("&Quit", qApp, SLOT(closeAllWindows()),'');
			//Display Menu
			QPopupMenu * affichage = new QPopupMenu(this);
			affichage->insertItem("&Wire Mode",this,SLOT(mode_wire()), 'r');
			affichage->insertItem("&Surface Mode",this,SLOT(mode_surface()),'t');
			affichage->insertItem("Show &Border",this,SLOT(show_contour()),'f');
			affichage->insertItem("&Hide Border",this,SLOT(hide_contour()),'g');
			//Move Menu
			QPopupMenu * move = new QPopupMenu(this);
			move->insertItem("Rotate &Left",this,SLOT(rotate_left()), 4116);
			move->insertItem("Rotate &Right",this,SLOT(rotate_right()), 4114);
			move->insertItem("Rotate Z Left",this,SLOT(rotate_tonneau_left()),'q');
			move->insertItem("Rotate Z Right",this,SLOT(rotate_tonneau_right()), 'e');
			move->insertItem("Rotate &Up",this,SLOT(rotate_up()),4115);
			move->insertItem("Rotate &Down",this,SLOT(rotate_down()), 4117);
			move->insertItem("Zoom &In",this,SLOT(zoom_in()), 4118);
			move->insertItem("Zoom &Out",this,SLOT(zoom_out()), 4119);
			//Menu bar
			this->menuBar()->insertItem("&File", file);
			this->menuBar()->insertItem("&Display", affichage);
			this->menuBar()->insertItem("&Move", move);
			//tool bar
			this->layers_toolbar = new QToolBar("Tools", this,
					QMainWindow::Top, TRUE, "Tools");
			//buttons in tool bar
			this->show_wire_button =
				new QToolButton(QPixmap((const char**)::button_wire_xpm),
						"Wire Mode",
						0,
						this,
						SLOT(mode_wire()),
						this->layers_toolbar,
						"Wire Mode");
			this->show_wire_button->setToggleButton(true);
			this->show_surface_button =
				new QToolButton(QPixmap((const char**)::button_solid_xpm),
						"Surface Mode",
						0,
						this,
						SLOT(mode_surface()),
						this->layers_toolbar,
						"Mode Surface");
			this->show_surface_button->setToggleButton(true);
			this->show_surface_button->toggle();
			this->layers_toolbar->addSeparator();
			this->rotate_left_button =
				new QToolButton(QPixmap((const char**)::arrow_01_left_xpm),
						"Rotate Left",
						0,
						this,
						SLOT(rotate_left()),
						this->layers_toolbar,
						"Rotate Left");
			this->rotate_right_button =
				new QToolButton(QPixmap((const char**)::arrow_01_right_xpm),
						"Rotate Right",
						0,
						this,
						SLOT(rotate_right()),
						layers_toolbar,
						"Rotate Right");

			this->rotate_up_button =
				new QToolButton(QPixmap((const char**)::arrow_01_up_xpm),
						"Rotate Up",
						0,
						this,
						SLOT(rotate_up()),
						this->layers_toolbar,
						"Rotate UP");
			this->rotate_down_button =
				new QToolButton(QPixmap((const char**)::arrow_01_down_xpm),
						"Rotate Down",
						0,
						this,
						SLOT(rotate_down()),
						this->layers_toolbar,
						"Rotate Down");
			this->zoom_in_button =
				new QToolButton(QPixmap((const char**)::zoom_in_xpm),
						"Zoom In",
						0,
						this,
						SLOT(zoom_in()),
						this->layers_toolbar,
						"Zoom In");
			this->zoom_out_button =
				new QToolButton(QPixmap((const char**)::zoom_out_xpm),
						"Zoom Out",
						0,
						this,
						SLOT(zoom_out()),
						this->layers_toolbar,
						"Zoom Out");
		}

	protected:

		void keyPressEvent (QKeyEvent* /*qke*/) {
			//std::cout << "Key ascii<" << qke->ascii() << "> key<" << qke->key() << ">" << std::endl;
		}

	public slots:

		inline void show_contour() {this->wGL->show_contour();}
		inline void hide_contour() {this->wGL->hide_contour();}

		void mode_wire () {
			std::cout << "Mode Wire" << std::endl;
			this->show_wire_button->setOn(true);
			this->show_surface_button->setOn(false);
			this->wGL->wire_display();
		}

		void mode_surface () {
			std::cout << "Mode Surface" << std::endl;
			this->show_wire_button->setOn(false);
			this->show_surface_button->setOn(true);
			this->wGL->surface_display();
		}

		inline void add_cercle (const Circle_3& p, int i) {
			this->wGL->add_cercle(p,i);
		}

		void rotate_tonneau_right () {
			//std::cout << "Rotation en tonneau a droite" << std::endl;
			this->rz -= 1.f;
			this->wGL->rotate(this->rx, this->ry, this->rz);
			this->wGL->paintGL();
		}
		void rotate_tonneau_left () {
			//std::cout << "Rotation en tonneau a gauche" << std::endl;
			this->rz += 1.f;
			this->wGL->rotate(this->rx, this->ry, this->rz);
			this->wGL->paintGL();
		}
		void rotate_right(){
			//std::cout << "Rotation a droite" << std::endl;
			this->ry += 1.f;
			this->wGL->rotate(this->rx, this->ry, this->rz);
			this->wGL->paintGL();
		}
		void rotate_left(){
			//std::cout << "Rotation a gauche" << std::endl;
			this->ry -= 1.f;
			this->wGL->rotate(this->rx, this->ry, this->rz);
			this->wGL->paintGL();
		}
		void rotate_up(){
			//std::cout << "Rotation vers le haut" << std::endl;
			this->rx -= 1.f;
			this->wGL->rotate(this->rx, this->ry, this->rz);
			this->wGL->paintGL();
		}
		void rotate_down(){
			//std::cout << "Rotation vers le bas" << std::endl;
			this->rx += 1.f;
			this->wGL->rotate(this->rx, this->ry, this->rz);
			this->wGL->paintGL();
		}
		void zoom_in(){
			//std::cout << "Augmentation du zoom" << std::endl;
			this->z += 0.1f;
			this->wGL->translate(this->x, this->y, this->z);
			this->wGL->paintGL();
		}
		void zoom_out(){
			//std::cout << "Diminution du zoom" << std::endl;
			this->z -= 0.1f;
			this->wGL->translate(this->x, this->y, this->z);
			this->wGL->paintGL();
		}

	private:
		MyGLDrawer* wGL;
		QToolBar* layers_toolbar;
		QToolButton* show_wire_button;
		QToolButton* show_surface_button;
		QToolButton* rotate_right_button;
		QToolButton* rotate_left_button;
		QToolButton* rotate_up_button;
		QToolButton* rotate_down_button;
		QToolButton* zoom_in_button;
		QToolButton* zoom_out_button;
		float x, y, z, rx, ry, rz; //coordinate of camera

};
