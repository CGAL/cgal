#include "CGAL/Nef_S2/Qt_widget_Nef_S2.h"

namespace CGAL {
    Qt_widget_Nef_S2::~Qt_widget_Nef_S2() {
        for (size_t i=0;i<object_.size();i++) {
            delete object_[i];
        }
        object_.clear();
    }

    Qt_widget_Nef_S2::Qt_widget_Nef_S2() :  
        Qt_widget_OpenGL(300,300,1.5) {

            resize(window_width, window_height);
            // this->setContextMenuPolicy(Qt::ActionsContextMenu);

            main = new QMenu;
            sub1 = main->addMenu("&Control");
            sub2 = main->addMenu("&Render");
            sub3 = main->addMenu("&Options");

            QAction * qa;
            qa = sub1->addAction("Reset");
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotControlMenuReset()));
            qa = sub1->addAction("Rotate");
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotControlMenuRotate()));
            qa = sub1->addAction("Scale");
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotControlMenuScale()));
            qa = sub1->addAction("Translate in XY");
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotControlMenuTranslate()));

            qa = sub2->addAction("Faces");
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotRenderMenuFaces()));
            qa = sub2->addAction("Skeleton");
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotRenderMenuSkeleton()));
            qa = sub2->addAction("Triangulation");
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotRenderMenuTriangulation()));

            qa = sub3->addAction("Toggle Axes");
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotOptionsMenuAxes()));
            qa = sub3->addAction("Toggle Unity Cube");
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotOptionsMenuUnityCube()));

            qa = main->addAction("&Next");
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotNextObject()));
            qa = main->addAction("&Previous");
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotPrevObject()));
            qa = main->addAction("&Toggle Fullscreen");
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotFullscreen()));
            qa = main->addAction("&Quit");
            QObject::connect(qa, SIGNAL(triggered()), qApp, SLOT(quit()));

            qa = new QAction("Reset",this);
            qa->setShortcut(Qt::Key_R);
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotControlMenuReset()));
            this->addAction(qa);

            qa = new QAction("Toggle Axes",this);
            qa->setShortcut(Qt::Key_A);
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotOptionsMenuAxes()));
            this->addAction(qa);

            qa = new QAction("Toggle Unity Cube",this);
            qa->setShortcut(Qt::Key_C);
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotOptionsMenuUnityCube()));
            this->addAction(qa);

            qa = new QAction("Next",this);
            qa->setShortcut(Qt::Key_Right);
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotNextObject()));
            this->addAction(qa);

            qa = new QAction("Previous",this);
            qa->setShortcut(Qt::Key_Left);
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotPrevObject()));
            this->addAction(qa);

            qa = new QAction("Toggle Fullscreen",this);
            qa->setShortcut(Qt::Key_F);
            QObject::connect(qa, SIGNAL(triggered()), this, SLOT(slotFullscreen()));
            this->addAction(qa);

            qa = new QAction("Quit",this);
            qa->setShortcut(Qt::Key_Q);
            QObject::connect(qa, SIGNAL(triggered()), qApp, SLOT(quit()));
            this->addAction(qa);
        }


} //namespace CGAL

