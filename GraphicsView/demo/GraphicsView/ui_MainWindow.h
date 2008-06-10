/********************************************************************************
** Form generated from reading ui file 'MainWindow.ui'
**
** Created: Tue Jun 10 17:19:28 2008
**      by: Qt User Interface Compiler version 4.3.2
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QGraphicsView>
#include <QtGui/QGridLayout>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>

class Ui_MainWindow
{
public:
    QAction *actionSfs;
    QAction *actionAbout;
    QAction *actionAboutCGAL;
    QAction *actionNew;
    QAction *actionExit;
    QAction *actionInsertRandomPoints;
    QAction *actionMovingPoint;
    QAction *actionInsertPolyline;
    QAction *actionClear;
    QAction *actionShowVoronoi;
    QAction *actionShowDelaunay;
    QAction *actionLoadConstraints;
    QAction *actionSaveConstraints;
    QAction *actionCircumcenter;
    QWidget *centralwidget;
    QGridLayout *gridLayout;
    QGraphicsView *graphicsView;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuEdit;
    QMenu *menuTools;
    QMenu *menuOptions;
    QMenu *menuHelp;
    QStatusBar *statusbar;
    QToolBar *fileToolBar;
    QToolBar *toolBar;

    void setupUi(QMainWindow *MainWindow)
    {
    if (MainWindow->objectName().isEmpty())
        MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
    MainWindow->resize(400, 300);
    MainWindow->setWindowIcon(QIcon(QString::fromUtf8(":/cgal/logos/cgal_logo.xpm")));
    actionSfs = new QAction(MainWindow);
    actionSfs->setObjectName(QString::fromUtf8("actionSfs"));
    actionAbout = new QAction(MainWindow);
    actionAbout->setObjectName(QString::fromUtf8("actionAbout"));
    actionAboutCGAL = new QAction(MainWindow);
    actionAboutCGAL->setObjectName(QString::fromUtf8("actionAboutCGAL"));
    actionNew = new QAction(MainWindow);
    actionNew->setObjectName(QString::fromUtf8("actionNew"));
    actionExit = new QAction(MainWindow);
    actionExit->setObjectName(QString::fromUtf8("actionExit"));
    actionInsertRandomPoints = new QAction(MainWindow);
    actionInsertRandomPoints->setObjectName(QString::fromUtf8("actionInsertRandomPoints"));
    actionMovingPoint = new QAction(MainWindow);
    actionMovingPoint->setObjectName(QString::fromUtf8("actionMovingPoint"));
    actionMovingPoint->setCheckable(true);
    actionMovingPoint->setIcon(QIcon(QString::fromUtf8(":/cgal/Actions/icons/movingPoint.png")));
    actionInsertPolyline = new QAction(MainWindow);
    actionInsertPolyline->setObjectName(QString::fromUtf8("actionInsertPolyline"));
    actionInsertPolyline->setCheckable(true);
    actionInsertPolyline->setChecked(false);
    actionInsertPolyline->setIcon(QIcon(QString::fromUtf8(":/cgal/Input/inputPolyline.png")));
    actionClear = new QAction(MainWindow);
    actionClear->setObjectName(QString::fromUtf8("actionClear"));
    actionClear->setIcon(QIcon(QString::fromUtf8(":/cgal/fileToolbar/fileNew.png")));
    actionShowVoronoi = new QAction(MainWindow);
    actionShowVoronoi->setObjectName(QString::fromUtf8("actionShowVoronoi"));
    actionShowVoronoi->setCheckable(true);
    actionShowVoronoi->setChecked(false);
    actionShowVoronoi->setIcon(QIcon(QString::fromUtf8(":/cgal/Triangulation_2/Voronoi_diagram_2.png")));
    actionShowDelaunay = new QAction(MainWindow);
    actionShowDelaunay->setObjectName(QString::fromUtf8("actionShowDelaunay"));
    actionShowDelaunay->setCheckable(true);
    actionShowDelaunay->setIcon(QIcon(QString::fromUtf8(":/cgal/Triangulation_2/Delaunay_triangulation_2.png")));
    actionLoadConstraints = new QAction(MainWindow);
    actionLoadConstraints->setObjectName(QString::fromUtf8("actionLoadConstraints"));
    actionLoadConstraints->setIcon(QIcon(QString::fromUtf8(":/cgal/fileToolbar/fileOpen.png")));
    actionSaveConstraints = new QAction(MainWindow);
    actionSaveConstraints->setObjectName(QString::fromUtf8("actionSaveConstraints"));
    actionSaveConstraints->setIcon(QIcon(QString::fromUtf8(":/cgal/fileToolbar/fileSave.png")));
    actionCircumcenter = new QAction(MainWindow);
    actionCircumcenter->setObjectName(QString::fromUtf8("actionCircumcenter"));
    actionCircumcenter->setCheckable(true);
    actionCircumcenter->setIcon(QIcon(QString::fromUtf8(":/cgal/Actions/icons/circumcenter.png")));
    centralwidget = new QWidget(MainWindow);
    centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
    gridLayout = new QGridLayout(centralwidget);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
    graphicsView = new QGraphicsView(centralwidget);
    graphicsView->setObjectName(QString::fromUtf8("graphicsView"));

    gridLayout->addWidget(graphicsView, 0, 0, 1, 1);

    MainWindow->setCentralWidget(centralwidget);
    menubar = new QMenuBar(MainWindow);
    menubar->setObjectName(QString::fromUtf8("menubar"));
    menubar->setGeometry(QRect(0, 0, 400, 19));
    menuFile = new QMenu(menubar);
    menuFile->setObjectName(QString::fromUtf8("menuFile"));
    menuEdit = new QMenu(menubar);
    menuEdit->setObjectName(QString::fromUtf8("menuEdit"));
    menuTools = new QMenu(menubar);
    menuTools->setObjectName(QString::fromUtf8("menuTools"));
    menuOptions = new QMenu(menubar);
    menuOptions->setObjectName(QString::fromUtf8("menuOptions"));
    menuHelp = new QMenu(menubar);
    menuHelp->setObjectName(QString::fromUtf8("menuHelp"));
    MainWindow->setMenuBar(menubar);
    statusbar = new QStatusBar(MainWindow);
    statusbar->setObjectName(QString::fromUtf8("statusbar"));
    MainWindow->setStatusBar(statusbar);
    fileToolBar = new QToolBar(MainWindow);
    fileToolBar->setObjectName(QString::fromUtf8("fileToolBar"));
    MainWindow->addToolBar(Qt::TopToolBarArea, fileToolBar);
    toolBar = new QToolBar(MainWindow);
    toolBar->setObjectName(QString::fromUtf8("toolBar"));
    MainWindow->addToolBar(Qt::TopToolBarArea, toolBar);

    menubar->addAction(menuFile->menuAction());
    menubar->addAction(menuEdit->menuAction());
    menubar->addAction(menuTools->menuAction());
    menubar->addAction(menuOptions->menuAction());
    menubar->addAction(menuHelp->menuAction());
    menuFile->addSeparator();
    menuFile->addAction(actionClear);
    menuFile->addAction(actionLoadConstraints);
    menuFile->addAction(actionSaveConstraints);
    menuFile->addSeparator();
    menuFile->addAction(actionExit);
    menuEdit->addAction(actionInsertRandomPoints);
    menuTools->addAction(actionMovingPoint);
    menuTools->addAction(actionInsertPolyline);
    menuTools->addAction(actionShowDelaunay);
    menuTools->addAction(actionShowVoronoi);
    menuHelp->addAction(actionAbout);
    menuHelp->addAction(actionAboutCGAL);
    fileToolBar->addAction(actionClear);
    fileToolBar->addAction(actionLoadConstraints);
    fileToolBar->addAction(actionSaveConstraints);
    toolBar->addAction(actionInsertPolyline);
    toolBar->addAction(actionMovingPoint);
    toolBar->addAction(actionShowDelaunay);
    toolBar->addAction(actionShowVoronoi);
    toolBar->addAction(actionCircumcenter);

    retranslateUi(MainWindow);

    QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
    MainWindow->setWindowTitle(QApplication::translate("MainWindow", "CGAL Constrained Delaunay Triangulation", 0, QApplication::UnicodeUTF8));
    actionSfs->setText(QApplication::translate("MainWindow", "sfs", 0, QApplication::UnicodeUTF8));
    actionAbout->setText(QApplication::translate("MainWindow", "About", 0, QApplication::UnicodeUTF8));
    actionAboutCGAL->setText(QApplication::translate("MainWindow", "About CGAL", 0, QApplication::UnicodeUTF8));
    actionNew->setText(QApplication::translate("MainWindow", "New", 0, QApplication::UnicodeUTF8));
    actionExit->setText(QApplication::translate("MainWindow", "Quit", 0, QApplication::UnicodeUTF8));
    actionInsertRandomPoints->setText(QApplication::translate("MainWindow", "Insert random points", 0, QApplication::UnicodeUTF8));
    actionMovingPoint->setText(QApplication::translate("MainWindow", "Moving Point", 0, QApplication::UnicodeUTF8));
    actionMovingPoint->setToolTip(QApplication::translate("MainWindow", "Simulate Insertion", "The comment", QApplication::UnicodeUTF8));
    actionMovingPoint->setStatusTip(QApplication::translate("MainWindow", "Move mouse with left button pressed to see where point would be inserted", "and its comment", QApplication::UnicodeUTF8));
    actionMovingPoint->setWhatsThis(QApplication::translate("MainWindow", "whats this", "what", QApplication::UnicodeUTF8));
    actionInsertPolyline->setText(QApplication::translate("MainWindow", "Insert Polyline", 0, QApplication::UnicodeUTF8));
    actionInsertPolyline->setToolTip(QApplication::translate("MainWindow", "Insert Point or Polyline", 0, QApplication::UnicodeUTF8));
    actionInsertPolyline->setStatusTip(QApplication::translate("MainWindow", "Left: Insert vtx | Right: Final vtx | Del: Delete vtx", 0, QApplication::UnicodeUTF8));
    actionClear->setText(QApplication::translate("MainWindow", "Clear", 0, QApplication::UnicodeUTF8));
    actionShowVoronoi->setText(QApplication::translate("MainWindow", "Show Voronoi Diagram", 0, QApplication::UnicodeUTF8));
    actionShowDelaunay->setText(QApplication::translate("MainWindow", "Show Delaunay Triangulation", 0, QApplication::UnicodeUTF8));
    actionLoadConstraints->setText(QApplication::translate("MainWindow", "Load Constraints", 0, QApplication::UnicodeUTF8));
    actionSaveConstraints->setText(QApplication::translate("MainWindow", "Save Constraints", 0, QApplication::UnicodeUTF8));
    actionCircumcenter->setText(QApplication::translate("MainWindow", "Circumcenter", 0, QApplication::UnicodeUTF8));
    actionCircumcenter->setToolTip(QApplication::translate("MainWindow", "Draw circumcenter", 0, QApplication::UnicodeUTF8));
    menuFile->setTitle(QApplication::translate("MainWindow", "File", 0, QApplication::UnicodeUTF8));
    menuEdit->setTitle(QApplication::translate("MainWindow", "Edit", 0, QApplication::UnicodeUTF8));
    menuTools->setTitle(QApplication::translate("MainWindow", "Tools", 0, QApplication::UnicodeUTF8));
    menuOptions->setTitle(QApplication::translate("MainWindow", "Options", 0, QApplication::UnicodeUTF8));
    menuHelp->setTitle(QApplication::translate("MainWindow", "Help", 0, QApplication::UnicodeUTF8));
    fileToolBar->setWindowTitle(QApplication::translate("MainWindow", "File Tools", 0, QApplication::UnicodeUTF8));
    toolBar->setWindowTitle(QApplication::translate("MainWindow", "Visualization Tools", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

#endif // UI_MAINWINDOW_H
