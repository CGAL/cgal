/********************************************************************************
** Form generated from reading ui file 'MainWindow.ui'
**
** Created: Mon Dec 20 14:14:32 2010
**      by: Qt User Interface Compiler version 4.4.1
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHBoxLayout>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>
#include "Viewer.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionGenerate_Points;
    QAction *actionLoad_Points;
    QAction *actionSave_Points;
    QAction *actionShow_Axis;
    QAction *actionQuit;
    QAction *actionClear_Scene;
    QAction *actionShow_Vertex;
    QAction *actionShow_DEdge;
    QAction *actionShow_VEdge;
    QAction *actionShow_Facet;
    QAction *actionFlat;
    QAction *actionPreferences;
    QAction *actionInsert_Vertex;
    QAction *actionInsert_Point;
    QAction *actionSelect_Vertex;
    QAction *actionMove_Vertex;
    QAction *actionFind_NearestNb;
    QAction *actionEmpty_Sphere;
    QAction *actionNormal_View;
    QAction *actionDemo_Help;
    QAction *actionIncremental_Construct;
    QAction *actionStop_Animation;
    QAction *actionAbout_T3_demo;
    QWidget *centralwidget;
    QHBoxLayout *horizontalLayout;
    Viewer *viewer;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuEdit;
    QMenu *menuMode;
    QMenu *menuShow;
    QMenu *menuHelp;
    QStatusBar *statusbar;
    QToolBar *toolBar;

    void setupUi(QMainWindow *MainWindow)
    {
    if (MainWindow->objectName().isEmpty())
        MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
    MainWindow->setWindowModality(Qt::NonModal);
    MainWindow->resize(1100, 500);
    QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(MainWindow->sizePolicy().hasHeightForWidth());
    MainWindow->setSizePolicy(sizePolicy);
    QFont font;
    font.setFamily(QString::fromUtf8("Arial"));
    MainWindow->setFont(font);
    MainWindow->setCursor(QCursor(Qt::PointingHandCursor));
    QIcon icon;
    icon.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/cgal_logo.xpm")), QIcon::Normal, QIcon::Off);
    MainWindow->setWindowIcon(icon);
    actionGenerate_Points = new QAction(MainWindow);
    actionGenerate_Points->setObjectName(QString::fromUtf8("actionGenerate_Points"));
    QIcon icon1;
    icon1.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/pointRandom.png")), QIcon::Normal, QIcon::Off);
    actionGenerate_Points->setIcon(icon1);
    actionLoad_Points = new QAction(MainWindow);
    actionLoad_Points->setObjectName(QString::fromUtf8("actionLoad_Points"));
    QIcon icon2;
    icon2.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/fileOpen.png")), QIcon::Normal, QIcon::Off);
    actionLoad_Points->setIcon(icon2);
    actionSave_Points = new QAction(MainWindow);
    actionSave_Points->setObjectName(QString::fromUtf8("actionSave_Points"));
    QIcon icon3;
    icon3.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/fileSave.png")), QIcon::Normal, QIcon::Off);
    actionSave_Points->setIcon(icon3);
    actionShow_Axis = new QAction(MainWindow);
    actionShow_Axis->setObjectName(QString::fromUtf8("actionShow_Axis"));
    actionShow_Axis->setCheckable(true);
    QIcon icon4;
    icon4.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/coordinates.jpeg")), QIcon::Normal, QIcon::Off);
    actionShow_Axis->setIcon(icon4);
    actionQuit = new QAction(MainWindow);
    actionQuit->setObjectName(QString::fromUtf8("actionQuit"));
    QIcon icon5;
    icon5.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/quit.jpeg")), QIcon::Normal, QIcon::Off);
    actionQuit->setIcon(icon5);
    actionClear_Scene = new QAction(MainWindow);
    actionClear_Scene->setObjectName(QString::fromUtf8("actionClear_Scene"));
    QIcon icon6;
    icon6.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/clear.jpeg")), QIcon::Normal, QIcon::Off);
    actionClear_Scene->setIcon(icon6);
    actionShow_Vertex = new QAction(MainWindow);
    actionShow_Vertex->setObjectName(QString::fromUtf8("actionShow_Vertex"));
    actionShow_Vertex->setCheckable(true);
    actionShow_Vertex->setChecked(true);
    QIcon icon7;
    icon7.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/show_point.jpeg")), QIcon::Normal, QIcon::Off);
    actionShow_Vertex->setIcon(icon7);
    actionShow_DEdge = new QAction(MainWindow);
    actionShow_DEdge->setObjectName(QString::fromUtf8("actionShow_DEdge"));
    actionShow_DEdge->setCheckable(true);
    actionShow_DEdge->setChecked(true);
    QIcon icon8;
    icon8.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/show_delaunay.jpeg")), QIcon::Normal, QIcon::Off);
    actionShow_DEdge->setIcon(icon8);
    actionShow_VEdge = new QAction(MainWindow);
    actionShow_VEdge->setObjectName(QString::fromUtf8("actionShow_VEdge"));
    actionShow_VEdge->setCheckable(true);
    QIcon icon9;
    icon9.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/show_voronoi.jpeg")), QIcon::Normal, QIcon::Off);
    actionShow_VEdge->setIcon(icon9);
    actionShow_Facet = new QAction(MainWindow);
    actionShow_Facet->setObjectName(QString::fromUtf8("actionShow_Facet"));
    actionShow_Facet->setCheckable(true);
    QIcon icon10;
    icon10.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/show_facet.jpeg")), QIcon::Normal, QIcon::Off);
    actionShow_Facet->setIcon(icon10);
    actionFlat = new QAction(MainWindow);
    actionFlat->setObjectName(QString::fromUtf8("actionFlat"));
    actionFlat->setCheckable(true);
    QIcon icon11;
    icon11.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/flat.png")), QIcon::Normal, QIcon::Off);
    icon11.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/stereo.png")), QIcon::Normal, QIcon::On);
    actionFlat->setIcon(icon11);
    actionPreferences = new QAction(MainWindow);
    actionPreferences->setObjectName(QString::fromUtf8("actionPreferences"));
    QIcon icon12;
    icon12.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/preferences.jpeg")), QIcon::Normal, QIcon::Off);
    actionPreferences->setIcon(icon12);
    actionInsert_Vertex = new QAction(MainWindow);
    actionInsert_Vertex->setObjectName(QString::fromUtf8("actionInsert_Vertex"));
    actionInsert_Vertex->setCheckable(true);
    QIcon icon13;
    icon13.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/insert.jpeg")), QIcon::Normal, QIcon::Off);
    actionInsert_Vertex->setIcon(icon13);
    actionInsert_Point = new QAction(MainWindow);
    actionInsert_Point->setObjectName(QString::fromUtf8("actionInsert_Point"));
    actionInsert_Point->setCheckable(true);
    QIcon icon14;
    icon14.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/insert_point.jpg")), QIcon::Normal, QIcon::Off);
    actionInsert_Point->setIcon(icon14);
    actionSelect_Vertex = new QAction(MainWindow);
    actionSelect_Vertex->setObjectName(QString::fromUtf8("actionSelect_Vertex"));
    actionSelect_Vertex->setCheckable(true);
    QIcon icon15;
    icon15.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/select_hand.jpeg")), QIcon::Normal, QIcon::Off);
    actionSelect_Vertex->setIcon(icon15);
    actionMove_Vertex = new QAction(MainWindow);
    actionMove_Vertex->setObjectName(QString::fromUtf8("actionMove_Vertex"));
    actionMove_Vertex->setCheckable(true);
    QIcon icon16;
    icon16.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/move_1.jpeg")), QIcon::Normal, QIcon::Off);
    actionMove_Vertex->setIcon(icon16);
    actionFind_NearestNb = new QAction(MainWindow);
    actionFind_NearestNb->setObjectName(QString::fromUtf8("actionFind_NearestNb"));
    actionFind_NearestNb->setCheckable(true);
    QIcon icon17;
    icon17.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/nearest_nb.png")), QIcon::Normal, QIcon::Off);
    actionFind_NearestNb->setIcon(icon17);
    actionEmpty_Sphere = new QAction(MainWindow);
    actionEmpty_Sphere->setObjectName(QString::fromUtf8("actionEmpty_Sphere"));
    actionEmpty_Sphere->setCheckable(true);
    QIcon icon18;
    icon18.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/empty_sphere.jpeg")), QIcon::Normal, QIcon::Off);
    actionEmpty_Sphere->setIcon(icon18);
    actionNormal_View = new QAction(MainWindow);
    actionNormal_View->setObjectName(QString::fromUtf8("actionNormal_View"));
    actionNormal_View->setCheckable(true);
    actionNormal_View->setChecked(true);
    QIcon icon19;
    icon19.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/normal_view.jpeg")), QIcon::Normal, QIcon::Off);
    actionNormal_View->setIcon(icon19);
    actionDemo_Help = new QAction(MainWindow);
    actionDemo_Help->setObjectName(QString::fromUtf8("actionDemo_Help"));
    actionIncremental_Construct = new QAction(MainWindow);
    actionIncremental_Construct->setObjectName(QString::fromUtf8("actionIncremental_Construct"));
    actionIncremental_Construct->setCheckable(true);
    QIcon icon20;
    icon20.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/play.jpeg")), QIcon::Normal, QIcon::Off);
    icon20.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/pause.jpeg")), QIcon::Normal, QIcon::On);
    actionIncremental_Construct->setIcon(icon20);
    actionStop_Animation = new QAction(MainWindow);
    actionStop_Animation->setObjectName(QString::fromUtf8("actionStop_Animation"));
    QIcon icon21;
    icon21.addPixmap(QPixmap(QString::fromUtf8(":/T3_demo/icons/icons/stop.jpeg")), QIcon::Normal, QIcon::Off);
    actionStop_Animation->setIcon(icon21);
    actionAbout_T3_demo = new QAction(MainWindow);
    actionAbout_T3_demo->setObjectName(QString::fromUtf8("actionAbout_T3_demo"));
    centralwidget = new QWidget(MainWindow);
    centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
    horizontalLayout = new QHBoxLayout(centralwidget);
    horizontalLayout->setSpacing(3);
    horizontalLayout->setMargin(1);
    horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
    viewer = new Viewer(centralwidget);
    viewer->setObjectName(QString::fromUtf8("viewer"));
    sizePolicy.setHeightForWidth(viewer->sizePolicy().hasHeightForWidth());
    viewer->setSizePolicy(sizePolicy);
    viewer->setCursor(QCursor(Qt::PointingHandCursor));

    horizontalLayout->addWidget(viewer);

    MainWindow->setCentralWidget(centralwidget);
    menubar = new QMenuBar(MainWindow);
    menubar->setObjectName(QString::fromUtf8("menubar"));
    menubar->setGeometry(QRect(0, 0, 1010, 22));
    menuFile = new QMenu(menubar);
    menuFile->setObjectName(QString::fromUtf8("menuFile"));
    menuEdit = new QMenu(menubar);
    menuEdit->setObjectName(QString::fromUtf8("menuEdit"));
    menuMode = new QMenu(menubar);
    menuMode->setObjectName(QString::fromUtf8("menuMode"));
    menuShow = new QMenu(menubar);
    menuShow->setObjectName(QString::fromUtf8("menuShow"));
    menuHelp = new QMenu(menubar);
    menuHelp->setObjectName(QString::fromUtf8("menuHelp"));
    MainWindow->setMenuBar(menubar);
    statusbar = new QStatusBar(MainWindow);
    statusbar->setObjectName(QString::fromUtf8("statusbar"));
    MainWindow->setStatusBar(statusbar);
    toolBar = new QToolBar(MainWindow);
    toolBar->setObjectName(QString::fromUtf8("toolBar"));
    MainWindow->addToolBar(Qt::TopToolBarArea, toolBar);

    menubar->addAction(menuFile->menuAction());
    menubar->addAction(menuEdit->menuAction());
    menubar->addAction(menuMode->menuAction());
    menubar->addAction(menuShow->menuAction());
    menubar->addAction(menuHelp->menuAction());
    menuFile->addAction(actionLoad_Points);
    menuFile->addAction(actionSave_Points);
    menuFile->addSeparator();
    menuFile->addAction(actionQuit);
    menuEdit->addAction(actionGenerate_Points);
    menuEdit->addSeparator();
    menuEdit->addAction(actionIncremental_Construct);
    menuEdit->addAction(actionStop_Animation);
    menuEdit->addSeparator();
    menuEdit->addAction(actionClear_Scene);
    menuMode->addAction(actionNormal_View);
    menuMode->addAction(actionInsert_Vertex);
    menuMode->addAction(actionInsert_Point);
    menuMode->addAction(actionSelect_Vertex);
    menuMode->addAction(actionMove_Vertex);
    menuMode->addAction(actionFind_NearestNb);
    menuMode->addAction(actionEmpty_Sphere);
    menuShow->addAction(actionShow_Axis);
    menuShow->addSeparator();
    menuShow->addAction(actionShow_Vertex);
    menuShow->addAction(actionShow_DEdge);
    menuShow->addAction(actionShow_VEdge);
    menuShow->addAction(actionShow_Facet);
    menuShow->addSeparator();
    menuShow->addAction(actionFlat);
    menuShow->addSeparator();
    menuShow->addAction(actionPreferences);
    menuHelp->addAction(actionDemo_Help);
    menuHelp->addAction(actionAbout_T3_demo);
    toolBar->addAction(actionLoad_Points);
    toolBar->addAction(actionSave_Points);
    toolBar->addSeparator();
    toolBar->addAction(actionGenerate_Points);
    toolBar->addAction(actionClear_Scene);
    toolBar->addSeparator();
    toolBar->addAction(actionIncremental_Construct);
    toolBar->addAction(actionStop_Animation);
    toolBar->addSeparator();
    toolBar->addAction(actionShow_Axis);
    toolBar->addSeparator();
    toolBar->addAction(actionFlat);
    toolBar->addSeparator();
    toolBar->addAction(actionShow_Vertex);
    toolBar->addAction(actionShow_DEdge);
    toolBar->addAction(actionShow_VEdge);
    toolBar->addAction(actionShow_Facet);
    toolBar->addSeparator();
    toolBar->addAction(actionNormal_View);
    toolBar->addAction(actionInsert_Vertex);
    toolBar->addAction(actionInsert_Point);
    toolBar->addAction(actionSelect_Vertex);
    toolBar->addAction(actionMove_Vertex);
    toolBar->addAction(actionFind_NearestNb);
    toolBar->addAction(actionEmpty_Sphere);
    toolBar->addSeparator();
    toolBar->addAction(actionPreferences);
    toolBar->addSeparator();
    toolBar->addAction(actionQuit);

    retranslateUi(MainWindow);

    QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
    MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Triangulation_demo_3", 0, QApplication::UnicodeUTF8));
    actionGenerate_Points->setText(QApplication::translate("MainWindow", "Generate Points", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionGenerate_Points->setToolTip(QApplication::translate("MainWindow", "Generate Points", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionGenerate_Points->setStatusTip(QApplication::translate("MainWindow", "Generate Points", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionGenerate_Points->setShortcut(QApplication::translate("MainWindow", "Ctrl+G", 0, QApplication::UnicodeUTF8));
    actionLoad_Points->setText(QApplication::translate("MainWindow", "Load Points...", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionLoad_Points->setToolTip(QApplication::translate("MainWindow", "Load Points from a file", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionLoad_Points->setStatusTip(QApplication::translate("MainWindow", "Load Points from a file", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionLoad_Points->setShortcut(QApplication::translate("MainWindow", "Ctrl+O", 0, QApplication::UnicodeUTF8));
    actionSave_Points->setText(QApplication::translate("MainWindow", "Save Points...", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionSave_Points->setToolTip(QApplication::translate("MainWindow", "Save points to a file", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionSave_Points->setStatusTip(QApplication::translate("MainWindow", "Save points to a file", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionSave_Points->setShortcut(QApplication::translate("MainWindow", "Ctrl+S", 0, QApplication::UnicodeUTF8));
    actionShow_Axis->setText(QApplication::translate("MainWindow", "Show Axis", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionShow_Axis->setToolTip(QApplication::translate("MainWindow", "Show/Hide Axis", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionShow_Axis->setStatusTip(QApplication::translate("MainWindow", "Show/Hide Axis", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionQuit->setText(QApplication::translate("MainWindow", "Quit", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionQuit->setToolTip(QApplication::translate("MainWindow", "Quit", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionQuit->setStatusTip(QApplication::translate("MainWindow", "Quit", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionClear_Scene->setText(QApplication::translate("MainWindow", "Clear Scene", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionClear_Scene->setToolTip(QApplication::translate("MainWindow", "Clear Scene", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionClear_Scene->setStatusTip(QApplication::translate("MainWindow", "Clear Scene", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionShow_Vertex->setText(QApplication::translate("MainWindow", "Show Vertices", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionShow_Vertex->setToolTip(QApplication::translate("MainWindow", "Show Vertices", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionShow_Vertex->setStatusTip(QApplication::translate("MainWindow", "Show Vertices", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionShow_DEdge->setText(QApplication::translate("MainWindow", "Show Delaunay Edges", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionShow_DEdge->setToolTip(QApplication::translate("MainWindow", "Show Delaunay edges", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionShow_DEdge->setStatusTip(QApplication::translate("MainWindow", "Show Delaunay edges", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionShow_VEdge->setText(QApplication::translate("MainWindow", "Show Voronoi Edges", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionShow_VEdge->setToolTip(QApplication::translate("MainWindow", "Show Voronoi edges", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionShow_VEdge->setStatusTip(QApplication::translate("MainWindow", "Show Voronoi edges", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionShow_Facet->setText(QApplication::translate("MainWindow", "Show Facets", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionShow_Facet->setToolTip(QApplication::translate("MainWindow", "Show Delaunay Facets", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionShow_Facet->setStatusTip(QApplication::translate("MainWindow", "Show Delaunay Facets", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionFlat->setText(QApplication::translate("MainWindow", "Flat", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionFlat->setToolTip(QApplication::translate("MainWindow", "Toggle 3D effect", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionFlat->setStatusTip(QApplication::translate("MainWindow", "Toggle 3D effect", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionPreferences->setText(QApplication::translate("MainWindow", "Preferences...", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionPreferences->setToolTip(QApplication::translate("MainWindow", "Change Colors, Transparency, etc.", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionPreferences->setStatusTip(QApplication::translate("MainWindow", "Change Colors, Transparency, etc.", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionInsert_Vertex->setText(QApplication::translate("MainWindow", "Insert Vertex", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionInsert_Vertex->setToolTip(QApplication::translate("MainWindow", "Insert a vertex and update the triangulation", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionInsert_Vertex->setStatusTip(QApplication::translate("MainWindow", "Insert a vertex and update the triangulation", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionInsert_Point->setText(QApplication::translate("MainWindow", "Insert Point", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionInsert_Point->setToolTip(QApplication::translate("MainWindow", "Insert a point and show its conflict region", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionInsert_Point->setStatusTip(QApplication::translate("MainWindow", "Insert a point and show its conflict region", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionSelect_Vertex->setText(QApplication::translate("MainWindow", "Select Vertex", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionSelect_Vertex->setToolTip(QApplication::translate("MainWindow", "Select vertices", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionSelect_Vertex->setStatusTip(QApplication::translate("MainWindow", "Select vertices", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionMove_Vertex->setText(QApplication::translate("MainWindow", "Move Vertex", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionMove_Vertex->setToolTip(QApplication::translate("MainWindow", "Move a vertex", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionMove_Vertex->setStatusTip(QApplication::translate("MainWindow", "Move a vertex", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionFind_NearestNb->setText(QApplication::translate("MainWindow", "Nearest Neighbor Search", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionFind_NearestNb->setToolTip(QApplication::translate("MainWindow", "Find the nearest neighbor of the query point", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionFind_NearestNb->setStatusTip(QApplication::translate("MainWindow", "Find the nearest neighbor of the query point", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionEmpty_Sphere->setText(QApplication::translate("MainWindow", "Show Empty Sphere", 0, QApplication::UnicodeUTF8));
    actionEmpty_Sphere->setIconText(QApplication::translate("MainWindow", "Click to select a cell and show empty sphere of that cell.", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionEmpty_Sphere->setToolTip(QApplication::translate("MainWindow", "Locate the query point in a cell and show empty sphere of that cell", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionEmpty_Sphere->setStatusTip(QApplication::translate("MainWindow", "Locate the query point in a cell and show empty sphere of that cell", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionNormal_View->setText(QApplication::translate("MainWindow", "Normal Mode", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionNormal_View->setToolTip(QApplication::translate("MainWindow", "Normal Mode", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionNormal_View->setStatusTip(QApplication::translate("MainWindow", "Normal Mode", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionDemo_Help->setText(QApplication::translate("MainWindow", "Triangulation_3D Demo Help", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionDemo_Help->setToolTip(QApplication::translate("MainWindow", "Triangulation_3D Demo Help", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionDemo_Help->setStatusTip(QApplication::translate("MainWindow", "Triangulation_3D Demo Help", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionDemo_Help->setShortcut(QApplication::translate("MainWindow", "H", 0, QApplication::UnicodeUTF8));
    actionIncremental_Construct->setText(QApplication::translate("MainWindow", "Insertion Animation", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionIncremental_Construct->setToolTip(QApplication::translate("MainWindow", "Animation of incremental Delaunay triangulation", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionIncremental_Construct->setStatusTip(QApplication::translate("MainWindow", "Animation of incremental Delaunay triangulation", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionStop_Animation->setText(QApplication::translate("MainWindow", "Stop Animation", 0, QApplication::UnicodeUTF8));

#ifndef QT_NO_TOOLTIP
    actionStop_Animation->setToolTip(QApplication::translate("MainWindow", "Stop Animation", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP


#ifndef QT_NO_STATUSTIP
    actionStop_Animation->setStatusTip(QApplication::translate("MainWindow", "Stop Animation", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP

    actionAbout_T3_demo->setText(QApplication::translate("MainWindow", "About T3_demo", 0, QApplication::UnicodeUTF8));
    menuFile->setTitle(QApplication::translate("MainWindow", "&File", 0, QApplication::UnicodeUTF8));
    menuEdit->setTitle(QApplication::translate("MainWindow", "Edit", 0, QApplication::UnicodeUTF8));
    menuMode->setTitle(QApplication::translate("MainWindow", "Mode", 0, QApplication::UnicodeUTF8));
    menuShow->setTitle(QApplication::translate("MainWindow", "Show", 0, QApplication::UnicodeUTF8));
    menuHelp->setTitle(QApplication::translate("MainWindow", "Help", 0, QApplication::UnicodeUTF8));
    toolBar->setWindowTitle(QApplication::translate("MainWindow", "toolBar", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
