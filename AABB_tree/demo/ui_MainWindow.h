/********************************************************************************
** Form generated from reading ui file 'MainWindow.ui'
**
** Created: Sun 28. Jun 21:42:55 2009
**      by: Qt User Interface Compiler version 4.4.3
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QLocale>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QGridLayout>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QStatusBar>
#include <QtGui/QWidget>
#include "Viewer.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionQuit;
    QAction *actionLoadPolyhedron;
    QAction *actionInside_points;
    QAction *actionDo_intersect;
    QAction *actionAny_intersection;
    QAction *actionAll_intersections;
    QAction *actionBench_distances;
    QAction *actionNb_intersections;
    QAction *actionAll_intersected_primitives;
    QAction *actionUnsigned_distance_function_to_facets;
    QAction *actionUnsigned_distance_function_to_edges;
    QAction *actionSigned_distance_function_to_facets;
    QAction *actionView_polyhedron;
    QAction *actionView_points;
    QAction *actionClear_points;
    QAction *actionBoundary_segments;
    QAction *actionBoundary_points;
    QAction *actionClear_segments;
    QAction *actionView_segments;
    QAction *actionEdge_points;
    QAction *actionBench_intersections;
    QWidget *centralwidget;
    QGridLayout *gridLayout;
    Viewer *viewer;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuView;
    QMenu *menuAlgorithms;
    QMenu *menuBenchmarks;
    QMenu *menuEdit;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *MainWindow)
    {
    if (MainWindow->objectName().isEmpty())
        MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
    MainWindow->resize(978, 594);
    QIcon icon;
    icon.addPixmap(QPixmap(QString::fromUtf8(":/cgal/icons/resources/cgal_logo.xpm")), QIcon::Normal, QIcon::Off);
    MainWindow->setWindowIcon(icon);
    actionQuit = new QAction(MainWindow);
    actionQuit->setObjectName(QString::fromUtf8("actionQuit"));
    actionLoadPolyhedron = new QAction(MainWindow);
    actionLoadPolyhedron->setObjectName(QString::fromUtf8("actionLoadPolyhedron"));
    actionInside_points = new QAction(MainWindow);
    actionInside_points->setObjectName(QString::fromUtf8("actionInside_points"));
    actionDo_intersect = new QAction(MainWindow);
    actionDo_intersect->setObjectName(QString::fromUtf8("actionDo_intersect"));
    actionAny_intersection = new QAction(MainWindow);
    actionAny_intersection->setObjectName(QString::fromUtf8("actionAny_intersection"));
    actionAll_intersections = new QAction(MainWindow);
    actionAll_intersections->setObjectName(QString::fromUtf8("actionAll_intersections"));
    actionBench_distances = new QAction(MainWindow);
    actionBench_distances->setObjectName(QString::fromUtf8("actionBench_distances"));
    actionNb_intersections = new QAction(MainWindow);
    actionNb_intersections->setObjectName(QString::fromUtf8("actionNb_intersections"));
    actionAll_intersected_primitives = new QAction(MainWindow);
    actionAll_intersected_primitives->setObjectName(QString::fromUtf8("actionAll_intersected_primitives"));
    actionUnsigned_distance_function_to_facets = new QAction(MainWindow);
    actionUnsigned_distance_function_to_facets->setObjectName(QString::fromUtf8("actionUnsigned_distance_function_to_facets"));
    actionUnsigned_distance_function_to_edges = new QAction(MainWindow);
    actionUnsigned_distance_function_to_edges->setObjectName(QString::fromUtf8("actionUnsigned_distance_function_to_edges"));
    actionSigned_distance_function_to_facets = new QAction(MainWindow);
    actionSigned_distance_function_to_facets->setObjectName(QString::fromUtf8("actionSigned_distance_function_to_facets"));
    actionView_polyhedron = new QAction(MainWindow);
    actionView_polyhedron->setObjectName(QString::fromUtf8("actionView_polyhedron"));
    actionView_points = new QAction(MainWindow);
    actionView_points->setObjectName(QString::fromUtf8("actionView_points"));
    actionClear_points = new QAction(MainWindow);
    actionClear_points->setObjectName(QString::fromUtf8("actionClear_points"));
    actionBoundary_segments = new QAction(MainWindow);
    actionBoundary_segments->setObjectName(QString::fromUtf8("actionBoundary_segments"));
    actionBoundary_points = new QAction(MainWindow);
    actionBoundary_points->setObjectName(QString::fromUtf8("actionBoundary_points"));
    actionClear_segments = new QAction(MainWindow);
    actionClear_segments->setObjectName(QString::fromUtf8("actionClear_segments"));
    actionView_segments = new QAction(MainWindow);
    actionView_segments->setObjectName(QString::fromUtf8("actionView_segments"));
    actionEdge_points = new QAction(MainWindow);
    actionEdge_points->setObjectName(QString::fromUtf8("actionEdge_points"));
    actionBench_intersections = new QAction(MainWindow);
    actionBench_intersections->setObjectName(QString::fromUtf8("actionBench_intersections"));
    centralwidget = new QWidget(MainWindow);
    centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
    gridLayout = new QGridLayout(centralwidget);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
    viewer = new Viewer(centralwidget);
    viewer->setObjectName(QString::fromUtf8("viewer"));
    viewer->setLocale(QLocale(QLocale::English, QLocale::UnitedStates));

    gridLayout->addWidget(viewer, 0, 1, 1, 1);

    MainWindow->setCentralWidget(centralwidget);
    menubar = new QMenuBar(MainWindow);
    menubar->setObjectName(QString::fromUtf8("menubar"));
    menubar->setGeometry(QRect(0, 0, 978, 24));
    menuFile = new QMenu(menubar);
    menuFile->setObjectName(QString::fromUtf8("menuFile"));
    menuView = new QMenu(menubar);
    menuView->setObjectName(QString::fromUtf8("menuView"));
    menuAlgorithms = new QMenu(menubar);
    menuAlgorithms->setObjectName(QString::fromUtf8("menuAlgorithms"));
    menuBenchmarks = new QMenu(menubar);
    menuBenchmarks->setObjectName(QString::fromUtf8("menuBenchmarks"));
    menuEdit = new QMenu(menubar);
    menuEdit->setObjectName(QString::fromUtf8("menuEdit"));
    MainWindow->setMenuBar(menubar);
    statusbar = new QStatusBar(MainWindow);
    statusbar->setObjectName(QString::fromUtf8("statusbar"));
    MainWindow->setStatusBar(statusbar);

    menubar->addAction(menuFile->menuAction());
    menubar->addAction(menuEdit->menuAction());
    menubar->addAction(menuAlgorithms->menuAction());
    menubar->addAction(menuBenchmarks->menuAction());
    menubar->addAction(menuView->menuAction());
    menuFile->addAction(actionLoadPolyhedron);
    menuFile->addSeparator();
    menuFile->addAction(actionQuit);
    menuView->addAction(actionView_polyhedron);
    menuView->addAction(actionView_points);
    menuView->addAction(actionView_segments);
    menuAlgorithms->addAction(actionEdge_points);
    menuAlgorithms->addAction(actionInside_points);
    menuAlgorithms->addAction(actionBoundary_points);
    menuAlgorithms->addAction(actionBoundary_segments);
    menuAlgorithms->addSeparator();
    menuAlgorithms->addAction(actionSigned_distance_function_to_facets);
    menuAlgorithms->addAction(actionUnsigned_distance_function_to_facets);
    menuAlgorithms->addAction(actionUnsigned_distance_function_to_edges);
    menuBenchmarks->addAction(actionBench_distances);
    menuBenchmarks->addAction(actionBench_intersections);
    menuEdit->addAction(actionClear_points);
    menuEdit->addAction(actionClear_segments);

    retranslateUi(MainWindow);

    QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
    MainWindow->setWindowTitle(QApplication::translate("MainWindow", "CGAL AABB tree demo", 0, QApplication::UnicodeUTF8));
    actionQuit->setText(QApplication::translate("MainWindow", "&Quit", 0, QApplication::UnicodeUTF8));
    actionQuit->setShortcut(QApplication::translate("MainWindow", "Ctrl+Q", 0, QApplication::UnicodeUTF8));
    actionLoadPolyhedron->setText(QApplication::translate("MainWindow", "Load polyhedron...", 0, QApplication::UnicodeUTF8));
    actionInside_points->setText(QApplication::translate("MainWindow", "Inside points...", 0, QApplication::UnicodeUTF8));
    actionDo_intersect->setText(QApplication::translate("MainWindow", "do_intersect", 0, QApplication::UnicodeUTF8));
    actionAny_intersection->setText(QApplication::translate("MainWindow", "any_intersection", 0, QApplication::UnicodeUTF8));
    actionAll_intersections->setText(QApplication::translate("MainWindow", "all_intersections", 0, QApplication::UnicodeUTF8));
    actionBench_distances->setText(QApplication::translate("MainWindow", "Distances", 0, QApplication::UnicodeUTF8));
    actionNb_intersections->setText(QApplication::translate("MainWindow", "nb_intersections", 0, QApplication::UnicodeUTF8));
    actionAll_intersected_primitives->setText(QApplication::translate("MainWindow", "all_intersected_primitives", 0, QApplication::UnicodeUTF8));
    actionUnsigned_distance_function_to_facets->setText(QApplication::translate("MainWindow", "Unsigned distance function to facets", 0, QApplication::UnicodeUTF8));
    actionUnsigned_distance_function_to_edges->setText(QApplication::translate("MainWindow", "Unsigned distance function to edges", 0, QApplication::UnicodeUTF8));
    actionSigned_distance_function_to_facets->setText(QApplication::translate("MainWindow", "Signed distance function to facets", 0, QApplication::UnicodeUTF8));
    actionView_polyhedron->setText(QApplication::translate("MainWindow", "Polyhedron", 0, QApplication::UnicodeUTF8));
    actionView_points->setText(QApplication::translate("MainWindow", "Points", 0, QApplication::UnicodeUTF8));
    actionClear_points->setText(QApplication::translate("MainWindow", "Clear points", 0, QApplication::UnicodeUTF8));
    actionBoundary_segments->setText(QApplication::translate("MainWindow", "Boundary segments...", 0, QApplication::UnicodeUTF8));
    actionBoundary_points->setText(QApplication::translate("MainWindow", "Boundary points...", 0, QApplication::UnicodeUTF8));
    actionClear_segments->setText(QApplication::translate("MainWindow", "Clear segments", 0, QApplication::UnicodeUTF8));
    actionView_segments->setText(QApplication::translate("MainWindow", "Segments", 0, QApplication::UnicodeUTF8));
    actionEdge_points->setText(QApplication::translate("MainWindow", "Edge points...", 0, QApplication::UnicodeUTF8));
    actionBench_intersections->setText(QApplication::translate("MainWindow", "Intersections", 0, QApplication::UnicodeUTF8));
    menuFile->setTitle(QApplication::translate("MainWindow", "&File", 0, QApplication::UnicodeUTF8));
    menuView->setTitle(QApplication::translate("MainWindow", "&View", 0, QApplication::UnicodeUTF8));
    menuAlgorithms->setTitle(QApplication::translate("MainWindow", "Algorithms", 0, QApplication::UnicodeUTF8));
    menuBenchmarks->setTitle(QApplication::translate("MainWindow", "Benchmark", 0, QApplication::UnicodeUTF8));
    menuEdit->setTitle(QApplication::translate("MainWindow", "Edit", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
