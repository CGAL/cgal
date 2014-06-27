/********************************************************************************
** Form generated from reading UI file 'pwsrec.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PWSREC_H
#define UI_PWSREC_H

#include <QtCore/QLocale>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QSlider>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QStatusBar>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "glviewer.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *actionQuit;
    QAction *actionInsertPoint;
    QAction *actionClear;
    QAction *actionLoadPoints;
    QAction *actionSave;
    QAction *actionCircle;
    QAction *actionHalf_circle;
    QAction *actionBox;
    QAction *actionLine;
    QAction *actionReconstruction_init;
    QAction *actionReconstruction_one_step;
    QAction *actionView_foot_points;
    QAction *actionView_points;
    QAction *actionView_edges;
    QAction *actionRecenter;
    QAction *actionView_vertices;
    QAction *actionBoxes;
    QAction *actionStair;
    QAction *actionSkyline;
    QAction *actionView_edge_priority;
    QAction *actionReconstruction_10_steps;
    QAction *actionReconstruction_100_steps;
    QAction *actionReconstruction_1000_steps;
    QAction *actionAdd_outliers;
    QAction *actionSnapshot;
    QAction *actionIncreasingly_sharp_angles;
    QAction *actionBox_with_boundaries;
    QAction *actionBox_with_missing_corners;
    QAction *actionStar;
    QAction *actionSpiral;
    QAction *actionSet_parameters;
    QAction *actionView_edge_cost;
    QAction *actionReconstruction_until;
    QAction *actionParallel_lines;
    QAction *actionNoise;
    QAction *actionActivate_simulation;
    QAction *actionView_simulation;
    QAction *actionRelocate_vertices;
    QAction *actionView_relocation;
    QAction *actionView_ghost;
    QAction *actionInvert_mass;
    QAction *actionView_relevance;
    QAction *actionView_tolerance;
    QAction *actionView_incolors;
    QAction *actionClamp_mass;
    QAction *actionView_bins;
    QAction *actionPrint_Stats;
    QAction *actionSubdivide;
    QAction *actionWidely_variable_sampling;
    QAction *actionDecimate;
    QAction *actionKeep_one_point_out_of_n;
    QAction *actionSet_MChoice;
    QWidget *centralwidget;
    QGridLayout *gridLayout;
    QVBoxLayout *vboxLayout;
    GlViewer *viewer;
    QVBoxLayout *vboxLayout1;
    QSlider *min_mass_slider;
    QHBoxLayout *hboxLayout;
    QSpacerItem *spacerItem;
    QHBoxLayout *hboxLayout1;
    QLabel *label;
    QSpinBox *discard_spinbox;
    QStatusBar *statusbar;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuPoint_set;
    QMenu *menuPredefined;
    QMenu *menuAlgorithms;
    QMenu *menuReconstruction;
    QMenu *menuView;
    QToolBar *toolBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(680, 680);
        actionQuit = new QAction(MainWindow);
        actionQuit->setObjectName(QString::fromUtf8("actionQuit"));
        actionInsertPoint = new QAction(MainWindow);
        actionInsertPoint->setObjectName(QString::fromUtf8("actionInsertPoint"));
        actionInsertPoint->setCheckable(true);
        actionInsertPoint->setChecked(false);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/icons/inputPoint.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionInsertPoint->setIcon(icon);
        actionClear = new QAction(MainWindow);
        actionClear->setObjectName(QString::fromUtf8("actionClear"));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/icons/fileNew.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionClear->setIcon(icon1);
        actionLoadPoints = new QAction(MainWindow);
        actionLoadPoints->setObjectName(QString::fromUtf8("actionLoadPoints"));
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/icons/fileOpen.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLoadPoints->setIcon(icon2);
        actionSave = new QAction(MainWindow);
        actionSave->setObjectName(QString::fromUtf8("actionSave"));
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/icons/fileSave.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSave->setIcon(icon3);
        actionCircle = new QAction(MainWindow);
        actionCircle->setObjectName(QString::fromUtf8("actionCircle"));
        actionHalf_circle = new QAction(MainWindow);
        actionHalf_circle->setObjectName(QString::fromUtf8("actionHalf_circle"));
        actionBox = new QAction(MainWindow);
        actionBox->setObjectName(QString::fromUtf8("actionBox"));
        actionLine = new QAction(MainWindow);
        actionLine->setObjectName(QString::fromUtf8("actionLine"));
        actionReconstruction_init = new QAction(MainWindow);
        actionReconstruction_init->setObjectName(QString::fromUtf8("actionReconstruction_init"));
        actionReconstruction_one_step = new QAction(MainWindow);
        actionReconstruction_one_step->setObjectName(QString::fromUtf8("actionReconstruction_one_step"));
        actionView_foot_points = new QAction(MainWindow);
        actionView_foot_points->setObjectName(QString::fromUtf8("actionView_foot_points"));
        actionView_foot_points->setCheckable(true);
        actionView_foot_points->setChecked(false);
        actionView_points = new QAction(MainWindow);
        actionView_points->setObjectName(QString::fromUtf8("actionView_points"));
        actionView_points->setCheckable(true);
        actionView_points->setChecked(true);
        actionView_edges = new QAction(MainWindow);
        actionView_edges->setObjectName(QString::fromUtf8("actionView_edges"));
        actionView_edges->setCheckable(true);
        actionView_edges->setChecked(false);
        actionRecenter = new QAction(MainWindow);
        actionRecenter->setObjectName(QString::fromUtf8("actionRecenter"));
        QIcon icon4;
        icon4.addFile(QString::fromUtf8(":/icons/fit-page-32.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionRecenter->setIcon(icon4);
        actionView_vertices = new QAction(MainWindow);
        actionView_vertices->setObjectName(QString::fromUtf8("actionView_vertices"));
        actionView_vertices->setCheckable(true);
        actionView_vertices->setChecked(true);
        actionBoxes = new QAction(MainWindow);
        actionBoxes->setObjectName(QString::fromUtf8("actionBoxes"));
        actionStair = new QAction(MainWindow);
        actionStair->setObjectName(QString::fromUtf8("actionStair"));
        actionSkyline = new QAction(MainWindow);
        actionSkyline->setObjectName(QString::fromUtf8("actionSkyline"));
        actionView_edge_priority = new QAction(MainWindow);
        actionView_edge_priority->setObjectName(QString::fromUtf8("actionView_edge_priority"));
        actionView_edge_priority->setCheckable(true);
        actionView_edge_priority->setChecked(false);
        actionReconstruction_10_steps = new QAction(MainWindow);
        actionReconstruction_10_steps->setObjectName(QString::fromUtf8("actionReconstruction_10_steps"));
        actionReconstruction_100_steps = new QAction(MainWindow);
        actionReconstruction_100_steps->setObjectName(QString::fromUtf8("actionReconstruction_100_steps"));
        actionReconstruction_1000_steps = new QAction(MainWindow);
        actionReconstruction_1000_steps->setObjectName(QString::fromUtf8("actionReconstruction_1000_steps"));
        actionAdd_outliers = new QAction(MainWindow);
        actionAdd_outliers->setObjectName(QString::fromUtf8("actionAdd_outliers"));
        actionSnapshot = new QAction(MainWindow);
        actionSnapshot->setObjectName(QString::fromUtf8("actionSnapshot"));
        QIcon icon5;
        icon5.addFile(QString::fromUtf8(":/icons/snapshot.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSnapshot->setIcon(icon5);
        actionIncreasingly_sharp_angles = new QAction(MainWindow);
        actionIncreasingly_sharp_angles->setObjectName(QString::fromUtf8("actionIncreasingly_sharp_angles"));
        actionBox_with_boundaries = new QAction(MainWindow);
        actionBox_with_boundaries->setObjectName(QString::fromUtf8("actionBox_with_boundaries"));
        actionBox_with_missing_corners = new QAction(MainWindow);
        actionBox_with_missing_corners->setObjectName(QString::fromUtf8("actionBox_with_missing_corners"));
        actionStar = new QAction(MainWindow);
        actionStar->setObjectName(QString::fromUtf8("actionStar"));
        actionSpiral = new QAction(MainWindow);
        actionSpiral->setObjectName(QString::fromUtf8("actionSpiral"));
        actionSet_parameters = new QAction(MainWindow);
        actionSet_parameters->setObjectName(QString::fromUtf8("actionSet_parameters"));
        actionView_edge_cost = new QAction(MainWindow);
        actionView_edge_cost->setObjectName(QString::fromUtf8("actionView_edge_cost"));
        actionView_edge_cost->setCheckable(true);
        actionView_edge_cost->setChecked(false);
        actionReconstruction_until = new QAction(MainWindow);
        actionReconstruction_until->setObjectName(QString::fromUtf8("actionReconstruction_until"));
        actionParallel_lines = new QAction(MainWindow);
        actionParallel_lines->setObjectName(QString::fromUtf8("actionParallel_lines"));
        actionNoise = new QAction(MainWindow);
        actionNoise->setObjectName(QString::fromUtf8("actionNoise"));
        actionActivate_simulation = new QAction(MainWindow);
        actionActivate_simulation->setObjectName(QString::fromUtf8("actionActivate_simulation"));
        actionActivate_simulation->setCheckable(true);
        QIcon icon6;
        icon6.addFile(QString::fromUtf8(":/icons/vertex.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionActivate_simulation->setIcon(icon6);
        actionView_simulation = new QAction(MainWindow);
        actionView_simulation->setObjectName(QString::fromUtf8("actionView_simulation"));
        actionRelocate_vertices = new QAction(MainWindow);
        actionRelocate_vertices->setObjectName(QString::fromUtf8("actionRelocate_vertices"));
        actionView_relocation = new QAction(MainWindow);
        actionView_relocation->setObjectName(QString::fromUtf8("actionView_relocation"));
        actionView_relocation->setCheckable(true);
        actionView_ghost = new QAction(MainWindow);
        actionView_ghost->setObjectName(QString::fromUtf8("actionView_ghost"));
        actionView_ghost->setCheckable(true);
        actionView_ghost->setChecked(false);
        actionInvert_mass = new QAction(MainWindow);
        actionInvert_mass->setObjectName(QString::fromUtf8("actionInvert_mass"));
        actionView_relevance = new QAction(MainWindow);
        actionView_relevance->setObjectName(QString::fromUtf8("actionView_relevance"));
        actionView_relevance->setCheckable(true);
        actionView_relevance->setChecked(true);
        actionView_tolerance = new QAction(MainWindow);
        actionView_tolerance->setObjectName(QString::fromUtf8("actionView_tolerance"));
        actionView_tolerance->setCheckable(true);
        actionView_incolors = new QAction(MainWindow);
        actionView_incolors->setObjectName(QString::fromUtf8("actionView_incolors"));
        actionView_incolors->setCheckable(true);
        actionView_incolors->setChecked(false);
        actionClamp_mass = new QAction(MainWindow);
        actionClamp_mass->setObjectName(QString::fromUtf8("actionClamp_mass"));
        actionView_bins = new QAction(MainWindow);
        actionView_bins->setObjectName(QString::fromUtf8("actionView_bins"));
        actionView_bins->setCheckable(true);
        actionPrint_Stats = new QAction(MainWindow);
        actionPrint_Stats->setObjectName(QString::fromUtf8("actionPrint_Stats"));
        actionSubdivide = new QAction(MainWindow);
        actionSubdivide->setObjectName(QString::fromUtf8("actionSubdivide"));
        actionWidely_variable_sampling = new QAction(MainWindow);
        actionWidely_variable_sampling->setObjectName(QString::fromUtf8("actionWidely_variable_sampling"));
        actionDecimate = new QAction(MainWindow);
        actionDecimate->setObjectName(QString::fromUtf8("actionDecimate"));
        actionKeep_one_point_out_of_n = new QAction(MainWindow);
        actionKeep_one_point_out_of_n->setObjectName(QString::fromUtf8("actionKeep_one_point_out_of_n"));
        actionSet_MChoice = new QAction(MainWindow);
        actionSet_MChoice->setObjectName(QString::fromUtf8("actionSet_MChoice"));
        actionSet_MChoice->setCheckable(true);
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        gridLayout = new QGridLayout(centralwidget);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        vboxLayout = new QVBoxLayout();
        vboxLayout->setObjectName(QString::fromUtf8("vboxLayout"));
        viewer = new GlViewer(centralwidget);
        viewer->setObjectName(QString::fromUtf8("viewer"));
        viewer->setFocusPolicy(Qt::StrongFocus);
        viewer->setLocale(QLocale(QLocale::English, QLocale::UnitedStates));

        vboxLayout->addWidget(viewer);

        vboxLayout1 = new QVBoxLayout();
        vboxLayout1->setObjectName(QString::fromUtf8("vboxLayout1"));
        min_mass_slider = new QSlider(centralwidget);
        min_mass_slider->setObjectName(QString::fromUtf8("min_mass_slider"));
        min_mass_slider->setOrientation(Qt::Horizontal);

        vboxLayout1->addWidget(min_mass_slider);

        hboxLayout = new QHBoxLayout();
        hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
        spacerItem = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        hboxLayout->addItem(spacerItem);

        hboxLayout1 = new QHBoxLayout();
        hboxLayout1->setObjectName(QString::fromUtf8("hboxLayout1"));
        label = new QLabel(centralwidget);
        label->setObjectName(QString::fromUtf8("label"));

        hboxLayout1->addWidget(label);

        discard_spinbox = new QSpinBox(centralwidget);
        discard_spinbox->setObjectName(QString::fromUtf8("discard_spinbox"));
        discard_spinbox->setMaximum(10000);

        hboxLayout1->addWidget(discard_spinbox);


        hboxLayout->addLayout(hboxLayout1);


        vboxLayout1->addLayout(hboxLayout);


        vboxLayout->addLayout(vboxLayout1);


        gridLayout->addLayout(vboxLayout, 0, 0, 1, 1);

        MainWindow->setCentralWidget(centralwidget);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        MainWindow->setStatusBar(statusbar);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 680, 22));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuPoint_set = new QMenu(menubar);
        menuPoint_set->setObjectName(QString::fromUtf8("menuPoint_set"));
        menuPredefined = new QMenu(menuPoint_set);
        menuPredefined->setObjectName(QString::fromUtf8("menuPredefined"));
        menuAlgorithms = new QMenu(menubar);
        menuAlgorithms->setObjectName(QString::fromUtf8("menuAlgorithms"));
        menuReconstruction = new QMenu(menuAlgorithms);
        menuReconstruction->setObjectName(QString::fromUtf8("menuReconstruction"));
        menuView = new QMenu(menubar);
        menuView->setObjectName(QString::fromUtf8("menuView"));
        MainWindow->setMenuBar(menubar);
        toolBar = new QToolBar(MainWindow);
        toolBar->setObjectName(QString::fromUtf8("toolBar"));
        toolBar->setLocale(QLocale(QLocale::English, QLocale::UnitedStates));
        MainWindow->addToolBar(Qt::TopToolBarArea, toolBar);
        QWidget::setTabOrder(min_mass_slider, discard_spinbox);

        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuPoint_set->menuAction());
        menubar->addAction(menuAlgorithms->menuAction());
        menubar->addAction(menuView->menuAction());
        menuFile->addAction(actionLoadPoints);
        menuFile->addAction(actionSave);
        menuFile->addAction(actionSnapshot);
        menuFile->addAction(actionQuit);
        menuPoint_set->addAction(menuPredefined->menuAction());
        menuPoint_set->addAction(actionNoise);
        menuPoint_set->addAction(actionAdd_outliers);
        menuPoint_set->addSeparator();
        menuPoint_set->addAction(actionDecimate);
        menuPoint_set->addAction(actionKeep_one_point_out_of_n);
        menuPoint_set->addAction(actionSubdivide);
        menuPoint_set->addSeparator();
        menuPoint_set->addAction(actionInvert_mass);
        menuPoint_set->addAction(actionClamp_mass);
        menuPoint_set->addSeparator();
        menuPoint_set->addAction(actionInsertPoint);
        menuPoint_set->addSeparator();
        menuPoint_set->addAction(actionClear);
        menuPredefined->addAction(actionLine);
        menuPredefined->addAction(actionParallel_lines);
        menuPredefined->addAction(actionCircle);
        menuPredefined->addAction(actionHalf_circle);
        menuPredefined->addAction(actionWidely_variable_sampling);
        menuPredefined->addAction(actionSpiral);
        menuPredefined->addAction(actionBox);
        menuPredefined->addAction(actionBoxes);
        menuPredefined->addAction(actionBox_with_boundaries);
        menuPredefined->addAction(actionBox_with_missing_corners);
        menuPredefined->addAction(actionStar);
        menuPredefined->addAction(actionStair);
        menuPredefined->addAction(actionSkyline);
        menuPredefined->addAction(actionIncreasingly_sharp_angles);
        menuAlgorithms->addAction(actionSet_parameters);
        menuAlgorithms->addAction(actionSet_MChoice);
        menuAlgorithms->addSeparator();
        menuAlgorithms->addAction(actionReconstruction_init);
        menuAlgorithms->addAction(menuReconstruction->menuAction());
        menuAlgorithms->addAction(actionRelocate_vertices);
        menuAlgorithms->addSeparator();
        menuAlgorithms->addAction(actionPrint_Stats);
        menuReconstruction->addAction(actionReconstruction_one_step);
        menuReconstruction->addAction(actionReconstruction_10_steps);
        menuReconstruction->addAction(actionReconstruction_100_steps);
        menuReconstruction->addAction(actionReconstruction_1000_steps);
        menuReconstruction->addAction(actionReconstruction_until);
        menuView->addAction(actionView_points);
        menuView->addSeparator();
        menuView->addAction(actionView_vertices);
        menuView->addAction(actionView_edges);
        menuView->addSeparator();
        menuView->addAction(actionView_ghost);
        menuView->addAction(actionView_relevance);
        menuView->addAction(actionView_incolors);
        menuView->addSeparator();
        menuView->addAction(actionView_edge_cost);
        menuView->addAction(actionView_edge_priority);
        menuView->addSeparator();
        menuView->addAction(actionView_bins);
        menuView->addAction(actionView_foot_points);
        menuView->addAction(actionView_relocation);
        menuView->addAction(actionView_tolerance);
        menuView->addSeparator();
        menuView->addAction(actionActivate_simulation);
        menuView->addAction(actionView_simulation);
        menuView->addSeparator();
        menuView->addAction(actionRecenter);
        toolBar->addAction(actionClear);
        toolBar->addAction(actionLoadPoints);
        toolBar->addAction(actionSave);
        toolBar->addAction(actionSnapshot);
        toolBar->addSeparator();
        toolBar->addAction(actionInsertPoint);
        toolBar->addAction(actionActivate_simulation);
        toolBar->addSeparator();
        toolBar->addAction(actionRecenter);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "Reconstruction_simplification_2 Demo", 0, QApplication::UnicodeUTF8));
        actionQuit->setText(QApplication::translate("MainWindow", "Quit", 0, QApplication::UnicodeUTF8));
        actionQuit->setShortcut(QApplication::translate("MainWindow", "Ctrl+Q", 0, QApplication::UnicodeUTF8));
        actionInsertPoint->setText(QApplication::translate("MainWindow", "Insert mode", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionInsertPoint->setToolTip(QApplication::translate("MainWindow", "Insert Point", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_STATUSTIP
        actionInsertPoint->setStatusTip(QApplication::translate("MainWindow", "Insert Point", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP
        actionClear->setText(QApplication::translate("MainWindow", "Clear", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_STATUSTIP
        actionClear->setStatusTip(QApplication::translate("MainWindow", "Clear", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP
        actionClear->setShortcut(QApplication::translate("MainWindow", "Space", 0, QApplication::UnicodeUTF8));
        actionLoadPoints->setText(QApplication::translate("MainWindow", "Load Points", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_STATUSTIP
        actionLoadPoints->setStatusTip(QApplication::translate("MainWindow", "Load Points", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP
        actionLoadPoints->setShortcut(QApplication::translate("MainWindow", "Ctrl+O", 0, QApplication::UnicodeUTF8));
        actionSave->setText(QApplication::translate("MainWindow", "Save", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_STATUSTIP
        actionSave->setStatusTip(QApplication::translate("MainWindow", "Save Points", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP
        actionSave->setShortcut(QApplication::translate("MainWindow", "Ctrl+S", 0, QApplication::UnicodeUTF8));
        actionCircle->setText(QApplication::translate("MainWindow", "Circle...", 0, QApplication::UnicodeUTF8));
        actionHalf_circle->setText(QApplication::translate("MainWindow", "Half circle...", 0, QApplication::UnicodeUTF8));
        actionBox->setText(QApplication::translate("MainWindow", "Box...", 0, QApplication::UnicodeUTF8));
        actionBox->setShortcut(QApplication::translate("MainWindow", "B", 0, QApplication::UnicodeUTF8));
        actionLine->setText(QApplication::translate("MainWindow", "Line...", 0, QApplication::UnicodeUTF8));
        actionReconstruction_init->setText(QApplication::translate("MainWindow", "Init", 0, QApplication::UnicodeUTF8));
        actionReconstruction_init->setShortcut(QApplication::translate("MainWindow", "I", 0, QApplication::UnicodeUTF8));
        actionReconstruction_one_step->setText(QApplication::translate("MainWindow", "1 step", 0, QApplication::UnicodeUTF8));
        actionReconstruction_one_step->setShortcut(QApplication::translate("MainWindow", "R", 0, QApplication::UnicodeUTF8));
        actionView_foot_points->setText(QApplication::translate("MainWindow", "Foot points", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_STATUSTIP
        actionView_foot_points->setStatusTip(QApplication::translate("MainWindow", "View foot points", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP
        actionView_foot_points->setShortcut(QApplication::translate("MainWindow", "T", 0, QApplication::UnicodeUTF8));
        actionView_points->setText(QApplication::translate("MainWindow", "Points", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_STATUSTIP
        actionView_points->setStatusTip(QApplication::translate("MainWindow", "View points", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP
        actionView_points->setShortcut(QApplication::translate("MainWindow", "P", 0, QApplication::UnicodeUTF8));
        actionView_edges->setText(QApplication::translate("MainWindow", "Edges", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_STATUSTIP
        actionView_edges->setStatusTip(QApplication::translate("MainWindow", "View edges", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_STATUSTIP
        actionView_edges->setShortcut(QApplication::translate("MainWindow", "E", 0, QApplication::UnicodeUTF8));
        actionRecenter->setText(QApplication::translate("MainWindow", "Recenter", 0, QApplication::UnicodeUTF8));
        actionView_vertices->setText(QApplication::translate("MainWindow", "Vertices", 0, QApplication::UnicodeUTF8));
        actionView_vertices->setShortcut(QApplication::translate("MainWindow", "V", 0, QApplication::UnicodeUTF8));
        actionBoxes->setText(QApplication::translate("MainWindow", "Two boxes...", 0, QApplication::UnicodeUTF8));
        actionStair->setText(QApplication::translate("MainWindow", "Stair...", 0, QApplication::UnicodeUTF8));
        actionSkyline->setText(QApplication::translate("MainWindow", "Skyline...", 0, QApplication::UnicodeUTF8));
        actionView_edge_priority->setText(QApplication::translate("MainWindow", "Edge priority", 0, QApplication::UnicodeUTF8));
        actionView_edge_priority->setShortcut(QApplication::translate("MainWindow", "Z", 0, QApplication::UnicodeUTF8));
        actionReconstruction_10_steps->setText(QApplication::translate("MainWindow", "10 steps", 0, QApplication::UnicodeUTF8));
        actionReconstruction_10_steps->setShortcut(QApplication::translate("MainWindow", "1", 0, QApplication::UnicodeUTF8));
        actionReconstruction_100_steps->setText(QApplication::translate("MainWindow", "100 steps", 0, QApplication::UnicodeUTF8));
        actionReconstruction_100_steps->setShortcut(QApplication::translate("MainWindow", "2", 0, QApplication::UnicodeUTF8));
        actionReconstruction_1000_steps->setText(QApplication::translate("MainWindow", "1000 steps", 0, QApplication::UnicodeUTF8));
        actionReconstruction_1000_steps->setShortcut(QApplication::translate("MainWindow", "3", 0, QApplication::UnicodeUTF8));
        actionAdd_outliers->setText(QApplication::translate("MainWindow", "Add outliers", 0, QApplication::UnicodeUTF8));
        actionAdd_outliers->setShortcut(QApplication::translate("MainWindow", "O", 0, QApplication::UnicodeUTF8));
        actionSnapshot->setText(QApplication::translate("MainWindow", "Snapshot", 0, QApplication::UnicodeUTF8));
        actionSnapshot->setShortcut(QApplication::translate("MainWindow", "Ctrl+C", 0, QApplication::UnicodeUTF8));
        actionIncreasingly_sharp_angles->setText(QApplication::translate("MainWindow", "Increasingly sharp angles...", 0, QApplication::UnicodeUTF8));
        actionBox_with_boundaries->setText(QApplication::translate("MainWindow", "Box with boundaries...", 0, QApplication::UnicodeUTF8));
        actionBox_with_missing_corners->setText(QApplication::translate("MainWindow", "Box with missing corners...", 0, QApplication::UnicodeUTF8));
        actionStar->setText(QApplication::translate("MainWindow", "Star...", 0, QApplication::UnicodeUTF8));
        actionSpiral->setText(QApplication::translate("MainWindow", "Spiral...", 0, QApplication::UnicodeUTF8));
        actionSet_parameters->setText(QApplication::translate("MainWindow", "Parameters", 0, QApplication::UnicodeUTF8));
        actionSet_parameters->setShortcut(QApplication::translate("MainWindow", "Ctrl+P", 0, QApplication::UnicodeUTF8));
        actionView_edge_cost->setText(QApplication::translate("MainWindow", "Edge cost", 0, QApplication::UnicodeUTF8));
        actionView_edge_cost->setShortcut(QApplication::translate("MainWindow", "C", 0, QApplication::UnicodeUTF8));
        actionReconstruction_until->setText(QApplication::translate("MainWindow", "until", 0, QApplication::UnicodeUTF8));
        actionReconstruction_until->setShortcut(QApplication::translate("MainWindow", "U", 0, QApplication::UnicodeUTF8));
        actionParallel_lines->setText(QApplication::translate("MainWindow", "Parallel lines...", 0, QApplication::UnicodeUTF8));
        actionNoise->setText(QApplication::translate("MainWindow", "Noise", 0, QApplication::UnicodeUTF8));
        actionNoise->setShortcut(QApplication::translate("MainWindow", "N", 0, QApplication::UnicodeUTF8));
        actionActivate_simulation->setText(QApplication::translate("MainWindow", "Simulation", 0, QApplication::UnicodeUTF8));
        actionView_simulation->setText(QApplication::translate("MainWindow", "Simulation stage", 0, QApplication::UnicodeUTF8));
        actionView_simulation->setShortcut(QApplication::translate("MainWindow", "A", 0, QApplication::UnicodeUTF8));
        actionRelocate_vertices->setText(QApplication::translate("MainWindow", "Relocate", 0, QApplication::UnicodeUTF8));
        actionRelocate_vertices->setShortcut(QApplication::translate("MainWindow", "L", 0, QApplication::UnicodeUTF8));
        actionView_relocation->setText(QApplication::translate("MainWindow", "Relocation", 0, QApplication::UnicodeUTF8));
        actionView_relocation->setShortcut(QApplication::translate("MainWindow", "Shift+L", 0, QApplication::UnicodeUTF8));
        actionView_ghost->setText(QApplication::translate("MainWindow", "Ghost edges", 0, QApplication::UnicodeUTF8));
        actionView_ghost->setShortcut(QApplication::translate("MainWindow", "G", 0, QApplication::UnicodeUTF8));
        actionInvert_mass->setText(QApplication::translate("MainWindow", "Invert mass", 0, QApplication::UnicodeUTF8));
        actionInvert_mass->setShortcut(QApplication::translate("MainWindow", "Shift+I", 0, QApplication::UnicodeUTF8));
        actionView_relevance->setText(QApplication::translate("MainWindow", "Relevance", 0, QApplication::UnicodeUTF8));
        actionView_relevance->setShortcut(QApplication::translate("MainWindow", "Shift+R", 0, QApplication::UnicodeUTF8));
        actionView_tolerance->setText(QApplication::translate("MainWindow", "Tolerance", 0, QApplication::UnicodeUTF8));
        actionView_tolerance->setShortcut(QApplication::translate("MainWindow", "Shift+T", 0, QApplication::UnicodeUTF8));
        actionView_incolors->setText(QApplication::translate("MainWindow", "In colors", 0, QApplication::UnicodeUTF8));
        actionClamp_mass->setText(QApplication::translate("MainWindow", "Clamp mass", 0, QApplication::UnicodeUTF8));
        actionView_bins->setText(QApplication::translate("MainWindow", "Bins", 0, QApplication::UnicodeUTF8));
        actionView_bins->setShortcut(QApplication::translate("MainWindow", "Shift+B", 0, QApplication::UnicodeUTF8));
        actionPrint_Stats->setText(QApplication::translate("MainWindow", "Print Stats", 0, QApplication::UnicodeUTF8));
        actionPrint_Stats->setShortcut(QApplication::translate("MainWindow", "Shift+S", 0, QApplication::UnicodeUTF8));
        actionSubdivide->setText(QApplication::translate("MainWindow", "Subdivide", 0, QApplication::UnicodeUTF8));
        actionWidely_variable_sampling->setText(QApplication::translate("MainWindow", "Widely variable sampling", 0, QApplication::UnicodeUTF8));
        actionDecimate->setText(QApplication::translate("MainWindow", "Decimate", 0, QApplication::UnicodeUTF8));
        actionKeep_one_point_out_of_n->setText(QApplication::translate("MainWindow", "One point out of n", 0, QApplication::UnicodeUTF8));
        actionSet_MChoice->setText(QApplication::translate("MainWindow", "Multiple Choice", 0, QApplication::UnicodeUTF8));
        actionSet_MChoice->setShortcut(QApplication::translate("MainWindow", "M", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "Discard", 0, QApplication::UnicodeUTF8));
        menuFile->setTitle(QApplication::translate("MainWindow", "&File", 0, QApplication::UnicodeUTF8));
        menuPoint_set->setTitle(QApplication::translate("MainWindow", "Data", 0, QApplication::UnicodeUTF8));
        menuPredefined->setTitle(QApplication::translate("MainWindow", "Predefined", 0, QApplication::UnicodeUTF8));
        menuAlgorithms->setTitle(QApplication::translate("MainWindow", "Algorithms", 0, QApplication::UnicodeUTF8));
        menuReconstruction->setTitle(QApplication::translate("MainWindow", "Decimate", 0, QApplication::UnicodeUTF8));
        menuView->setTitle(QApplication::translate("MainWindow", "View", 0, QApplication::UnicodeUTF8));
        toolBar->setWindowTitle(QApplication::translate("MainWindow", "toolBar", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PWSREC_H
