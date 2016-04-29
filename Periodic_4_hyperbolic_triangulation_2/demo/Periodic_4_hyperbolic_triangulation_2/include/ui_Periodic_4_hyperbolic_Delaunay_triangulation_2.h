/********************************************************************************
** Form generated from reading UI file 'Hyperbolic_translations_2.ui'
**
** Created by: Qt User Interface Compiler version 4.8.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_H
#define UI_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_H

#if QT_VERSION >= 0x050000

    #include <QVariant>
    #include <QAction>
    #include <QApplication>
    #include <QButtonGroup>
    #include <QGraphicsView>
    #include <QGridLayout>
    #include <QHeaderView>
    #include <QMainWindow>
    #include <QMenu>
    #include <QMenuBar>
    #include <QStatusBar>
    #include <QToolBar>
    #include <QWidget>

#else

    #include <QtCore/QVariant>
    #include <QtGui/QAction>
    #include <QtGui/QApplication>
    #include <QtGui/QButtonGroup>
    #include <QtGui/QGraphicsView>
    #include <QtGui/QGridLayout>
    #include <QtGui/QHeaderView>
    #include <QtGui/QMainWindow>
    #include <QtGui/QMenu>
    #include <QtGui/QMenuBar>
    #include <QtGui/QStatusBar>
    #include <QtGui/QToolBar>
    #include <QtGui/QWidget>

#endif

QT_BEGIN_NAMESPACE

class Ui_Periodic_4_hyperbolic_Delaunay_triangulation_2
{
public:
    QAction *actionAbout;
    QAction *actionAboutCGAL;
    QAction *actionQuit;
    QAction *actionInsertRandomPoints;
    QAction *actionMovingPoint;
    QAction *actionInsertPoint;
    QAction *actionClear;
    QAction *actionInsertOrigin;
    QAction *actionInsertDummyPoints;
    QAction *actionModifyDepth;
    QAction *actionShowVoronoi;
    QAction *actionShowDelaunay;
    QAction *actionLoadPoints;
    QAction *actionSavePoints;
    QAction *actionCircumcenter;
    QAction *actionRecenter;
    QAction *actionShowConflictZone;
    QAction *actionDo_translation_a;
    QAction *actionDo_translation_b;
    QAction *actionDo_translation_c;
    QAction *actionDo_translation_d;
    QAction *actionG;
    QAction *actionG2;
    QAction *actionG4;
    QAction *actionG8;
    QAction *actionG16;
    QWidget *centralwidget;
    QGridLayout *gridLayout;
    QGraphicsView *graphicsView;
    QStatusBar *statusbar;
    QToolBar *fileToolBar;
    QToolBar *toolBar;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuEdit;
    QMenu *menuTools;

    void setupUi(QMainWindow *Periodic_4_hyperbolic_Delaunay_triangulation_2)
    {
        if (Periodic_4_hyperbolic_Delaunay_triangulation_2->objectName().isEmpty())
            Periodic_4_hyperbolic_Delaunay_triangulation_2->setObjectName(QString::fromUtf8("Periodic_4_hyperbolic_Delaunay_triangulation_2"));
        Periodic_4_hyperbolic_Delaunay_triangulation_2->resize(800, 600);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/cgal/logos/cgal_icon"), QSize(), QIcon::Normal, QIcon::Off);
        Periodic_4_hyperbolic_Delaunay_triangulation_2->setWindowIcon(icon);
        actionAbout = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionAbout->setObjectName(QString::fromUtf8("actionAbout"));
        actionAboutCGAL = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionAboutCGAL->setObjectName(QString::fromUtf8("actionAboutCGAL"));
        actionQuit = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionQuit->setObjectName(QString::fromUtf8("actionQuit"));
        actionInsertRandomPoints = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionInsertRandomPoints->setObjectName(QString::fromUtf8("actionInsertRandomPoints"));
        actionInsertOrigin = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionInsertOrigin->setObjectName(QString::fromUtf8("actionInsertOrigin"));
        actionInsertDummyPoints = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionInsertDummyPoints->setObjectName(QString::fromUtf8("actionInsertDummyPoints"));
        actionModifyDepth = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionModifyDepth->setObjectName(QString::fromUtf8("actionModifyDepth"));
        actionMovingPoint = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionMovingPoint->setObjectName(QString::fromUtf8("actionMovingPoint"));
        actionMovingPoint->setCheckable(true);
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/cgal/Actions/icons/moving_point.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionMovingPoint->setIcon(icon1);
        actionInsertPoint = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionInsertPoint->setObjectName(QString::fromUtf8("actionInsertPoint"));
        actionInsertPoint->setCheckable(true);
        actionInsertPoint->setChecked(false);
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/cgal/Input/inputPoint.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionInsertPoint->setIcon(icon2);
        actionClear = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionClear->setObjectName(QString::fromUtf8("actionClear"));
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/cgal/fileToolbar/fileNew.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionClear->setIcon(icon3);
        actionShowVoronoi = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionShowVoronoi->setObjectName(QString::fromUtf8("actionShowVoronoi"));
        actionShowVoronoi->setCheckable(true);
        actionShowVoronoi->setChecked(false);
        QIcon icon4;
        icon4.addFile(QString::fromUtf8(":/cgal/Triangulation_2/Voronoi_diagram_2.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionShowVoronoi->setIcon(icon4);
        actionShowDelaunay = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionShowDelaunay->setObjectName(QString::fromUtf8("actionShowDelaunay"));
        actionShowDelaunay->setCheckable(true);
        QIcon icon5;
        icon5.addFile(QString::fromUtf8(":/cgal/Actions/icons/triangulation.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionShowDelaunay->setIcon(icon5);
        actionLoadPoints = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionLoadPoints->setObjectName(QString::fromUtf8("actionLoadPoints"));
        QIcon icon6;
        icon6.addFile(QString::fromUtf8(":/cgal/fileToolbar/fileOpen.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLoadPoints->setIcon(icon6);
        actionSavePoints = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionSavePoints->setObjectName(QString::fromUtf8("actionSavePoints"));
        QIcon icon7;
        icon7.addFile(QString::fromUtf8(":/cgal/fileToolbar/fileSave.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSavePoints->setIcon(icon7);
        actionCircumcenter = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionCircumcenter->setObjectName(QString::fromUtf8("actionCircumcenter"));
        actionCircumcenter->setCheckable(true);
        QIcon icon8;
        icon8.addFile(QString::fromUtf8(":/cgal/Actions/icons/circumcenter.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionCircumcenter->setIcon(icon8);
        actionRecenter = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionRecenter->setObjectName(QString::fromUtf8("actionRecenter"));
        QIcon icon9;
        icon9.addFile(QString::fromUtf8(":/cgal/Input/zoom-best-fit"), QSize(), QIcon::Normal, QIcon::Off);
        actionRecenter->setIcon(icon9);
        actionShowConflictZone = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionShowConflictZone->setObjectName(QString::fromUtf8("actionShowConflictZone"));
        actionShowConflictZone->setCheckable(true);
        QIcon icon10;
        icon10.addFile(QString::fromUtf8(":/cgal/Actions/icons/conflict_zone.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionShowConflictZone->setIcon(icon10);
        actionDo_translation_a = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionDo_translation_a->setObjectName(QString::fromUtf8("actionDo_translation_a"));
        actionDo_translation_a->setCheckable(true);
        QIcon icon11;
        icon11.addFile(QString::fromUtf8(":/cgal/Actions/icons/a.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionDo_translation_a->setIcon(icon11);
        actionDo_translation_b = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionDo_translation_b->setObjectName(QString::fromUtf8("actionDo_translation_b"));
        actionDo_translation_b->setCheckable(true);
        QIcon icon12;
        icon12.addFile(QString::fromUtf8(":/cgal/Actions/icons/b.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionDo_translation_b->setIcon(icon12);
        actionDo_translation_c = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionDo_translation_c->setObjectName(QString::fromUtf8("actionDo_translation_c"));
        actionDo_translation_c->setCheckable(true);
        QIcon icon13;
        icon13.addFile(QString::fromUtf8(":/cgal/Actions/icons/c.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionDo_translation_c->setIcon(icon13);
        actionDo_translation_d = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionDo_translation_d->setObjectName(QString::fromUtf8("actionDo_translation_d"));
        actionDo_translation_d->setCheckable(true);
        QIcon icon14;
        icon14.addFile(QString::fromUtf8(":/cgal/Actions/icons/d.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionDo_translation_d->setIcon(icon14);
        actionG = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionG->setObjectName(QString::fromUtf8("actionG"));
        actionG->setCheckable(true);
        QIcon icon15;
        icon15.addFile(QString::fromUtf8(":/cgal/Actions/icons/G.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionG->setIcon(icon15);
        actionG2 = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionG2->setObjectName(QString::fromUtf8("actionG2"));
        actionG2->setCheckable(true);
        QIcon icon16;
        icon16.addFile(QString::fromUtf8(":/cgal/Actions/icons/G2.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionG2->setIcon(icon16);
        actionG4 = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionG4->setObjectName(QString::fromUtf8("actionG4"));
        actionG4->setCheckable(true);
        QIcon icon17;
        icon17.addFile(QString::fromUtf8(":/cgal/Actions/icons/G4.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionG4->setIcon(icon17);
        actionG8 = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionG8->setObjectName(QString::fromUtf8("actionG8"));
        actionG8->setCheckable(true);
        QIcon icon18;
        icon18.addFile(QString::fromUtf8(":/cgal/Actions/icons/G8.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionG8->setIcon(icon18);
        actionG16 = new QAction(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        actionG16->setObjectName(QString::fromUtf8("actionG16"));
        actionG16->setCheckable(true);
        QIcon icon19;
        icon19.addFile(QString::fromUtf8(":/cgal/Actions/icons/G16.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionG16->setIcon(icon19);
        centralwidget = new QWidget(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        gridLayout = new QGridLayout(centralwidget);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        graphicsView = new QGraphicsView(centralwidget);
        graphicsView->setObjectName(QString::fromUtf8("graphicsView"));
        graphicsView->setFocusPolicy(Qt::StrongFocus);
        graphicsView->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
        graphicsView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
        graphicsView->setTransformationAnchor(QGraphicsView::NoAnchor);

        gridLayout->addWidget(graphicsView, 0, 0, 1, 1);

        Periodic_4_hyperbolic_Delaunay_triangulation_2->setCentralWidget(centralwidget);
        statusbar = new QStatusBar(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        Periodic_4_hyperbolic_Delaunay_triangulation_2->setStatusBar(statusbar);
        fileToolBar = new QToolBar(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        fileToolBar->setObjectName(QString::fromUtf8("fileToolBar"));
        Periodic_4_hyperbolic_Delaunay_triangulation_2->addToolBar(Qt::TopToolBarArea, fileToolBar);
        toolBar = new QToolBar(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        toolBar->setObjectName(QString::fromUtf8("toolBar"));
        Periodic_4_hyperbolic_Delaunay_triangulation_2->addToolBar(Qt::TopToolBarArea, toolBar);
        menubar = new QMenuBar(Periodic_4_hyperbolic_Delaunay_triangulation_2);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 800, 22));
        menuFile = new QMenu(menubar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        menuEdit = new QMenu(menubar);
        menuEdit->setObjectName(QString::fromUtf8("menuEdit"));
        menuTools = new QMenu(menubar);
        menuTools->setObjectName(QString::fromUtf8("menuTools"));
        Periodic_4_hyperbolic_Delaunay_triangulation_2->setMenuBar(menubar);

        fileToolBar->addAction(actionClear);
        fileToolBar->addAction(actionLoadPoints);
        fileToolBar->addAction(actionSavePoints);
        toolBar->addAction(actionInsertPoint);
        toolBar->addAction(actionDo_translation_a);
        toolBar->addAction(actionDo_translation_b);
        toolBar->addAction(actionDo_translation_c);
        toolBar->addAction(actionDo_translation_d);
        toolBar->addAction(actionMovingPoint);
        toolBar->addAction(actionCircumcenter);
        toolBar->addAction(actionShowConflictZone);
        toolBar->addSeparator();
        toolBar->addAction(actionShowDelaunay);
        toolBar->addAction(actionShowVoronoi);
        toolBar->addSeparator();
        toolBar->addAction(actionRecenter);
        menubar->addAction(menuFile->menuAction());
        menubar->addAction(menuEdit->menuAction());
        menubar->addAction(menuTools->menuAction());
        menuFile->addSeparator();
        menuFile->addAction(actionClear);
        menuFile->addAction(actionLoadPoints);
        menuFile->addAction(actionSavePoints);
        menuFile->addSeparator();
        menuFile->addAction(actionQuit);
        menuEdit->addAction(actionInsertRandomPoints);
        menuEdit->addAction(actionInsertOrigin);
        menuEdit->addAction(actionInsertDummyPoints);
        menuEdit->addAction(actionModifyDepth);
        menuTools->addAction(actionInsertPoint);
        menuTools->addAction(actionDo_translation_a);
        menuTools->addAction(actionDo_translation_b);
        menuTools->addAction(actionDo_translation_c);
        menuTools->addAction(actionDo_translation_d);
        menuTools->addAction(actionG);
        menuTools->addAction(actionG2);
        menuTools->addAction(actionG4);
        menuTools->addAction(actionG8);
        menuTools->addAction(actionG16);
        menuTools->addAction(actionMovingPoint);
        menuTools->addAction(actionCircumcenter);
        menuTools->addAction(actionShowConflictZone);
        menuTools->addSeparator();
        menuTools->addAction(actionShowDelaunay);
        menuTools->addAction(actionShowVoronoi);
        menuTools->addSeparator();
        menuTools->addAction(actionRecenter);

        retranslateUi(Periodic_4_hyperbolic_Delaunay_triangulation_2);

        QMetaObject::connectSlotsByName(Periodic_4_hyperbolic_Delaunay_triangulation_2);
    } // setupUi

    void retranslateUi(QMainWindow *Periodic_4_hyperbolic_Delaunay_triangulation_2)
    {
        Periodic_4_hyperbolic_Delaunay_triangulation_2->setWindowTitle(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "CGAL Periodic hyperbolic Delaunay Triangulation", 0));
        actionAbout->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "&About", 0));
        actionAboutCGAL->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "About &CGAL", 0));
        actionQuit->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "&Quit", 0));
        actionQuit->setShortcut(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Ctrl+Q", 0));
        actionInsertRandomPoints->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "&Insert random points", 0));
        actionInsertRandomPoints->setShortcut(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Ctrl+I", 0));
        actionInsertOrigin->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Insert Origin", 0));
        actionInsertOrigin->setShortcut(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Ctrl+Alt+O", 0));
        actionInsertDummyPoints->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Insert dummy points", 0));
        actionInsertDummyPoints->setShortcut(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Ctrl+Alt+D", 0));
        actionModifyDepth->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "&Modify recursion depth", 0));
        actionModifyDepth->setShortcut(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Ctrl+M", 0));
        actionMovingPoint->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "&Simulate insertion", 0));
#ifndef QT_NO_TOOLTIP
        actionMovingPoint->setToolTip(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Simulate Insertion", "The comment"));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_STATUSTIP
        actionMovingPoint->setStatusTip(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Move mouse with left button pressed to see where point would be inserted", "and its comment"));
#endif // QT_NO_STATUSTIP
#ifndef QT_NO_WHATSTHIS
        actionMovingPoint->setWhatsThis(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "whats this", "what"));
#endif // QT_NO_WHATSTHIS
        actionInsertPoint->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "&Insert Point", 0));
#ifndef QT_NO_TOOLTIP
        actionInsertPoint->setToolTip(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Insert Point", 0));
#endif // QT_NO_TOOLTIP
#ifndef QT_NO_STATUSTIP
        actionInsertPoint->setStatusTip(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Left: Insert vtx", 0));
#endif // QT_NO_STATUSTIP
        actionClear->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "&Clear", 0));
        actionClear->setShortcut(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Ctrl+C", 0));
        actionShowVoronoi->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Show &Voronoi Diagram", 0));
        actionShowVoronoi->setShortcut(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Ctrl+V", 0));
        actionShowDelaunay->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Show &Delaunay Triangulation", 0));
        actionShowDelaunay->setShortcut(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Ctrl+D", 0));
        actionLoadPoints->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "&Load Points...", 0));
        actionLoadPoints->setShortcut(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Ctrl+L", 0));
        actionSavePoints->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "&Save Points...", 0));
        actionSavePoints->setShortcut(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Ctrl+S", 0));
        actionCircumcenter->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "&Circumcenter", 0));
#ifndef QT_NO_TOOLTIP
        actionCircumcenter->setToolTip(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Draw circumcenter", 0));
#endif // QT_NO_TOOLTIP
        actionRecenter->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Re&center the viewport", 0));
        actionRecenter->setShortcut(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Ctrl+R", 0));
        actionShowConflictZone->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Show Conflict Zone", 0));
        actionDo_translation_a->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Do translation \"a\"", 0));
        actionDo_translation_b->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Do translation \"b\"", 0));
        actionDo_translation_c->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Do translation \"c\"", 0));
        actionDo_translation_d->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Do translation \"d\"", 0));
        actionG->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "G", 0));
        actionG2->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "G2", 0));
        actionG4->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "G4", 0));
        actionG8->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "G8", 0));
        actionG16->setText(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "G16", 0));
        fileToolBar->setWindowTitle(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "File Tools", 0));
        toolBar->setWindowTitle(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "Visualization Tools", 0));
        menuFile->setTitle(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "&File", 0));
        menuEdit->setTitle(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "&Edit", 0));
        menuTools->setTitle(QApplication::translate("Periodic_4_hyperbolic_Delaunay_triangulation_2", "&Tools", 0));
    } // retranslateUi

};

namespace Ui {
    class Periodic_4_hyperbolic_Delaunay_triangulation_2: public Ui_Periodic_4_hyperbolic_Delaunay_triangulation_2 {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PERIODIC_4_HYPERBOLIC_DELAUNAY_TRIANGULATION_H
