/****************************************************************************
** TrFrame meta object code from reading C++ file 'viewer.h'
**
** Created: Wed Nov 7 10:13:25 2001
**      by: The Qt MOC ($Id$)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#if !defined(Q_MOC_OUTPUT_REVISION)
#define Q_MOC_OUTPUT_REVISION 9
#elif Q_MOC_OUTPUT_REVISION != 9
#error "Moc format conflict - please regenerate all moc files"
#endif

#include "viewer.h"
#include <qmetaobject.h>
#include <qapplication.h>



const char *TrFrame::className() const
{
    return "TrFrame";
}

QMetaObject *TrFrame::metaObj = 0;

void TrFrame::initMetaObject()
{
    if ( metaObj )
	return;
    if ( qstrcmp(QMainWindow::className(), "QMainWindow") != 0 )
	badSuperclassWarning("TrFrame","QMainWindow");
    (void) staticMetaObject();
}

#ifndef QT_NO_TRANSLATION

QString TrFrame::tr(const char* s)
{
    return qApp->translate( "TrFrame", s, 0 );
}

QString TrFrame::tr(const char* s, const char * c)
{
    return qApp->translate( "TrFrame", s, c );
}

#endif // QT_NO_TRANSLATION

QMetaObject* TrFrame::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    (void) QMainWindow::staticMetaObject();
#ifndef QT_NO_PROPERTIES
#endif // QT_NO_PROPERTIES
    typedef void (TrFrame::*m1_t0)();
    typedef void (QObject::*om1_t0)();
    typedef void (TrFrame::*m1_t1)();
    typedef void (QObject::*om1_t1)();
    typedef void (TrFrame::*m1_t2)();
    typedef void (QObject::*om1_t2)();
    typedef void (TrFrame::*m1_t3)();
    typedef void (QObject::*om1_t3)();
    typedef void (TrFrame::*m1_t4)();
    typedef void (QObject::*om1_t4)();
    m1_t0 v1_0 = &TrFrame::clearTriangulation;
    om1_t0 ov1_0 = (om1_t0)v1_0;
    m1_t1 v1_1 = &TrFrame::openTriangulation;
    om1_t1 ov1_1 = (om1_t1)v1_1;
    m1_t2 v1_2 = &TrFrame::saveTriangulation;
    om1_t2 ov1_2 = (om1_t2)v1_2;
    m1_t3 v1_3 = &TrFrame::onChangeSizes;
    om1_t3 ov1_3 = (om1_t3)v1_3;
    m1_t4 v1_4 = &TrFrame::mesh;
    om1_t4 ov1_4 = (om1_t4)v1_4;
    QMetaData *slot_tbl = QMetaObject::new_metadata(5);
    QMetaData::Access *slot_tbl_access = QMetaObject::new_metaaccess(5);
    slot_tbl[0].name = "clearTriangulation()";
    slot_tbl[0].ptr = (QMember)ov1_0;
    slot_tbl_access[0] = QMetaData::Private;
    slot_tbl[1].name = "openTriangulation()";
    slot_tbl[1].ptr = (QMember)ov1_1;
    slot_tbl_access[1] = QMetaData::Private;
    slot_tbl[2].name = "saveTriangulation()";
    slot_tbl[2].ptr = (QMember)ov1_2;
    slot_tbl_access[2] = QMetaData::Private;
    slot_tbl[3].name = "onChangeSizes()";
    slot_tbl[3].ptr = (QMember)ov1_3;
    slot_tbl_access[3] = QMetaData::Private;
    slot_tbl[4].name = "mesh()";
    slot_tbl[4].ptr = (QMember)ov1_4;
    slot_tbl_access[4] = QMetaData::Private;
    metaObj = QMetaObject::new_metaobject(
	"TrFrame", "QMainWindow",
	slot_tbl, 5,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    metaObj->set_slot_access( slot_tbl_access );
#ifndef QT_NO_PROPERTIES
#endif // QT_NO_PROPERTIES
    return metaObj;
}


const char *TrViewer::className() const
{
    return "TrViewer";
}

QMetaObject *TrViewer::metaObj = 0;

void TrViewer::initMetaObject()
{
    if ( metaObj )
	return;
    if ( qstrcmp(QWidget::className(), "QWidget") != 0 )
	badSuperclassWarning("TrViewer","QWidget");
    (void) staticMetaObject();
}

#ifndef QT_NO_TRANSLATION

QString TrViewer::tr(const char* s)
{
    return qApp->translate( "TrViewer", s, 0 );
}

QString TrViewer::tr(const char* s, const char * c)
{
    return qApp->translate( "TrViewer", s, c );
}

#endif // QT_NO_TRANSLATION

QMetaObject* TrViewer::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    (void) QWidget::staticMetaObject();
#ifndef QT_NO_PROPERTIES
#endif // QT_NO_PROPERTIES
    typedef void (TrViewer::*m1_t0)();
    typedef void (QObject::*om1_t0)();
    typedef void (TrViewer::*m1_t1)();
    typedef void (QObject::*om1_t1)();
    m1_t0 v1_0 = &TrViewer::onStep;
    om1_t0 ov1_0 = (om1_t0)v1_0;
    m1_t1 v1_1 = &TrViewer::onMesh;
    om1_t1 ov1_1 = (om1_t1)v1_1;
    QMetaData *slot_tbl = QMetaObject::new_metadata(2);
    QMetaData::Access *slot_tbl_access = QMetaObject::new_metaaccess(2);
    slot_tbl[0].name = "onStep()";
    slot_tbl[0].ptr = (QMember)ov1_0;
    slot_tbl_access[0] = QMetaData::Public;
    slot_tbl[1].name = "onMesh()";
    slot_tbl[1].ptr = (QMember)ov1_1;
    slot_tbl_access[1] = QMetaData::Public;
    metaObj = QMetaObject::new_metaobject(
	"TrViewer", "QWidget",
	slot_tbl, 2,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    metaObj->set_slot_access( slot_tbl_access );
#ifndef QT_NO_PROPERTIES
#endif // QT_NO_PROPERTIES
    return metaObj;
}
