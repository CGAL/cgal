/****************************************************************************
** ToolApp meta object code from reading C++ file 'tool.h'
**
** Created: Thu May 22 16:06:12 2003
**      by: The Qt MOC ($Id$)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#if !defined(Q_MOC_OUTPUT_REVISION)
#define Q_MOC_OUTPUT_REVISION 9
#elif Q_MOC_OUTPUT_REVISION != 9
#error "Moc format conflict - please regenerate all moc files"
#endif

#include "tool.h"
#include <qmetaobject.h>
#include <qapplication.h>



const char *ToolApp::className() const
{
    return "ToolApp";
}

QMetaObject *ToolApp::metaObj = 0;

void ToolApp::initMetaObject()
{
    if ( metaObj )
	return;
    if ( qstrcmp(QMainWindow::className(), "QMainWindow") != 0 )
	badSuperclassWarning("ToolApp","QMainWindow");
    (void) staticMetaObject();
}

#ifndef QT_NO_TRANSLATION

QString ToolApp::tr(const char* s)
{
    return qApp->translate( "ToolApp", s, 0 );
}

QString ToolApp::tr(const char* s, const char * c)
{
    return qApp->translate( "ToolApp", s, c );
}

#endif // QT_NO_TRANSLATION

QMetaObject* ToolApp::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    (void) QMainWindow::staticMetaObject();
#ifndef QT_NO_PROPERTIES
#endif // QT_NO_PROPERTIES
    typedef void (ToolApp::*m1_t0)();
    typedef void (QObject::*om1_t0)();
    typedef void (ToolApp::*m1_t1)();
    typedef void (QObject::*om1_t1)();
    typedef void (ToolApp::*m1_t2)();
    typedef void (QObject::*om1_t2)();
    typedef void (ToolApp::*m1_t3)();
    typedef void (QObject::*om1_t3)();
    typedef void (ToolApp::*m1_t4)();
    typedef void (QObject::*om1_t4)();
    typedef void (ToolApp::*m1_t5)();
    typedef void (QObject::*om1_t5)();
    typedef void (ToolApp::*m1_t6)();
    typedef void (QObject::*om1_t6)();
    typedef void (ToolApp::*m1_t7)();
    typedef void (QObject::*om1_t7)();
    typedef void (ToolApp::*m1_t8)();
    typedef void (QObject::*om1_t8)();
    typedef void (ToolApp::*m1_t9)();
    typedef void (QObject::*om1_t9)();
    typedef void (ToolApp::*m1_t10)();
    typedef void (QObject::*om1_t10)();
    typedef void (ToolApp::*m1_t11)(bool);
    typedef void (QObject::*om1_t11)(bool);
    typedef void (ToolApp::*m1_t12)(bool);
    typedef void (QObject::*om1_t12)(bool);
    typedef void (ToolApp::*m1_t13)();
    typedef void (QObject::*om1_t13)();
    typedef void (ToolApp::*m1_t14)();
    typedef void (QObject::*om1_t14)();
    typedef void (ToolApp::*m1_t15)();
    typedef void (QObject::*om1_t15)();
    typedef void (ToolApp::*m1_t16)();
    typedef void (QObject::*om1_t16)();
    typedef void (ToolApp::*m1_t17)();
    typedef void (QObject::*om1_t17)();
    typedef void (ToolApp::*m1_t18)();
    typedef void (QObject::*om1_t18)();
    typedef void (ToolApp::*m1_t19)();
    typedef void (QObject::*om1_t19)();
    typedef void (ToolApp::*m1_t20)();
    typedef void (QObject::*om1_t20)();
    typedef void (ToolApp::*m1_t21)();
    typedef void (QObject::*om1_t21)();
    typedef void (ToolApp::*m1_t22)();
    typedef void (QObject::*om1_t22)();
    typedef void (ToolApp::*m1_t23)();
    typedef void (QObject::*om1_t23)();
    typedef void (ToolApp::*m1_t24)();
    typedef void (QObject::*om1_t24)();
    typedef void (ToolApp::*m1_t25)();
    typedef void (QObject::*om1_t25)();
    typedef void (ToolApp::*m1_t26)();
    typedef void (QObject::*om1_t26)();
    typedef void (ToolApp::*m1_t27)();
    typedef void (QObject::*om1_t27)();
    typedef void (ToolApp::*m1_t28)(const QString&);
    typedef void (QObject::*om1_t28)(const QString&);
    typedef void (ToolApp::*m1_t29)();
    typedef void (QObject::*om1_t29)();
    typedef void (ToolApp::*m1_t30)(int);
    typedef void (QObject::*om1_t30)(int);
    m1_t0 v1_0 = &ToolApp::slotFileNew;
    om1_t0 ov1_0 = (om1_t0)v1_0;
    m1_t1 v1_1 = &ToolApp::slotFileOpen;
    om1_t1 ov1_1 = (om1_t1)v1_1;
    m1_t2 v1_2 = &ToolApp::slotFileSave;
    om1_t2 ov1_2 = (om1_t2)v1_2;
    m1_t3 v1_3 = &ToolApp::slotFileSaveAs;
    om1_t3 ov1_3 = (om1_t3)v1_3;
    m1_t4 v1_4 = &ToolApp::slotFileClose;
    om1_t4 ov1_4 = (om1_t4)v1_4;
    m1_t5 v1_5 = &ToolApp::slotFilePrint;
    om1_t5 ov1_5 = (om1_t5)v1_5;
    m1_t6 v1_6 = &ToolApp::slotFileQuit;
    om1_t6 ov1_6 = (om1_t6)v1_6;
    m1_t7 v1_7 = &ToolApp::slotEditUndo;
    om1_t7 ov1_7 = (om1_t7)v1_7;
    m1_t8 v1_8 = &ToolApp::slotEditCut;
    om1_t8 ov1_8 = (om1_t8)v1_8;
    m1_t9 v1_9 = &ToolApp::slotEditCopy;
    om1_t9 ov1_9 = (om1_t9)v1_9;
    m1_t10 v1_10 = &ToolApp::slotEditPaste;
    om1_t10 ov1_10 = (om1_t10)v1_10;
    m1_t11 v1_11 = &ToolApp::slotViewToolBar;
    om1_t11 ov1_11 = (om1_t11)v1_11;
    m1_t12 v1_12 = &ToolApp::slotViewStatusBar;
    om1_t12 ov1_12 = (om1_t12)v1_12;
    m1_t13 v1_13 = &ToolApp::slotViewCopyViewpoint;
    om1_t13 ov1_13 = (om1_t13)v1_13;
    m1_t14 v1_14 = &ToolApp::slotViewPasteViewpoint;
    om1_t14 ov1_14 = (om1_t14)v1_14;
    m1_t15 v1_15 = &ToolApp::slotRenderWhite;
    om1_t15 ov1_15 = (om1_t15)v1_15;
    m1_t16 v1_16 = &ToolApp::slotRenderBlack;
    om1_t16 ov1_16 = (om1_t16)v1_16;
    m1_t17 v1_17 = &ToolApp::slotRenderWireframe;
    om1_t17 ov1_17 = (om1_t17)v1_17;
    m1_t18 v1_18 = &ToolApp::slotRenderEdges;
    om1_t18 ov1_18 = (om1_t18)v1_18;
    m1_t19 v1_19 = &ToolApp::slotRenderFill;
    om1_t19 ov1_19 = (om1_t19)v1_19;
    m1_t20 v1_20 = &ToolApp::slotRenderSmooth;
    om1_t20 ov1_20 = (om1_t20)v1_20;
    m1_t21 v1_21 = &ToolApp::slotRenderCulling;
    om1_t21 ov1_21 = (om1_t21)v1_21;
    m1_t22 v1_22 = &ToolApp::slotRenderLighting;
    om1_t22 ov1_22 = (om1_t22)v1_22;
    m1_t23 v1_23 = &ToolApp::slotRenderVertex;
    om1_t23 ov1_23 = (om1_t23)v1_23;
    m1_t24 v1_24 = &ToolApp::slotRenderSuperimposing;
    om1_t24 ov1_24 = (om1_t24)v1_24;
    m1_t25 v1_25 = &ToolApp::slotRenderAntialiasing;
    om1_t25 ov1_25 = (om1_t25)v1_25;
    m1_t26 v1_26 = &ToolApp::slotWindowNewWindow;
    om1_t26 ov1_26 = (om1_t26)v1_26;
    m1_t27 v1_27 = &ToolApp::slotHelpAbout;
    om1_t27 ov1_27 = (om1_t27)v1_27;
    m1_t28 v1_28 = &ToolApp::slotStatusHelpMsg;
    om1_t28 ov1_28 = (om1_t28)v1_28;
    m1_t29 v1_29 = &ToolApp::windowMenuAboutToShow;
    om1_t29 ov1_29 = (om1_t29)v1_29;
    m1_t30 v1_30 = &ToolApp::windowMenuActivated;
    om1_t30 ov1_30 = (om1_t30)v1_30;
    QMetaData *slot_tbl = QMetaObject::new_metadata(31);
    QMetaData::Access *slot_tbl_access = QMetaObject::new_metaaccess(31);
    slot_tbl[0].name = "slotFileNew()";
    slot_tbl[0].ptr = (QMember)ov1_0;
    slot_tbl_access[0] = QMetaData::Private;
    slot_tbl[1].name = "slotFileOpen()";
    slot_tbl[1].ptr = (QMember)ov1_1;
    slot_tbl_access[1] = QMetaData::Private;
    slot_tbl[2].name = "slotFileSave()";
    slot_tbl[2].ptr = (QMember)ov1_2;
    slot_tbl_access[2] = QMetaData::Private;
    slot_tbl[3].name = "slotFileSaveAs()";
    slot_tbl[3].ptr = (QMember)ov1_3;
    slot_tbl_access[3] = QMetaData::Private;
    slot_tbl[4].name = "slotFileClose()";
    slot_tbl[4].ptr = (QMember)ov1_4;
    slot_tbl_access[4] = QMetaData::Private;
    slot_tbl[5].name = "slotFilePrint()";
    slot_tbl[5].ptr = (QMember)ov1_5;
    slot_tbl_access[5] = QMetaData::Private;
    slot_tbl[6].name = "slotFileQuit()";
    slot_tbl[6].ptr = (QMember)ov1_6;
    slot_tbl_access[6] = QMetaData::Private;
    slot_tbl[7].name = "slotEditUndo()";
    slot_tbl[7].ptr = (QMember)ov1_7;
    slot_tbl_access[7] = QMetaData::Private;
    slot_tbl[8].name = "slotEditCut()";
    slot_tbl[8].ptr = (QMember)ov1_8;
    slot_tbl_access[8] = QMetaData::Private;
    slot_tbl[9].name = "slotEditCopy()";
    slot_tbl[9].ptr = (QMember)ov1_9;
    slot_tbl_access[9] = QMetaData::Private;
    slot_tbl[10].name = "slotEditPaste()";
    slot_tbl[10].ptr = (QMember)ov1_10;
    slot_tbl_access[10] = QMetaData::Private;
    slot_tbl[11].name = "slotViewToolBar(bool)";
    slot_tbl[11].ptr = (QMember)ov1_11;
    slot_tbl_access[11] = QMetaData::Private;
    slot_tbl[12].name = "slotViewStatusBar(bool)";
    slot_tbl[12].ptr = (QMember)ov1_12;
    slot_tbl_access[12] = QMetaData::Private;
    slot_tbl[13].name = "slotViewCopyViewpoint()";
    slot_tbl[13].ptr = (QMember)ov1_13;
    slot_tbl_access[13] = QMetaData::Private;
    slot_tbl[14].name = "slotViewPasteViewpoint()";
    slot_tbl[14].ptr = (QMember)ov1_14;
    slot_tbl_access[14] = QMetaData::Private;
    slot_tbl[15].name = "slotRenderWhite()";
    slot_tbl[15].ptr = (QMember)ov1_15;
    slot_tbl_access[15] = QMetaData::Private;
    slot_tbl[16].name = "slotRenderBlack()";
    slot_tbl[16].ptr = (QMember)ov1_16;
    slot_tbl_access[16] = QMetaData::Private;
    slot_tbl[17].name = "slotRenderWireframe()";
    slot_tbl[17].ptr = (QMember)ov1_17;
    slot_tbl_access[17] = QMetaData::Private;
    slot_tbl[18].name = "slotRenderEdges()";
    slot_tbl[18].ptr = (QMember)ov1_18;
    slot_tbl_access[18] = QMetaData::Private;
    slot_tbl[19].name = "slotRenderFill()";
    slot_tbl[19].ptr = (QMember)ov1_19;
    slot_tbl_access[19] = QMetaData::Private;
    slot_tbl[20].name = "slotRenderSmooth()";
    slot_tbl[20].ptr = (QMember)ov1_20;
    slot_tbl_access[20] = QMetaData::Private;
    slot_tbl[21].name = "slotRenderCulling()";
    slot_tbl[21].ptr = (QMember)ov1_21;
    slot_tbl_access[21] = QMetaData::Private;
    slot_tbl[22].name = "slotRenderLighting()";
    slot_tbl[22].ptr = (QMember)ov1_22;
    slot_tbl_access[22] = QMetaData::Private;
    slot_tbl[23].name = "slotRenderVertex()";
    slot_tbl[23].ptr = (QMember)ov1_23;
    slot_tbl_access[23] = QMetaData::Private;
    slot_tbl[24].name = "slotRenderSuperimposing()";
    slot_tbl[24].ptr = (QMember)ov1_24;
    slot_tbl_access[24] = QMetaData::Private;
    slot_tbl[25].name = "slotRenderAntialiasing()";
    slot_tbl[25].ptr = (QMember)ov1_25;
    slot_tbl_access[25] = QMetaData::Private;
    slot_tbl[26].name = "slotWindowNewWindow()";
    slot_tbl[26].ptr = (QMember)ov1_26;
    slot_tbl_access[26] = QMetaData::Private;
    slot_tbl[27].name = "slotHelpAbout()";
    slot_tbl[27].ptr = (QMember)ov1_27;
    slot_tbl_access[27] = QMetaData::Private;
    slot_tbl[28].name = "slotStatusHelpMsg(const QString&)";
    slot_tbl[28].ptr = (QMember)ov1_28;
    slot_tbl_access[28] = QMetaData::Private;
    slot_tbl[29].name = "windowMenuAboutToShow()";
    slot_tbl[29].ptr = (QMember)ov1_29;
    slot_tbl_access[29] = QMetaData::Private;
    slot_tbl[30].name = "windowMenuActivated(int)";
    slot_tbl[30].ptr = (QMember)ov1_30;
    slot_tbl_access[30] = QMetaData::Private;
    metaObj = QMetaObject::new_metaobject(
	"ToolApp", "QMainWindow",
	slot_tbl, 31,
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
