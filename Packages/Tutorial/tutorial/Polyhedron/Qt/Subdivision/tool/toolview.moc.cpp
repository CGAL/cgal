/****************************************************************************
** ToolView meta object code from reading C++ file 'toolview.h'
**
** Created: Thu May 22 16:04:07 2003
**      by: The Qt MOC ($Id$)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#if !defined(Q_MOC_OUTPUT_REVISION)
#define Q_MOC_OUTPUT_REVISION 9
#elif Q_MOC_OUTPUT_REVISION != 9
#error "Moc format conflict - please regenerate all moc files"
#endif

#include "toolview.h"
#include <qmetaobject.h>
#include <qapplication.h>



const char *ToolView::className() const
{
    return "ToolView";
}

QMetaObject *ToolView::metaObj = 0;

void ToolView::initMetaObject()
{
    if ( metaObj )
	return;
    if ( qstrcmp(QGLView::className(), "QGLView") != 0 )
	badSuperclassWarning("ToolView","QGLView");
    (void) staticMetaObject();
}

#ifndef QT_NO_TRANSLATION

QString ToolView::tr(const char* s)
{
    return qApp->translate( "ToolView", s, 0 );
}

QString ToolView::tr(const char* s, const char * c)
{
    return qApp->translate( "ToolView", s, c );
}

#endif // QT_NO_TRANSLATION

QMetaObject* ToolView::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    (void) QGLView::staticMetaObject();
#ifndef QT_NO_PROPERTIES
#endif // QT_NO_PROPERTIES
    QMetaData::Access *slot_tbl_access = 0;
    metaObj = QMetaObject::new_metaobject(
	"ToolView", "QGLView",
	0, 0,
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
