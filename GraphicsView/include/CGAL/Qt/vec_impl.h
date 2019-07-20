/****************************************************************************

 Copyright (c) 2018  GeometryFactory Sarl (France).
 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of a fork of the QGLViewer library version 2.7.0.
 http://www.libqglviewer.com - contact@libqglviewer.com

 This file may be used under the terms of the GNU General Public License 
 version 3.0 as published by the Free Software Foundation and
 appearing in the LICENSE file included in the packaging of this file.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

*****************************************************************************/
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline

#include <CGAL/license/GraphicsView.h>

#else
#define CGAL_INLINE_FUNCTION
#endif

#include <CGAL/Qt/vec.h>
#include <CGAL/Qt/domUtils.h>

// Most of the methods are declared inline in vec.h

namespace CGAL{
namespace qglviewer{


/*! Projects the Vec on the axis of direction \p direction that passes through
the origin.

\p direction does not need to be normalized (but must be non null). */
CGAL_INLINE_FUNCTION
void Vec::projectOnAxis(const Vec &direction) {
#ifndef QT_NO_DEBUG
  if (direction.squaredNorm() < 1.0E-10)
    qWarning("Vec::projectOnAxis: axis direction is not normalized (norm=%f).",
             direction.norm());
#endif

  *this = (((*this) * direction) / direction.squaredNorm()) * direction;
}

/*! Projects the Vec on the plane whose normal is \p normal that passes through
the origin.

\p normal does not need to be normalized (but must be non null). */
CGAL_INLINE_FUNCTION
void Vec::projectOnPlane(const Vec &normal) {
#ifndef QT_NO_DEBUG
  if (normal.squaredNorm() < 1.0E-10)
    qWarning("Vec::projectOnPlane: plane normal is not normalized (norm=%f).",
             normal.norm());
#endif

  *this -= (((*this) * normal) / normal.squaredNorm()) * normal;
}

/*! Returns a Vec orthogonal to the Vec. Its norm() depends on the Vec, but is
 zero only for a null Vec. Note that the function that associates an
 orthogonalVec() to a Vec is not continous. */
CGAL_INLINE_FUNCTION
Vec Vec::orthogonalVec() const {
  // Find smallest component. Keep equal case for null values.
  if ((fabs(y) >= 0.9 * fabs(x)) && (fabs(z) >= 0.9 * fabs(x)))
    return Vec(0.0, -z, y);
  else if ((fabs(x) >= 0.9 * fabs(y)) && (fabs(z) >= 0.9 * fabs(y)))
    return Vec(-z, 0.0, x);
  else
    return Vec(-y, x, 0.0);
}

/*! Constructs a Vec from a \c QDomElement representing an XML code of the form
 \code< anyTagName x=".." y=".." z=".." />\endcode

If one of these attributes is missing or is not a number, a warning is displayed
and the associated value is set to 0.0.

See also domElement() and initFromDOMElement(). */
CGAL_INLINE_FUNCTION
Vec::Vec(const QDomElement &element) {
  QStringList attribute;
  attribute << "x"
            << "y"
            << "z";
  for (int i = 0; i < attribute.size(); ++i)
#ifdef QGLVIEWER_UNION_NOT_SUPPORTED
    this->operator[](i) = DomUtils::qrealFromDom(element, attribute[i], 0.0);
#else
    v_[i] = DomUtils::qrealFromDom(element, attribute[i], 0.0);
#endif
}

/*! Returns an XML \c QDomElement that represents the Vec.

 \p name is the name of the QDomElement tag. \p doc is the \c QDomDocument
 factory used to create QDomElement.

 When output to a file, the resulting QDomElement will look like:
 \code
 <name x=".." y=".." z=".." />
 \endcode

 Use initFromDOMElement() to restore the Vec state from the resulting \c
 QDomElement. See also the Vec(const QDomElement&) constructor.

 Here is complete example that creates a QDomDocument and saves it into a file:
 \code
 Vec sunPos;
 QDomDocument document("myDocument");
 QDomElement sunElement = document.createElement("Sun");
 document.appendChild(sunElement);
 sunElement.setAttribute("brightness", sunBrightness());
 sunElement.appendChild(sunPos.domElement("sunPosition", document));
 // Other additions to the document hierarchy...

 // Save doc document
 QFile f("myFile.xml");
 if (f.open(IO_WriteOnly))
 {
   QTextStream out(&f);
   document.save(out, 2);
   f.close();
 }
 \endcode

CGAL_INLINE_FUNCTION
 See also Quaternion::domElement(), Frame::domElement(), Camera::domElement()...
 */
CGAL_INLINE_FUNCTION
QDomElement Vec::domElement(const QString &name, QDomDocument &document) const {
  QDomElement de = document.createElement(name);
  de.setAttribute("x", QString::number(x));
  de.setAttribute("y", QString::number(y));
  de.setAttribute("z", QString::number(z));
  return de;
}

/*! Restores the Vec state from a \c QDomElement created by domElement().

 The \c QDomElement should contain \c x, \c y and \c z attributes. If one of
 these attributes is missing or is not a number, a warning is displayed and the
 associated value is set to 0.0.

 To restore the Vec state from an xml file, use:
 \code
 // Load DOM from file
 QDomDocument doc;
 QFile f("myFile.xml");
 if (f.open(IO_ReadOnly))
 {
   doc.setContent(&f);
   f.close();
 }
 // Parse the DOM tree and initialize
 QDomElement main=doc.documentElement();
 myVec.initFromDOMElement(main);
 \endcode

 See also the Vec(const QDomElement&) constructor. */
CGAL_INLINE_FUNCTION
void Vec::initFromDOMElement(const QDomElement &element) {
  const Vec v(element);
  *this = v;
}

CGAL_INLINE_FUNCTION
std::ostream &operator<<(std::ostream &o, const Vec &v) {
  return o << v.x << '\t' << v.y << '\t' << v.z;
}  

}} // namespace CGAL::qglviewer
