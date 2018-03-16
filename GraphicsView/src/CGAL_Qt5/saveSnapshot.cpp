/****************************************************************************

 Copyright (C) 2002-2014 Gilles Debunne. All rights reserved.

 This file is part of the QGLViewer library version 2.7.0.

 http://www.libqglviewer.com - contact@libqglviewer.com

 This file may be used under the terms of the GNU General Public License 
 versions 2.0 or 3.0 as published by the Free Software Foundation and
 appearing in the LICENSE file included in the packaging of this file.
 In addition, as a special exception, Gilles Debunne gives you certain 
 additional rights, described in the file GPL_EXCEPTION in this package.

 libQGLViewer uses dual licensing. Commercial/proprietary software must
 purchase a libQGLViewer Commercial License.

 This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

*****************************************************************************/

#include <CGAL/Qt/qglviewer.h>

// Output format list
#include <QImageWriter>

#include <qapplication.h>
#include <qcursor.h>
#include <qfiledialog.h>
#include <qfileinfo.h>
#include <qinputdialog.h>
#include <qmap.h>
#include <qmessagebox.h>
#include <qprogressdialog.h>

using namespace std;

////// Static global variables - local to this file //////
// List of available output file formats, formatted for QFileDialog.
static QString formats;
// Converts QFileDialog resulting format to Qt snapshotFormat.
static QMap<QString, QString> Qtformat;
// Converts Qt snapshotFormat to QFileDialog menu string.
static QMap<QString, QString> FDFormatString;
// Converts snapshotFormat to file extension
static QMap<QString, QString> extension;

/*! Sets snapshotFileName(). */
void QGLViewer::setSnapshotFileName(const QString &name) {
  snapshotFileName_ = QFileInfo(name).absoluteFilePath();
}

#ifndef DOXYGEN
const QString &QGLViewer::snapshotFilename() const {
  qWarning("snapshotFilename is deprecated. Use snapshotFileName() (uppercase "
           "N) instead.");
  return snapshotFileName();
}
#endif

/*! Opens a dialog that displays the different available snapshot formats.

Then calls setSnapshotFormat() with the selected one (unless the user cancels).

Returns \c false if the user presses the Cancel button and \c true otherwise. */
bool QGLViewer::openSnapshotFormatDialog() {
  bool ok = false;
  QStringList list = formats.split(";;", QString::SkipEmptyParts);
  int current = list.indexOf(FDFormatString[snapshotFormat()]);
  QString format =
      QInputDialog::getItem(this, "Snapshot format", "Select a snapshot format",
                            list, current, false, &ok);
  if (ok)
    setSnapshotFormat(Qtformat[format]);
  return ok;
}

// Finds all available Qt output formats, so that they can be available in
// saveSnapshot dialog. Initialize snapshotFormat() to the first one.
void QGLViewer::initializeSnapshotFormats() {
  QList<QByteArray> list = QImageWriter::supportedImageFormats();
  QStringList formatList;
  for (int i = 0; i < list.size(); ++i)
    formatList << QString(list.at(i).toUpper());
//        qWarning("Available image formats: ");
//        QStringList::Iterator it = formatList.begin();
//        while( it != formatList.end() )
//  	      qWarning((*it++).);  QT4 change this. qWarning no longer accepts
//  QString


  // Check that the interesting formats are available and add them in "formats"
  // Unused formats: XPM XBM PBM PGM
  QStringList QtText, MenuText, Ext;
  QtText += "JPEG";
  MenuText += "JPEG (*.jpg)";
  Ext += "jpg";
  QtText += "PNG";
  MenuText += "PNG (*.png)";
  Ext += "png";
  QtText += "EPS";
  MenuText += "Encapsulated Postscript (*.eps)";
  Ext += "eps";
  QtText += "PS";
  MenuText += "Postscript (*.ps)";
  Ext += "ps";
  QtText += "PPM";
  MenuText += "24bit RGB Bitmap (*.ppm)";
  Ext += "ppm";
  QtText += "BMP";
  MenuText += "Windows Bitmap (*.bmp)";
  Ext += "bmp";
  QtText += "XFIG";
  MenuText += "XFig (*.fig)";
  Ext += "fig";

  QStringList::iterator itText = QtText.begin();
  QStringList::iterator itMenu = MenuText.begin();
  QStringList::iterator itExt = Ext.begin();

  while (itText != QtText.end()) {
    // QMessageBox::information(this, "Snapshot ", "Trying format\n"+(*itText));
    if (formatList.contains((*itText))) {
      // QMessageBox::information(this, "Snapshot ", "Recognized
      // format\n"+(*itText));
      if (formats.isEmpty())
        setSnapshotFormat(*itText);
      else
        formats += ";;";
      formats += (*itMenu);
      Qtformat[(*itMenu)] = (*itText);
      FDFormatString[(*itText)] = (*itMenu);
      extension[(*itText)] = (*itExt);
    }
    // Synchronize parsing
    itText++;
    itMenu++;
    itExt++;
  }
}

// Returns false if the user refused to use the fileName
static bool checkFileName(QString &fileName, QWidget *widget,
                          const QString &snapshotFormat) {
  if (fileName.isEmpty())
    return false;

  // Check that extension has been provided
  QFileInfo info(fileName);

  if (info.suffix().isEmpty()) {
    // No extension given. Silently add one
    if (fileName.right(1) != ".")
      fileName += ".";
    fileName += extension[snapshotFormat];
    info.setFile(fileName);
  } else if (info.suffix() != extension[snapshotFormat]) {
    // Extension is not appropriate. Propose a modification
    QString modifiedName = info.absolutePath() + '/' + info.baseName() + "." +
                           extension[snapshotFormat];
    QFileInfo modifInfo(modifiedName);
    int i = (QMessageBox::warning(
        widget, "Wrong extension",
        info.fileName() + " has a wrong extension.\nSave as " +
            modifInfo.fileName() + " instead ?",
        QMessageBox::Yes, QMessageBox::No, QMessageBox::Cancel));
    if (i == QMessageBox::Cancel)
      return false;

    if (i == QMessageBox::Yes) {
      fileName = modifiedName;
      info.setFile(fileName);
    }
  }

  return true;
}

class ImageInterface : public QDialog{
public:
  ImageInterface(QWidget *parent) : QDialog(parent) {}
};

// Pops-up an image settings dialog box and save to fileName.
// Returns false in case of problem.
bool QGLViewer::saveImageSnapshot(const QString &) {
    return false;
  }

  

/*! Saves a snapshot of the current image displayed by the widget.

 Options are set using snapshotFormat(), snapshotFileName() and
 snapshotQuality(). For non vectorial image formats, the image size is equal to
 the current viewer's dimensions (see width() and height()). See
 snapshotFormat() for details on supported formats.

 If \p automatic is \c false (or if snapshotFileName() is empty), a file dialog
 is opened to ask for the file name.

 When \p automatic is \c true, the file name is set to \c NAME-NUMBER, where \c
 NAME is snapshotFileName() and \c NUMBER is snapshotCounter(). The
 snapshotCounter() is automatically incremented after each snapshot saving. This
 is useful to create videos from your application: \code void Viewer::init()
 {
   resize(720, 576); // PAL DV format (use 720x480 for NTSC DV)
   connect(this, SIGNAL(drawFinished(bool)), SLOT(saveSnapshot(bool)));
 }
 \endcode
 Then call draw() in a loop (for instance using animate() and/or a camera()
 KeyFrameInterpolator replay) to create your image sequence.

 If you want to create a Quicktime VR panoramic sequence, simply use code like
 this: \code void Viewer::createQuicktime()
 {
   const int nbImages = 36;
   for (int i=0; i<nbImages; ++i)
         {
           camera()->setOrientation(2.0*M_PI/nbImages, 0.0); // Theta-Phi
 orientation showEntireScene(); update(); // calls draw(), which emits
 drawFinished(), which calls saveSnapshot()
         }
 }
 \endcode

 If snapshotCounter() is negative, no number is appended to snapshotFileName()
 and the snapshotCounter() is not incremented. This is useful to force the
 creation of a file, overwriting the previous one.

 When \p overwrite is set to \c false (default), a window asks for confirmation
 if the file already exists. In \p automatic mode, the snapshotCounter() is
 incremented (if positive) until a non-existing file name is found instead.
 Otherwise the file is overwritten without confirmation.

 The VRender library was written by Cyril Soler (Cyril dot Soler at imag dot
 fr). If the generated PS or EPS file is not properly displayed, remove the
 anti-aliasing option in your postscript viewer.

 \note In order to correctly grab the frame buffer, the QGLViewer window is
 raised in front of other windows by this method. */
void QGLViewer::saveSnapshot(bool automatic, bool overwrite) {
  // Ask for file name
  if (snapshotFileName().isEmpty() || !automatic) {
    QString fileName;
    QString selectedFormat = FDFormatString[snapshotFormat()];
    fileName = QFileDialog::getSaveFileName(
        this, "Choose a file name to save under", snapshotFileName(), formats,
        &selectedFormat,
        overwrite ? QFileDialog::DontConfirmOverwrite
                  : QFlags<QFileDialog::Option>(0));
    setSnapshotFormat(Qtformat[selectedFormat]);

    if (checkFileName(fileName, this, snapshotFormat()))
      setSnapshotFileName(fileName);
    else
      return;
  }

  QFileInfo fileInfo(snapshotFileName());

  if ((automatic) && (snapshotCounter() >= 0)) {
    // In automatic mode, names have a number appended
    const QString baseName = fileInfo.baseName();
    QString count;
    count.sprintf("%.04d", snapshotCounter_++);
    QString suffix;
    suffix = fileInfo.suffix();
    if (suffix.isEmpty())
      suffix = extension[snapshotFormat()];
    fileInfo.setFile(fileInfo.absolutePath() + '/' + baseName + '-' + count +
                     '.' + suffix);

    if (!overwrite)
      while (fileInfo.exists()) {
        count.sprintf("%.04d", snapshotCounter_++);
        fileInfo.setFile(fileInfo.absolutePath() + '/' + baseName + '-' +
                         count + '.' + fileInfo.suffix());
      }
  }

  bool saveOK;
      if (automatic) {
    QImage snapshot = frameBufferSnapshot();
    saveOK = snapshot.save(fileInfo.filePath(),
                           snapshotFormat().toLatin1().constData(),
                           snapshotQuality());
  } else
    saveOK = saveImageSnapshot(fileInfo.filePath());

  if (!saveOK)
    QMessageBox::warning(this, "Snapshot problem",
                         "Unable to save snapshot in\n" + fileInfo.filePath());
}

QImage QGLViewer::frameBufferSnapshot() {
  // Viewer must be on top of other windows.
  makeCurrent();
  raise();
  // Hack: Qt has problems if the frame buffer is grabbed after QFileDialog is
  // displayed. We grab the frame buffer before, even if it might be not
  // necessary (vectorial rendering). The problem could not be reproduced on a
  // simple example to submit a Qt bug. However, only grabs the backgroundImage
  // in the eponym example. May come from the driver.
  return QOpenGLWidget::grabFramebuffer();
}

/*! Same as saveSnapshot(), except that it uses \p fileName instead of
 snapshotFileName().

 If \p fileName is empty, opens a file dialog to select the name.

 Snapshot settings are set from snapshotFormat() and snapshotQuality().

 Asks for confirmation when \p fileName already exists and \p overwrite is \c
 false (default).

 \attention If \p fileName is a char* (as is "myFile.jpg"), it may be casted
 into a \c bool, and the other saveSnapshot() method may be used instead. Pass
 QString("myFile.jpg") as a parameter to prevent this. */
void QGLViewer::saveSnapshot(const QString &fileName, bool overwrite) {
  const QString previousName = snapshotFileName();
  const int previousCounter = snapshotCounter();
  setSnapshotFileName(fileName);
  setSnapshotCounter(-1);
  saveSnapshot(true, overwrite);
  setSnapshotFileName(previousName);
  setSnapshotCounter(previousCounter);
}

/*! Takes a snapshot of the current display and pastes it to the clipboard.

This action is activated by the KeyboardAction::SNAPSHOT_TO_CLIPBOARD enum,
binded to \c Ctrl+C by default.
*/
void QGLViewer::snapshotToClipboard() {
  QClipboard *cb = QApplication::clipboard();
  cb->setImage(frameBufferSnapshot());
}
