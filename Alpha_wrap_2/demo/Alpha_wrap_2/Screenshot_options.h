// Copyright (c) 2019-2020 X, The Moonshot Factory (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Author(s)     : Pierre Alliez pierre.alliez@inria.fr
//               : Michael Hemmer mhsaar@gmail.com
//               : Cedric Portaneri cportaneri@gmail.com
//
#ifndef CGAL_ALPHA_WRAP_2_SCREENSHOT_OPTIONS_H
#define CGAL_ALPHA_WRAP_2_SCREENSHOT_OPTIONS_H

#include <QFileDialog>

#include "ui_Screenshot_options.h"

/**
 * @class ScreenshotOptions
 * @brief Interface to tune the screenshot options
 */
class ScreenshotOptions
  : public QDialog, private Ui::ScreenshotOptions
{
  Q_OBJECT

  private:
  QString screenshot_dir_path = QDir::currentPath();

public:
  ScreenshotOptions(QWidget*,
                    QString folder,
                    QString filename,
                    int num)
  {
    setupUi(this);
    connect(pushButtonFolder, SIGNAL(clicked()), this,
            SLOT(on_pushButtonFolderTriggered()));
    connect(buttonBox, SIGNAL(accepted()), this,
            SLOT(accept()));
    connect(buttonBox, SIGNAL(rejected()), this,
            SLOT(close()));
    textEditFolder->setPlainText(folder);
    textEditFilename->setPlainText(filename);
    spinBoxNumbering->setValue(num);
  }

  QString get_current_screenshot_dir() {
    return screenshot_dir_path;
  }

  QString get_current_screenshot_filename() {
    return textEditFilename->toPlainText();
  }

  int get_current_screenshot_numbering_start() {
    return spinBoxNumbering->value();
  }

public Q_SLOTS:
  void on_pushButtonFolderTriggered()
  {
    screenshot_dir_path = QFileDialog::getExistingDirectory(0, ("Select screenshot output folder"), QDir::currentPath());
    QStringList all_dir = screenshot_dir_path.split(QDir::separator());
    textEditFolder->setPlainText(QDir::separator()+all_dir[1]+QDir::separator()+".."+QDir::separator()+all_dir[all_dir.size()-1]);
  }
};

#endif  // CGAL_ALPHA_WRAP_2_SCREENSHOT_OPTIONS_H
