// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s) : Saar Katz <kats.saar@gmail.com>
//             Efi Fogel <efifogel@gmain.com>

#ifndef CGAL_CHOISEDIALOGPOLYGON_H
#define CGAL_CHOISEDIALOGPOLYGON_H

#include <QDialog>

#include "PolygonTableModel.h"

namespace Ui {
  class ChoiseDialogPolygon;
}

class ChoiseDialogPolygon : public QDialog
{
  Q_OBJECT

public:
  explicit ChoiseDialogPolygon(PolygonTableModel* options, QWidget* parent = 0);
  ~ChoiseDialogPolygon();

  PolygonWithHoles* getChoise();

private:
  Ui::ChoiseDialogPolygon* ui;
  PolygonTableModel* m_options;
  PolygonWithHoles* m_choise;

private slots:
  void on_buttonBox_accepted();
};

#endif // CGAL_CHOISEDIALOGPOLYGON_H
