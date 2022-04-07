// Copyright (c) 2012, 2020  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#include "OverlayDialog.h"
#include "ArrangementDemoWindow.h"
#include "ArrangementTypesUtils.h"
#include "QtMetaTypes.h"

#include <QListWidgetItem>
#include <QString>
#include <QItemSelectionModel>

#include "ui_OverlayDialog.h"

OverlayDialog::OverlayDialog(
  QWidget* parent, const std::vector<ArrangementInfo>& arr_infos_) :
    QDialog(parent),
    arr_infos(arr_infos_),
    ui(new Ui::OverlayDialog)
{
  using namespace demo_types;

  // An extra parenthesis around QColor to avoid the
  // https://en.wikipedia.org/wiki/Most_vexing_parse
  // on clang
  QBrush segColor( ( QColor( ::Qt::red ) ) );
  QBrush polyColor( ( QColor( ::Qt::darkGreen ) ) );
  QBrush conicColor( ( QColor( ::Qt::blue ) ) );
  this->ui->setupUi( this );

  for (std::size_t i = 0; i < arr_infos.size(); i++)
  {
    auto& arr = arr_infos[i];

    QListWidgetItem* item =
      new QListWidgetItem( this->ui->arrangementsListWidget );
    item->setText(arr.label);
    item->setData(ARRANGEMENT, QVariant::fromValue(i));

    static constexpr std::array<const char*, 5> icons = {
      ":/cgal/icons/green_icon.xpm", ":/cgal/icons/yellow_icon.xpm",
      ":/cgal/icons/blue_icon.xpm", ":/cgal/icons/red_icon.xpm",
      ":/cgal/icons/pink_icon.xpm"};

    QIcon icon;
    icon.addFile(
      QString::fromUtf8(icons[static_cast<int>(arr.ttype) % icons.size()]));

    item->setIcon(icon);
  }
}

std::vector<int> OverlayDialog::selectedArrangements() const
{
  std::vector<int> res;
  for (int i = 0; i < this->ui->overlayListWidget->count(); ++i)
  {
    QListWidgetItem* item = this->ui->overlayListWidget->item(i);
    size_t arr_idx = item->data(ARRANGEMENT).value<std::size_t>();
    res.push_back(this->arr_infos[arr_idx].id);
  }
  return res;
}

void OverlayDialog::on_pickPushButton_pressed( )
{
  int currentIndex = this->ui->arrangementsListWidget->currentRow( );
  if ( currentIndex == -1 )
    return;

  if ( !(this->ui->arrangementsListWidget->item( currentIndex )->flags( ) &
         Qt::ItemIsEnabled) )
  {
    return;
  }
  QListWidgetItem* takenItem =
    this->ui->arrangementsListWidget->takeItem( currentIndex );
  this->ui->overlayListWidget->addItem( takenItem );

  if ( this->ui->overlayListWidget->count( ) == 2 )
  {
    this->ui->pickPushButton->setEnabled( false );
  }
  else
  {
    this->ui->pickPushButton->setEnabled( true );
  }
  QItemSelectionModel* selectionModel =
    this->ui->arrangementsListWidget->selectionModel( );
  selectionModel->clearSelection( );
  this->restrictSelection( takenItem );
}

void OverlayDialog::on_unpickPushButton_pressed( )
{
  int currentIndex = this->ui->overlayListWidget->currentRow( );
  if ( currentIndex == -1 )
    return;

  QListWidgetItem* takenItem =
    this->ui->overlayListWidget->takeItem( currentIndex );
  this->ui->arrangementsListWidget->addItem( takenItem );

  if ( this->ui->overlayListWidget->count( ) == 2 )
  {
    this->ui->pickPushButton->setEnabled( false );
  }
  else
  {
    this->ui->pickPushButton->setEnabled( true );
  }
  if ( this->ui->overlayListWidget->count( ) == 0 )
    this->unrestrictSelection( );
}

void OverlayDialog::restrictSelection( QListWidgetItem* item )
{
  using namespace demo_types;

  auto& arr1 = arr_infos[item->data(ARRANGEMENT).value<std::size_t>()];

  for (int i = 0; i < this->ui->arrangementsListWidget->count(); ++i)
  {
    auto* otherItem = this->ui->arrangementsListWidget->item(i);
    auto& arr2 = arr_infos[otherItem->data(ARRANGEMENT).value<std::size_t>()];
    bool enabled = (arr1.ttype == arr2.ttype);

    Qt::ItemFlags flags = otherItem->flags();
    if (!enabled)
      flags &= ~(Qt::ItemIsEnabled);
    else
      flags |= Qt::ItemIsEnabled;
    otherItem->setFlags(flags);
  }
}

void OverlayDialog::unrestrictSelection( )
{
  for ( int i = 0; i < this->ui->arrangementsListWidget->count( ); ++i )
  {
    QListWidgetItem* otherItem = this->ui->arrangementsListWidget->item( i );
    Qt::ItemFlags flags = otherItem->flags( );
    flags |= Qt::ItemIsEnabled;
    otherItem->setFlags( flags );
  }
}
