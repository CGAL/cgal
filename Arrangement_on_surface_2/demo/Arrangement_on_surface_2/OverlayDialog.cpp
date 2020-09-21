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

OverlayDialog::OverlayDialog( ArrangementDemoWindow* parent ) :
  QDialog( parent ),
  ui( new Ui::OverlayDialog )
{
  using namespace demo_types;

  // An extra parenthesis around QColor to avoid the
  // https://en.wikipedia.org/wiki/Most_vexing_parse
  // on clang
  QBrush segColor( ( QColor( ::Qt::red ) ) );
  QBrush polyColor( ( QColor( ::Qt::darkGreen ) ) );
  QBrush conicColor( ( QColor( ::Qt::blue ) ) );
  this->ui->setupUi( this );

  std::vector< QString > labels = parent->getTabLabels( );
  std::vector< CGAL::Object > arrangements = parent->getArrangements( );

  for ( unsigned int i = 0; i < labels.size( ); ++i )
  {
    QListWidgetItem* item =
      new QListWidgetItem( this->ui->arrangementsListWidget );
    item->setText( labels[ i ] );
    item->setData( ARRANGEMENT, QVariant::fromValue( arrangements[ i ] ) );
    QIcon icon;

    TraitsType traitsType = TraitsType::NONE;
    forEachArrangementType([&](auto type_holder) {
      using Arrangement = typename decltype(type_holder)::type;
      Arrangement* arr;
      if (CGAL::assign(arr, arrangements[i]))
        traitsType = enumFromArrType<Arrangement>();
    });

    if (traitsType == TraitsType::NONE)
      CGAL_error();

    static constexpr std::array<const char*, 5> icons = {
      ":/cgal/icons/green_icon.xpm", ":/cgal/icons/yellow_icon.xpm",
      ":/cgal/icons/blue_icon.xpm", ":/cgal/icons/red_icon.xpm",
      ":/cgal/icons/pink_icon.xpm"};

    icon.addFile(
      QString::fromUtf8(icons[static_cast<int>(traitsType) % icons.size()]));

    item->setIcon( icon );
  }
}

std::vector< CGAL::Object >
OverlayDialog::selectedArrangements( ) const
{
  std::vector< CGAL::Object > res;
  for ( int i = 0; i < this->ui->overlayListWidget->count( ); ++i )
  {
    QListWidgetItem* item = this->ui->overlayListWidget->item( i );
    QVariant data = item->data( ARRANGEMENT );
    CGAL::Object arr = data.value< CGAL::Object >( );
    res.push_back( arr );
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

  CGAL::Object o = item->data( ARRANGEMENT ).value< CGAL::Object >( );

  forEachArrangementType([&](auto type_holder) {
    using Arrangement = typename decltype(type_holder)::type;
    Arrangement* arr;
    if (CGAL::assign(arr, o))
    {
      for (int i = 0; i < this->ui->arrangementsListWidget->count(); ++i)
      {
        auto* otherItem = this->ui->arrangementsListWidget->item(i);
        auto o2 = otherItem->data(ARRANGEMENT).value<CGAL::Object>();
        bool enabled = CGAL::assign(arr, o2);
        Qt::ItemFlags flags = otherItem->flags();
        if (!enabled) { flags &= ~(Qt::ItemIsEnabled); }
        else { flags |= Qt::ItemIsEnabled; }

        otherItem->setFlags(flags);
      }
    }
  });
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
