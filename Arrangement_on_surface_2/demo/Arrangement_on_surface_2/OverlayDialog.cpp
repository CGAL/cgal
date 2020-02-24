// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#include "OverlayDialog.h"
#include "ArrangementDemoWindow.h"
#include "ArrangementTypes.h"

#include <QListWidgetItem>
#include <QString>
#include <QItemSelectionModel>

#include "ui_OverlayDialog.h"

// TODO: Don't color the text, but set an icon for each arrangement type...
// TODO: Or maybe don't bother
OverlayDialog::OverlayDialog( ArrangementDemoWindow* parent,
                              Qt::WindowFlags f ) :
  QDialog( parent, f ),
  ui( new Ui::OverlayDialog )
{
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
    Seg_arr* seg;
    Pol_arr* pol;

#ifdef CGAL_USE_CORE
    Conic_arr* conic;
#endif

    Lin_arr* lin;
    Arc_arr* arc;
    // Alg_seg_arr* alg;
    if ( CGAL::assign( seg, arrangements[ i ] ) )
    {
      icon.addFile(QString::fromUtf8(":/icons/green_icon.xpm"), QSize(),
                   QIcon::Normal, QIcon::Off);
    }
    else if ( CGAL::assign( pol, arrangements[ i ] ) )
    {
      icon.addFile(QString::fromUtf8(":/icons/yellow_icon.xpm"), QSize(),
                   QIcon::Normal, QIcon::Off);
    }

#ifdef CGAL_USE_CORE
    else if ( CGAL::assign( conic, arrangements[ i ] ) )
    {
      icon.addFile(QString::fromUtf8(":/icons/red_icon.xpm"), QSize(),
                   QIcon::Normal, QIcon::Off);
    }
#endif

    else if ( CGAL::assign( lin, arrangements[ i ] ) )
    {
      icon.addFile(QString::fromUtf8(":/icons/blue_icon.xpm"), QSize(),
                   QIcon::Normal, QIcon::Off);
    }
    else if ( CGAL::assign( arc, arrangements[ i ] ) )
    {
      icon.addFile(QString::fromUtf8(":/icons/green_icon.xpm"), QSize(),
                   QIcon::Normal, QIcon::Off);
    }
    // else if ( CGAL::assign( alg, arrangements[ i ] ) )
    // {
    //   icon.addFile(QString::fromUtf8(":/icons/yellow_icon.xpm"), QSize(),
    //                QIcon::Normal, QIcon::Off);
    // }
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
  CGAL::Object o = item->data( ARRANGEMENT ).value< CGAL::Object >( );
  Seg_arr* seg;
  Pol_arr* pol;

#ifdef CGAL_USE_CORE
  Conic_arr* conic;
#endif

  Lin_arr* lin;
  Arc_arr* arc;
  // Alg_seg_arr* alg;
  if ( CGAL::assign( seg, o ) )
  {
    for ( int i = 0; i < this->ui->arrangementsListWidget->count( ); ++i )
    {
      QListWidgetItem* otherItem = this->ui->arrangementsListWidget->item( i );
      CGAL::Object o2 = otherItem->data( ARRANGEMENT ).value< CGAL::Object >( );
      bool enabled = CGAL::assign( seg, o2 );
      Qt::ItemFlags flags = otherItem->flags( );
      if ( ! enabled )
      {
        flags &= ~( Qt::ItemIsEnabled );
      }
      else
      {
        flags |= Qt::ItemIsEnabled;
      }
      otherItem->setFlags( flags );
    }
  }
  else if ( CGAL::assign( pol, o ) )
  {
    for ( int i = 0; i < this->ui->arrangementsListWidget->count( ); ++i )
    {
      QListWidgetItem* otherItem = this->ui->arrangementsListWidget->item( i );
      CGAL::Object o2 = otherItem->data( ARRANGEMENT ).value< CGAL::Object >( );
      bool enabled = CGAL::assign( pol, o2 );
      Qt::ItemFlags flags = otherItem->flags( );
      if ( ! enabled )
      {
        flags &= ~( Qt::ItemIsEnabled );
      }
      else
      {
        flags |= Qt::ItemIsEnabled;
      }
      otherItem->setFlags( flags );
    }
  }

#ifdef CGAL_USE_CORE
  else if ( CGAL::assign( conic, o ) )
  {
    for ( int i = 0; i < this->ui->arrangementsListWidget->count( ); ++i )
    {
      QListWidgetItem* otherItem = this->ui->arrangementsListWidget->item( i );
      CGAL::Object o2 = otherItem->data( ARRANGEMENT ).value< CGAL::Object >( );
      bool enabled = CGAL::assign( conic, o2 );
      Qt::ItemFlags flags = otherItem->flags( );
      if ( ! enabled )
      {
        flags &= ~( Qt::ItemIsEnabled );
      }
      else
      {
        flags |= Qt::ItemIsEnabled;
      }
      otherItem->setFlags( flags );
    }
  }
#endif

  else if ( CGAL::assign( lin, o ) )
  {
    for ( int i = 0; i < this->ui->arrangementsListWidget->count( ); ++i )
    {
      QListWidgetItem* otherItem = this->ui->arrangementsListWidget->item( i );
      CGAL::Object o2 = otherItem->data( ARRANGEMENT ).value< CGAL::Object >( );
      bool enabled = CGAL::assign( lin, o2 );
      Qt::ItemFlags flags = otherItem->flags( );
      if ( ! enabled )
      {
        flags &= ~( Qt::ItemIsEnabled );
      }
      else
      {
        flags |= Qt::ItemIsEnabled;
      }
      otherItem->setFlags( flags );
    }
  }
  else if ( CGAL::assign( arc, o ) )
  {
    for ( int i = 0; i < this->ui->arrangementsListWidget->count( ); ++i )
    {
      QListWidgetItem* otherItem = this->ui->arrangementsListWidget->item( i );
      CGAL::Object o2 = otherItem->data( ARRANGEMENT ).value< CGAL::Object >( );
      bool enabled = CGAL::assign( arc, o2 );
      Qt::ItemFlags flags = otherItem->flags( );
      if ( ! enabled )
      {
        flags &= ~( Qt::ItemIsEnabled );
      }
      else
      {
        flags |= Qt::ItemIsEnabled;
      }
      otherItem->setFlags( flags );
    }
  }
  // else if ( CGAL::assign( alg, o ) )
  // {
  //   for ( int i = 0; i < this->ui->arrangementsListWidget->count( ); ++i )
  //   {
  //     QListWidgetItem* otherItem = this->ui->arrangementsListWidget->item( i );
  //     CGAL::Object o2 = otherItem->data( ARRANGEMENT ).value< CGAL::Object >( );
  //     bool enabled = CGAL::assign( alg, o2 );
  //     Qt::ItemFlags flags = otherItem->flags( );
  //     if ( ! enabled )
  //     {
  //       flags &= ~( Qt::ItemIsEnabled );
  //     }
  //     else
  //     {
  //       flags |= Qt::ItemIsEnabled;
  //     }
  //     otherItem->setFlags( flags );
  //   }
  // }
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
