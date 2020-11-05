// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#include "ArrangementDemoPropertiesDialog.h"
#include "ui_ArrangementDemoPropertiesDialog.h"
#include "ArrangementDemoWindow.h"
#include "PropertyValueDelegate.h"
#include "ArrangementDemoGraphicsView.h"
#include "VerticalRayShootCallback.h"
#include "ArrangementGraphicsItem.h"
#include "ArrangementDemoTab.h"
#include "EnvelopeCallback.h"
#include "GridGraphicsItem.h"


ArrangementDemoPropertiesDialog::
ArrangementDemoPropertiesDialog( ArrangementDemoWindow* parent_ ) :
  QDialog( parent_ ),
  parent( parent_ ),
  ui( new Ui::ArrangementDemoPropertiesDialog )
{
  this->setupUi( );
}

QVariant ArrangementDemoPropertiesDialog::property( int index )
{
  // return invalid variant if something goes wrong
  QVariant res;
  if ( index < 0 || index >= this->ui->tableWidget->rowCount( ) )
  {
    return res;
  }

  QTableWidgetItem* item = this->ui->tableWidget->item( index, 0 );
  if ( item == 0 )
  {
    return res;
  }

  // return user data, if it is set
  QVariant myData = item->data( Qt::UserRole );
  res = item->data( Qt::DisplayRole );

  return res;
}

//! function to setup the user interface that contains the buttons to draw and manipulate like segments and cones
/*!
*/
void ArrangementDemoPropertiesDialog::setupUi( )
{
  this->ui->setupUi( this );
  PropertyValueDelegate* myDelegate = new PropertyValueDelegate;
  this->ui->tableWidget->setItemDelegate( myDelegate );

  // populate the table widget with items
  QTableWidgetItem* edgeColorItem = new QTableWidgetItem;
  QTableWidgetItem* vertexColorItem = new QTableWidgetItem;
  QTableWidgetItem* edgeWidthItem = new QTableWidgetItem;
  QTableWidgetItem* vertexRadiusItem = new QTableWidgetItem;
  QTableWidgetItem* envelopeEdgeColorItem = new QTableWidgetItem;
  QTableWidgetItem* envelopeVertexColorItem = new QTableWidgetItem;
  QTableWidgetItem* envelopeEdgeWidthItem = new QTableWidgetItem;
  QTableWidgetItem* envelopeVertexRadiusItem = new QTableWidgetItem;
  QTableWidgetItem* verticalRayEdgeColorItem = new QTableWidgetItem;
  QTableWidgetItem* verticalRayEdgeWidthItem = new QTableWidgetItem;
  QTableWidgetItem* gridColorItem = new QTableWidgetItem;

  this->ui->tableWidget->setItem( int(EDGE_COLOR_KEY), 0, edgeColorItem );
  this->ui->tableWidget->setItem( int(VERTEX_COLOR_KEY), 0, vertexColorItem );
  this->ui->tableWidget->setItem( int(EDGE_WIDTH_KEY), 0, edgeWidthItem );
  this->ui->tableWidget->setItem( int(VERTEX_RADIUS_KEY), 0, vertexRadiusItem );
  this->ui->tableWidget->setItem( int(ENVELOPE_EDGE_COLOR_KEY), 0,
                  envelopeEdgeColorItem );
  this->ui->tableWidget->setItem( int(ENVELOPE_VERTEX_COLOR_KEY), 0,
                  envelopeVertexColorItem );
  this->ui->tableWidget->setItem( int(ENVELOPE_EDGE_WIDTH_KEY), 0,
                  envelopeEdgeWidthItem );
  this->ui->tableWidget->setItem( int(ENVELOPE_VERTEX_RADIUS_KEY), 0,
                  envelopeVertexRadiusItem );
  this->ui->tableWidget->setItem( int(VERTICAL_RAY_EDGE_COLOR_KEY), 0,
                  verticalRayEdgeColorItem );
  this->ui->tableWidget->setItem( int(VERTICAL_RAY_EDGE_WIDTH_KEY), 0,
                  verticalRayEdgeWidthItem );
  this->ui->tableWidget->setItem( int(GRID_COLOR_KEY), 0, gridColorItem );

  // fill in the items with data
  this->updateUi( );
}

//! Adding a new tab whenever the user creates a new file and updating the menu options accordingly
/*!
*/
void ArrangementDemoPropertiesDialog::updateUi( )
{
  if ( this->parent == NULL )
  {
    return;
  }
  ArrangementDemoTab* currentTab = this->parent->getCurrentTab();
  if ( currentTab == NULL )
  {
    return;
  }
  CGAL::Qt::ArrangementGraphicsItemBase* agi = currentTab->getArrangementGraphicsItem( );
  if ( agi == NULL )
  {
    return;
  }

  EnvelopeCallbackBase* envelopeCallback = currentTab->getEnvelopeCallback( );
  VerticalRayShootCallbackBase* verticalRayShootCallback =
  currentTab->getVerticalRayShootCallback( );
  GridGraphicsItem* gridGraphicsItem = currentTab->getGridGraphicsItem();

  QPen vertexPen = agi->getVerticesPen( );
  QPen edgePen = agi->getEdgesPen( );
  QBrush vertexPenBrush = vertexPen.brush( );
  QBrush edgePenBrush = edgePen.brush( );
  QColor vertexColor = vertexPenBrush.color( );
  QColor edgeColor = edgePenBrush.color( );
  unsigned int edgeWidth = edgePen.width( );
  unsigned int vertexRadius = vertexPen.width( );
  QColor gridColor = gridGraphicsItem->getAxesColor( );
  unsigned int envelopeEdgeWidth = envelopeCallback->getEnvelopeEdgeWidth( );
  unsigned int envelopeVertexRadius =
  envelopeCallback->getEnvelopeVertexRadius( );
  QColor envelopeEdgeColor = envelopeCallback->getEnvelopeEdgeColor( );
  QColor envelopeVertexColor = envelopeCallback->getEnvelopeVertexColor( );
  unsigned int verticalRayEdgeWidth = verticalRayShootCallback->edgeWidth( );
  QColor verticalRayEdgeColor = verticalRayShootCallback->edgeColor( );

  QTableWidgetItem* edgeColorItem =
  this->ui->tableWidget->item( int(EDGE_COLOR_KEY), 0 );
  QTableWidgetItem* edgeWidthItem =
  this->ui->tableWidget->item( int(EDGE_WIDTH_KEY), 0 );
  QTableWidgetItem* vertexColorItem =
  this->ui->tableWidget->item( int(VERTEX_COLOR_KEY), 0 );
  QTableWidgetItem* vertexRadiusItem =
  this->ui->tableWidget->item( int(VERTEX_RADIUS_KEY), 0 );
  QTableWidgetItem* gridColorItem =
  this->ui->tableWidget->item( int( GRID_COLOR_KEY ), 0 );
  QTableWidgetItem* envelopeEdgeColorItem =
  this->ui->tableWidget->item( int(ENVELOPE_EDGE_COLOR_KEY), 0 );
  QTableWidgetItem* envelopeEdgeWidthItem =
  this->ui->tableWidget->item( int(ENVELOPE_EDGE_WIDTH_KEY), 0 );
  QTableWidgetItem* envelopeVertexColorItem =
  this->ui->tableWidget->item( int(ENVELOPE_VERTEX_COLOR_KEY), 0 );
  QTableWidgetItem* envelopeVertexRadiusItem =
  this->ui->tableWidget->item( int(ENVELOPE_VERTEX_RADIUS_KEY), 0 );
  QTableWidgetItem* verticalRayEdgeColorItem =
  this->ui->tableWidget->item( int(VERTICAL_RAY_EDGE_COLOR_KEY), 0 );
  QTableWidgetItem* verticalRayEdgeWidthItem =
  this->ui->tableWidget->item( int(VERTICAL_RAY_EDGE_WIDTH_KEY), 0 );

  // arrangement properties
  edgeColorItem->setData( Qt::DisplayRole, edgeColor );
  edgeColorItem->setData( Qt::DecorationRole, edgeColor );
  edgeColorItem->setData( Qt::UserRole, QVariant::fromValue( edgeColor ) );

  edgeWidthItem->setData( Qt::DisplayRole, edgeWidth );

  vertexColorItem->setData( Qt::DisplayRole, vertexColor );
  vertexColorItem->setData( Qt::DecorationRole, vertexColor );
  vertexColorItem->setData( Qt::UserRole, QVariant::fromValue( vertexColor ) );

  vertexRadiusItem->setData( Qt::DisplayRole, vertexRadius );

  // envelope properties
  envelopeEdgeColorItem->setData( Qt::DisplayRole, envelopeEdgeColor );
  envelopeEdgeColorItem->setData( Qt::DecorationRole, envelopeEdgeColor );
  envelopeEdgeColorItem->setData( Qt::UserRole,
                  QVariant::fromValue( envelopeEdgeColor ) );

  envelopeEdgeWidthItem->setData( Qt::DisplayRole, envelopeEdgeWidth );

  envelopeVertexColorItem->setData( Qt::DisplayRole, envelopeVertexColor );
  envelopeVertexColorItem->setData( Qt::DecorationRole, envelopeVertexColor );
  envelopeVertexColorItem->setData( Qt::UserRole,
                  QVariant::fromValue(envelopeVertexColor));

  envelopeVertexRadiusItem->setData( Qt::DisplayRole, envelopeVertexRadius );

  // vertical ray properties
  verticalRayEdgeColorItem->setData( Qt::DisplayRole, verticalRayEdgeColor );
  verticalRayEdgeColorItem->setData( Qt::DecorationRole, verticalRayEdgeColor );
  verticalRayEdgeColorItem->setData( Qt::UserRole,
                   QVariant::fromValue(verticalRayEdgeColor));

  verticalRayEdgeWidthItem->setData( Qt::DisplayRole, verticalRayEdgeWidth );

  // grid properties
  gridColorItem->setData( Qt::DisplayRole, gridColor );
  gridColorItem->setData( Qt::DecorationRole, gridColor );
  gridColorItem->setData( Qt::UserRole, QVariant::fromValue( gridColor ) );
}
