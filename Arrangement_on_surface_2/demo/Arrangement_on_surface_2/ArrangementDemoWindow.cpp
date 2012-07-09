#include "ArrangementDemoWindow.h"
#include "NewTabDialog.h"
#include "OverlayDialog.h"
#include "ArrangementDemoPropertiesDialog.h"
#include "ArrangementDemoTab.h"
#include "Conic_reader.h"

#include <QActionGroup>
#include <QFileDialog>
#include <QMessageBox>

#include <CGAL/IO/Arr_with_history_iostream.h>
#include <CGAL/IO/Arr_text_formatter.h>
#include <CGAL/IO/Arr_with_history_text_formatter.h>


ArrangementDemoWindow::
ArrangementDemoWindow(QWidget* parent) :
    CGAL::Qt::DemosMainWindow( parent ),
    lastTabIndex( -1 ),
    ui( new Ui::ArrangementDemoWindow )
{
    this->setupUi( );

    // set up the demo window
    ArrangementDemoTabBase* demoTab = this->makeTab( SEGMENT_TRAITS ); 
    this->setupStatusBar( );
    this->setupOptionsMenu( );
    this->addAboutDemo( ":/help/about.html" );
    this->addAboutCGAL( );

    // set up callbacks
    QObject::connect( this->modeGroup, SIGNAL( triggered( QAction* ) ),
        this, SLOT( updateMode( QAction* ) ) );
    QObject::connect( this->envelopeGroup, SIGNAL( triggered( QAction* ) ),
        this, SLOT( updateEnvelope( QAction* ) ) );
    QObject::connect( this->snapGroup, SIGNAL( triggered( QAction* ) ),
        this, SLOT( updateSnapping( QAction* ) ) );
    QObject::connect( this->conicTypeGroup, SIGNAL( triggered( QAction* ) ),
        this, SLOT( updateConicType( QAction* ) ) );
}

ArrangementDemoWindow::
~ArrangementDemoWindow( )
{ }

ArrangementDemoTabBase*
ArrangementDemoWindow::
makeTab( TraitsType tt )
{
    static int tabLabelCounter = 1;
    QString tabLabel;

    ArrangementDemoTabBase* demoTab;
    Seg_arr* seg_arr;
    Pol_arr* pol_arr;
    Conic_arr* conic_arr;
    CGAL::Object arr;

    switch ( tt )
    {
    default:
    case SEGMENT_TRAITS:
        seg_arr = new Seg_arr;
        demoTab = new ArrangementDemoTab< Seg_arr >( seg_arr, 0 );
        arr = CGAL::make_object( seg_arr );
        tabLabel = QString( "%1 - Segment" ).arg( tabLabelCounter++ );
        break;
    case POLYLINE_TRAITS:
        pol_arr = new Pol_arr;
        demoTab = new ArrangementDemoTab< Pol_arr >( pol_arr, 0 );
        arr = CGAL::make_object( pol_arr );
        tabLabel = QString( "%1 - Polyline" ).arg( tabLabelCounter++ );
        break;
    case CONIC_TRAITS:
        conic_arr = new Conic_arr;
        demoTab = new ArrangementDemoTab< Conic_arr >( conic_arr, 0 );
        arr = CGAL::make_object( conic_arr );
        tabLabel = QString( "%1 - Conic" ).arg( tabLabelCounter++ );
        break;
    }

    this->arrangements.push_back( arr );
    this->tabs.push_back( demoTab );

    QGraphicsView* view = demoTab->getView( );
    this->addNavigation( view );
    this->ui->tabWidget->addTab( demoTab, tabLabel );
    this->lastTabIndex = this->ui->tabWidget->currentIndex( );
    this->ui->tabWidget->setCurrentWidget( demoTab );

    return demoTab;
}

ArrangementDemoTabBase*
ArrangementDemoWindow::
getTab( int tabIndex ) const
{
    if ( tabIndex < 0 || tabIndex > this->tabs.size( ) )
    {
        return NULL;
    }

    ArrangementDemoTabBase* tab = this->tabs[ tabIndex ];
    return tab;
}

ArrangementDemoTabBase*
ArrangementDemoWindow::
getCurrentTab( ) const
{
    int currentIndex = this->ui->tabWidget->currentIndex( );
    if ( currentIndex == -1 )
        return NULL;

    ArrangementDemoTabBase* res = this->tabs[ currentIndex ];
    return res;
}

std::vector< QString > 
ArrangementDemoWindow::
getTabLabels( ) const
{
    std::vector< QString > res;
    for ( int i = 0; i < this->ui->tabWidget->count( ); ++i )
    {
        res.push_back( this->ui->tabWidget->tabText( i ) );
    }
    return res;
}

std::vector< CGAL::Object > 
ArrangementDemoWindow::
getArrangements( ) const
{
    std::vector< CGAL::Object > res;
    for ( int i = 0; i < this->arrangements.size( ); ++i )
    {
        res.push_back( this->arrangements[ i ] );
    }
    return res;
}

void
ArrangementDemoWindow::
setupUi( )
{
    this->ui->setupUi( this );

    this->modeGroup = new QActionGroup( this );
    this->modeGroup->addAction( this->ui->actionDrag );
    this->modeGroup->addAction( this->ui->actionInsert );
    this->modeGroup->addAction( this->ui->actionDelete );
    this->modeGroup->addAction( this->ui->actionPointLocation );
    this->modeGroup->addAction( this->ui->actionRayShootingUp );
    this->modeGroup->addAction( this->ui->actionRayShootingDown );
    this->modeGroup->addAction( this->ui->actionMerge );
    this->modeGroup->addAction( this->ui->actionSplit );
    this->activeModes.push_back( this->ui->actionInsert );

    this->envelopeGroup = new QActionGroup( this );
    this->envelopeGroup->addAction( this->ui->actionLowerEnvelope );
    this->envelopeGroup->addAction( this->ui->actionUpperEnvelope );
    this->envelopeGroup->setExclusive( false );

    this->snapGroup = new QActionGroup( this );
    this->snapGroup->addAction( this->ui->actionSnapMode );
    this->snapGroup->addAction( this->ui->actionGridSnapMode );
    this->snapGroup->setExclusive( false );
    this->ui->actionGridSnapMode->setEnabled( false );

    this->conicTypeGroup = new QActionGroup( this );
    this->conicTypeGroup->addAction( this->ui->actionConicSegment );
    this->conicTypeGroup->addAction( this->ui->actionConicCircle );
    this->conicTypeGroup->addAction( this->ui->actionConicEllipse );
    this->conicTypeGroup->addAction( this->ui->actionConicThreePoint );
    this->conicTypeGroup->addAction( this->ui->actionConicFivePoint );
}

void
ArrangementDemoWindow::
updateMode( QAction* newMode )
{
    //QWidget* widget = this->ui->tabWidget->currentWidget( );
    //ArrangementDemoTabBase* demoTab = static_cast< ArrangementDemoTabBase* >( widget );
    const int TabIndex = this->ui->tabWidget->currentIndex( );
    if ( TabIndex == -1 ) return;
    ArrangementDemoTabBase* activeTab = this->tabs[ TabIndex ];
    QGraphicsScene* activeScene = activeTab->getScene( );
    QGraphicsView* activeView = activeTab->getView( );

    this->resetCallbackState( TabIndex );
    this->removeCallback( TabIndex );

    // update the active mode
    this->activeModes[ TabIndex ] = newMode;

    // hook up the new active mode
    if ( newMode == this->ui->actionInsert )
    {
        activeScene->installEventFilter( activeTab->getCurveInputCallback( ) );
    }
    else if ( newMode == this->ui->actionDrag )
    {
        activeView->setDragMode( QGraphicsView::ScrollHandDrag );
    }
    else if ( newMode == this->ui->actionDelete )
    {
        activeScene->installEventFilter( activeTab->getDeleteCurveCallback( ) );
    }
    else if ( newMode == this->ui->actionPointLocation )
    {
        activeScene->installEventFilter( activeTab->getPointLocationCallback( ) );
    }
    else if ( newMode == this->ui->actionRayShootingUp )
    {
        // -y is up for Qt, so we shoot down
        activeTab->getVerticalRayShootCallback( )->setShootingUp( false );
        activeScene->installEventFilter( activeTab->getVerticalRayShootCallback( ) );
    }
    else if ( newMode == this->ui->actionRayShootingDown )
    {
        // the bottom of the viewport for Qt is +y, so we shoot up
        activeTab->getVerticalRayShootCallback( )->setShootingUp( true );
        activeScene->installEventFilter( activeTab->getVerticalRayShootCallback( ) );
    }
    else if ( newMode == this->ui->actionMerge )
    {
        activeScene->installEventFilter( activeTab->getMergeEdgeCallback( ) );
    }
    else if ( newMode == this->ui->actionSplit )
    {
        activeScene->installEventFilter( activeTab->getSplitEdgeCallback( ) );
    }
}

void
ArrangementDemoWindow::
resetCallbackState( int tabIndex )
{
    if ( tabIndex == -1 )
    {
        return;
    }

    ArrangementDemoTabBase* activeTab = this->tabs[ tabIndex ];
    QAction* activeMode = this->activeModes[ tabIndex ];

    // unhook the old active mode
    if ( activeMode == this->ui->actionInsert )
    {
    
    }
    else if ( activeMode == this->ui->actionDrag )
    {
   
    }
    else if ( activeMode == this->ui->actionDelete )
    {
        activeTab->getDeleteCurveCallback( )->reset( );
    }
    else if ( activeMode == this->ui->actionPointLocation )
    {
        activeTab->getPointLocationCallback( )->reset( );
    }
    else if ( activeMode == this->ui->actionRayShootingUp )
    {
        activeTab->getVerticalRayShootCallback( )->reset( );
    }
    else if ( activeMode == this->ui->actionRayShootingDown )
    {
        activeTab->getVerticalRayShootCallback( )->reset( );
    }
    else if ( activeMode == this->ui->actionMerge )
    {
        activeTab->getMergeEdgeCallback( )->reset( );
    }
    else if ( activeMode == this->ui->actionSplit )
    {
        activeTab->getSplitEdgeCallback( )->reset( );
    }
}

void
ArrangementDemoWindow::
removeCallback( int tabIndex )
{
    if ( tabIndex == -1 )
    {
        return;
    }

    ArrangementDemoTabBase* activeTab = this->tabs[ tabIndex ];
    QGraphicsScene* activeScene = activeTab->getScene( );
    QGraphicsView* activeView = activeTab->getView( );
    QAction* activeMode = this->activeModes[ tabIndex ];

    // unhook the old active mode
    if ( activeMode == this->ui->actionInsert )
    {
        activeScene->removeEventFilter( activeTab->getCurveInputCallback( ) );
    }
    else if ( activeMode == this->ui->actionDrag )
    {
        activeView->setDragMode( QGraphicsView::NoDrag );
    }
    else if ( activeMode == this->ui->actionDelete )
    {
        activeScene->removeEventFilter( activeTab->getDeleteCurveCallback( ) );
    }
    else if ( activeMode == this->ui->actionPointLocation )
    {
        activeScene->removeEventFilter( activeTab->getPointLocationCallback( ) );
    }
    else if ( activeMode == this->ui->actionRayShootingUp )
    {
        activeScene->removeEventFilter( activeTab->getVerticalRayShootCallback( ) );
    }
    else if ( activeMode == this->ui->actionRayShootingDown )
    {
        activeScene->removeEventFilter( activeTab->getVerticalRayShootCallback( ) );
    }
    else if ( activeMode == this->ui->actionMerge )
    {
        activeScene->removeEventFilter( activeTab->getMergeEdgeCallback( ) );
    }
    else if ( activeMode == this->ui->actionSplit )
    {
        activeScene->removeEventFilter( activeTab->getSplitEdgeCallback( ) );
    }
}

void
ArrangementDemoWindow::
updateEnvelope( QAction* newMode )
{
    if ( this->ui->tabWidget->currentIndex( ) == -1 ) return;
    ArrangementDemoTabBase* activeTab = this->tabs[ this->ui->tabWidget->currentIndex( ) ];
    QGraphicsScene* activeScene = activeTab->getScene( );
    QGraphicsView* activeView = activeTab->getView( );

    bool show = newMode->isChecked( );
    if ( newMode == this->ui->actionLowerEnvelope )
    {
        activeTab->getEnvelopeCallback( )->showUpperEnvelope( show );
    }
    else if ( newMode == this->ui->actionUpperEnvelope )
    {
        activeTab->getEnvelopeCallback( )->showLowerEnvelope( show );
    }
}

void
ArrangementDemoWindow::
updateSnapping( QAction* newMode )
{
    ArrangementDemoTabBase* activeTab = this->tabs[ this->ui->tabWidget->currentIndex( ) ];
    QGraphicsScene* activeScene = activeTab->getScene( );
    ArrangementDemoGraphicsView* activeView = activeTab->getView( );

    bool enabled = newMode->isChecked( );
    if ( newMode == this->ui->actionSnapMode )
    {
        activeTab->getCurveInputCallback( )->setSnappingEnabled( enabled );
        activeTab->getSplitEdgeCallback( )->setSnappingEnabled( enabled );
        if ( ! enabled )
        {
            this->ui->actionGridSnapMode->setChecked( false );
            this->ui->actionGridSnapMode->setEnabled( false );
            activeTab->getCurveInputCallback( )->setSnapToGridEnabled( false );
            activeTab->getSplitEdgeCallback( )->setSnapToGridEnabled( false );
        }
        else
        {
            this->ui->actionGridSnapMode->setEnabled( true );
        }
    }
    else if ( newMode == this->ui->actionGridSnapMode )
    {
        activeTab->getCurveInputCallback( )->setSnapToGridEnabled( enabled );
        activeTab->getSplitEdgeCallback( )->setSnapToGridEnabled( enabled );
        activeView->setShowGrid( enabled );
    }
    activeScene->update( );
}

void
ArrangementDemoWindow::
updateConicType( QAction* newType )
{
    ArrangementDemoTabBase* activeTab = this->tabs[ this->ui->tabWidget->currentIndex( ) ];
    QGraphicsScene* activeScene = activeTab->getScene( );
    ArrangementDemoGraphicsView* activeView = activeTab->getView( );
    Conic_arr* conic_arr;
    bool isConicArr = CGAL::assign( conic_arr, this->arrangements[ this->ui->tabWidget->currentIndex( ) ] );
    if ( isConicArr )
    {
        std::cout << "do something conic arr related" << std::endl;
        typedef CGAL::Qt::GraphicsViewCurveInput< typename Conic_arr::Geometry_traits_2 > ConicCurveInputCallback;
        ConicCurveInputCallback* curveInputCallback = ( ConicCurveInputCallback* ) activeTab->getCurveInputCallback( );
        if ( newType == this->ui->actionConicSegment )
        {
            curveInputCallback->setConicType( ConicCurveInputCallback::CONIC_SEGMENT );
        }
        else if ( newType == this->ui->actionConicCircle )
        {
            curveInputCallback->setConicType( ConicCurveInputCallback::CONIC_CIRCLE );
        }
        else if ( newType == this->ui->actionConicEllipse )
        {
            curveInputCallback->setConicType( ConicCurveInputCallback::CONIC_ELLIPSE );
        }
        else if ( newType == this->ui->actionConicThreePoint )
        {
            curveInputCallback->setConicType( ConicCurveInputCallback::CONIC_THREE_POINT );
        }
        else if ( newType == this->ui->actionConicFivePoint )
        {
            curveInputCallback->setConicType( ConicCurveInputCallback::CONIC_FIVE_POINT );
        }
    }
    else
    {
        //std::cout << "do nothing" << std::endl;
    }
}

void
ArrangementDemoWindow::
on_actionSaveAs_triggered( )
{
    int index = this->ui->tabWidget->currentIndex( );
    if ( index == -1 )
        return;
    QString filename = 
        QFileDialog::getSaveFileName( this, tr( "Save file" ),
            "", "Arrangement (*.arr)" );
    if ( filename.isNull( ) )
        return;

    std::ofstream ofs( filename.toStdString( ).c_str( ) );
    CGAL::Object arr = this->arrangements[ index ];
    Seg_arr* seg;
    Pol_arr* pol;
    Conic_arr* conic;
    if ( CGAL::assign( seg, arr ) )
    {
        typedef CGAL::Arr_with_history_text_formatter< CGAL::Arr_text_formatter< Seg_arr > > ArrFormatter;
        ArrFormatter arrFormatter;
        CGAL::write( *seg, ofs, arrFormatter );
    }
    else if ( CGAL::assign( pol, arr ) )
    {
        typedef CGAL::Arr_with_history_text_formatter< CGAL::Arr_text_formatter< Pol_arr > > ArrFormatter;
        ArrFormatter arrFormatter;
        CGAL::write( *pol, ofs, arrFormatter );
    }
    else if ( CGAL::assign( conic, arr ) )
    {
#if 0
        typedef CGAL::Arr_with_history_text_formatter< CGAL::Arr_text_formatter< Conic_arr > > ArrFormatter;
        ArrFormatter arrFormatter;
        CGAL::write( *conic, ofs, arrFormatter );
#endif
        ofs << conic->number_of_curves( ) << std::endl;
        for ( typename Conic_arr::Curve_iterator it = conic->curves_begin( ); it != conic->curves_end( ); ++it )
        {
            if ( it->is_full_conic( ) )
            {
                ofs << "F ";
                ofs << it->r( ) << " ";
                ofs << it->s( ) << " ";
                ofs << it->t( ) << " ";
                ofs << it->u( ) << " ";
                ofs << it->v( ) << " ";
                ofs << it->w( ) << " ";
                ofs << std::endl;
            }
            else if ( it->orientation( ) == CGAL::COLLINEAR )
            {
                ofs << "S ";
                ofs << it->source( ) << " ";
                ofs << it->target( ) << " ";
                ofs << std::endl;
            }
            else
            {
                ofs << "A ";
                ofs << it->r( ) << " ";
                ofs << it->s( ) << " ";
                ofs << it->t( ) << " ";
                ofs << it->u( ) << " ";
                ofs << it->v( ) << " ";
                ofs << it->w( ) << " ";
                if ( it->orientation( ) == CGAL::COUNTERCLOCKWISE )
                    ofs << "1 ";
                else if ( it->orientation( ) == CGAL::CLOCKWISE )
                    ofs << "-1 ";
                else
                    ofs << "0 ";
                ofs << it->source( ) << " ";
                ofs << it->target( ) << " ";
                ofs << std::endl;
            }
        }
    }
    ofs.close( );
}

void
ArrangementDemoWindow::
on_actionOpen_triggered( )
{
    int index = this->ui->tabWidget->currentIndex( );
    if ( index == -1 )
    {
        QMessageBox::information( this, "Oops", "Create a new tab first" );
        return;
    }
    QString filename = 
        QFileDialog::getOpenFileName( this, tr( "Open file" ),
            "", "Arrangement file - *.arr (*.arr);;All files (*.*)" );
    if ( filename.isNull( ) )
        return;
    std::ifstream ifs( filename.toStdString( ).c_str( ) );
    CGAL::Object arr = this->arrangements[ index ];
    Seg_arr* seg;
    Pol_arr* pol;
    Conic_arr* conic;
    if ( CGAL::assign( seg, arr ) )
    {
        typedef CGAL::Arr_with_history_text_formatter< CGAL::Arr_text_formatter< Seg_arr > > ArrFormatter;
        typedef ArrangementDemoTab< Seg_arr > TabType;

        ArrFormatter arrFormatter;
        CGAL::read( *seg, ifs, arrFormatter );
        this->arrangements[ index ] = CGAL::make_object( seg );
        TabType* tab = static_cast< TabType* >( this->tabs[ index ] );
        tab->setArrangement( seg );
    }
    else if ( CGAL::assign( pol, arr ) )
    {
        typedef CGAL::Arr_with_history_text_formatter< CGAL::Arr_text_formatter< Pol_arr > > ArrFormatter;
        typedef ArrangementDemoTab< Pol_arr > TabType;

        ArrFormatter arrFormatter;
        CGAL::read( *pol, ifs, arrFormatter );
        this->arrangements[ index ] = CGAL::make_object( pol );
        TabType* tab = static_cast< TabType* >( this->tabs[ index ] );
        tab->setArrangement( pol );
    }
    else if ( CGAL::assign( conic, arr ) )
    {
#if 0
        typedef CGAL::Arr_with_history_text_formatter< CGAL::Arr_text_formatter< Conic_arr > > ArrFormatter;
        ArrFormatter arrFormatter;
        CGAL::read( *conic, ifs, arrFormatter );
        this->arrangements[ index ] = CGAL::make_object( conic );
        tab->setArrangement( conic );
#endif
        typedef ArrangementDemoTab< Conic_arr > TabType;
        Conic_reader< typename Conic_arr::Geometry_traits_2 > conicReader;
        std::vector< typename Conic_arr::Curve_2 > curve_list;
        CGAL::Bbox_2 bbox;
        conicReader.read_data( filename.toStdString( ).c_str( ),
            std::back_inserter( curve_list ), bbox );
        this->arrangements[ index ] = CGAL::make_object( conic );
        TabType* tab = static_cast< TabType* >( this->tabs[ index ] );

        CGAL::insert( *conic, curve_list.begin(), curve_list.end() );
        tab->setArrangement( conic );
        //QMessageBox::information( this, "Oops", "Reading conic arrangement not supported" );
    }
    ifs.close( );
}

void 
ArrangementDemoWindow::
on_actionQuit_triggered( )
{
    qApp->exit( ); 
}

void
ArrangementDemoWindow::
on_actionNewTab_triggered( )
{
    NewTabDialog* newTabDialog = new NewTabDialog;
    if ( newTabDialog->exec( ) == QDialog::Accepted )
    {
        int id = newTabDialog->checkedId( );
        if ( id == SEGMENT_TRAITS )
        {
            this->makeTab( SEGMENT_TRAITS );
        }
        else if ( id == POLYLINE_TRAITS )
        {
            this->makeTab( POLYLINE_TRAITS );
        }
        else if ( id == CONIC_TRAITS )
        {
            this->makeTab( CONIC_TRAITS );
        }
        else
        {
            std::cout << "Sorry, this trait is not yet supported" << std::endl;
        }
    }
    delete newTabDialog;
}

void
ArrangementDemoWindow::
on_tabWidget_currentChanged( )
{
    std::cout << "Tab changed" << std::endl;
    // disable the callback for the previously active tab
    this->resetCallbackState( this->lastTabIndex );
    this->removeCallback( this->lastTabIndex );
    this->lastTabIndex = this->ui->tabWidget->currentIndex( );

    this->updateMode( this->modeGroup->checkedAction( ) );

    CGAL::Object arr;
    if ( this->ui->tabWidget->currentIndex( ) != -1 )
        arr = this->arrangements[ this->ui->tabWidget->currentIndex( ) ];
    Seg_arr* seg;
    Pol_arr* pol;
    Conic_arr* conic;
    if ( CGAL::assign( conic, arr ) )
    {
        this->conicTypeGroup->setEnabled( true );
    }
    else
    {
        this->conicTypeGroup->setEnabled( false );
    }
}

void
ArrangementDemoWindow::
on_actionOverlay_triggered( )
{
    OverlayDialog* overlayDialog = new OverlayDialog( this );
    if ( overlayDialog->exec( ) == QDialog::Accepted )
    {
        std::vector< CGAL::Object > arrs = overlayDialog->selectedArrangements( );
        if ( arrs.size( ) == 2 )
        {
            Seg_arr* seg_arr;
            Pol_arr* pol_arr;
            Conic_arr* conic_arr;
            Seg_arr* seg_arr2;
            Pol_arr* pol_arr2;
            Conic_arr* conic_arr2;
            if ( CGAL::assign( seg_arr, arrs[ 0 ] ) && CGAL::assign( seg_arr2, arrs[ 1 ] ) )
            {
                this->makeOverlayTab( seg_arr, seg_arr2 );
            }
            if ( CGAL::assign( pol_arr, arrs[ 0 ] ) && CGAL::assign( pol_arr2, arrs[ 1 ] ) )
            {
                this->makeOverlayTab( pol_arr, pol_arr2 );
            }
            if ( CGAL::assign( conic_arr, arrs[ 0 ] ) && CGAL::assign( conic_arr2, arrs[ 1 ] ) )
            {
                this->makeOverlayTab( conic_arr, conic_arr2 );
            }


    #if 0
            if ( CGAL::assign( conic_arr, arrs[ 0 ] ) || CGAL::assign( conic_arr, arrs[ 1 ] ) )
            {
                this->makeOverlayTab< Conic_arr >( arrs[ 0 ], arrs[ 1 ] );
            }
            else if ( CGAL::assign( pol_arr, arrs[ 0 ] ) || CGAL::assign( pol_arr, arrs[ 1 ] ) )
            {
                this->makeOverlayTab< Pol_arr >( arrs[ 0 ], arrs[ 1 ] );
            }
            else
            {
                this->makeOverlayTab< Seg_arr >( arrs[ 0 ], arrs[ 1 ] );
            }
#endif
        }
    }
    delete overlayDialog;
}

void
ArrangementDemoWindow::
on_actionCloseTab_triggered( )
{
    int currentTabIndex = this->ui->tabWidget->currentIndex( );
    if ( ! this->ui->tabWidget->count( ) || currentTabIndex == -1 )
    {
        return;
    }

    // delete the tab
    this->ui->tabWidget->removeTab( currentTabIndex );
    this->tabs.erase( this->tabs.begin( ) + currentTabIndex );

    // delete the arrangement
    this->arrangements.erase( this->arrangements.begin( ) + currentTabIndex );
}

void
ArrangementDemoWindow::
on_actionPrintConicCurves_triggered( )
{
    int currentTabIndex = this->ui->tabWidget->currentIndex( );
    Conic_arr* arr;
    if ( currentTabIndex == -1 )
        return;
    CGAL::Object o = this->arrangements[ currentTabIndex ];
    if ( ! CGAL::assign( arr, o ) )
        return;
    typedef typename Conic_arr::Curve_iterator Curve_iterator;
    std::cout << arr->number_of_curves( ) << std::endl;
    for ( Curve_iterator it = arr->curves_begin( ); it != arr->curves_end( ); ++it )
    {
        std::cout << *it << std::endl;
        if ( it->is_full_conic( ) )
        {
            std::cout << "F ";
            std::cout << it->r( ) << " ";
            std::cout << it->s( ) << " ";
            std::cout << it->t( ) << " ";
            std::cout << it->u( ) << " ";
            std::cout << it->v( ) << " ";
            std::cout << it->w( ) << " ";
            std::cout << std::endl;
        }
        else if ( it->orientation( ) == CGAL::COLLINEAR )
        {
            std::cout << "S ";
            std::cout << it->source( ) << " ";
            std::cout << it->target( ) << " ";
            std::cout << std::endl;
        }
        else
        {
            std::cout << "A ";
            std::cout << it->r( ) << " ";
            std::cout << it->s( ) << " ";
            std::cout << it->t( ) << " ";
            std::cout << it->u( ) << " ";
            std::cout << it->v( ) << " ";
            std::cout << it->w( ) << " ";
            if ( it->orientation( ) == CGAL::COUNTERCLOCKWISE )
                std::cout << "1 ";
            else if ( it->orientation( ) == CGAL::CLOCKWISE )
                std::cout << "-1 ";
            else
                std::cout << "0 ";
            std::cout << it->source( ) << " ";
            std::cout << it->target( ) << " ";
            std::cout << std::endl;
        }
    }
}

void
ArrangementDemoWindow::
on_actionPreferences_triggered( )
{
    int currentTabIndex = this->ui->tabWidget->currentIndex( );
    if ( currentTabIndex == -1 )
        return;
    ArrangementDemoTabBase* currentTab = this->tabs[ currentTabIndex ];
    CGAL::Qt::ArrangementGraphicsItemBase* agi = currentTab->getArrangementGraphicsItem( );
    QPen vertexPen = agi->getVerticesPen( );
    QPen edgePen = agi->getEdgesPen( );
    QBrush vertexPenBrush = vertexPen.brush( );
    QColor vertexColor = vertexPenBrush.color( );
    QBrush edgePenBrush = edgePen.brush( );
    QColor edgeColor = edgePenBrush.color( );

    ArrangementDemoPropertiesDialog* dialog = new ArrangementDemoPropertiesDialog( this );
    if ( dialog->exec( ) == QDialog::Accepted )
    {

    }
}
