#include "ArrangementDemoWindow.h"
#include <QActionGroup>
#include "NewTabDialog.h"

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
}

ArrangementDemoWindow::
~ArrangementDemoWindow( )
{ }

ArrangementDemoTabBase*
ArrangementDemoWindow::
makeTab( TraitsType tt )
{
    static int tabLabelCounter = 1;
    QString tabLabel = QString( "Tab %1" ).arg( tabLabelCounter++ );

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
        break;
    case POLYLINE_TRAITS:
        pol_arr = new Pol_arr;
        demoTab = new ArrangementDemoTab< Pol_arr >( pol_arr, 0 );
        arr = CGAL::make_object( pol_arr );
        break;
    case CONIC_TRAITS:
        conic_arr = new Conic_arr;
        demoTab = new ArrangementDemoTab< Conic_arr >( conic_arr, 0 );
        arr = CGAL::make_object( conic_arr );
        break;
    }

    this->arrangements.push_back( arr );
    this->tabs.push_back( demoTab );

    QGraphicsView* view = demoTab->getView( );
    this->addNavigation( view );
    this->ui->tabWidget->addTab( demoTab, tabLabel );
    this->lastTabIndex = this->ui->tabWidget->currentIndex( );

    return demoTab;
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
}

void
ArrangementDemoWindow::
updateMode( QAction* newMode )
{
    //QWidget* widget = this->ui->tabWidget->currentWidget( );
    //ArrangementDemoTabBase* demoTab = static_cast< ArrangementDemoTabBase* >( widget );
    const int TabIndex = this->ui->tabWidget->currentIndex( );
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
}
