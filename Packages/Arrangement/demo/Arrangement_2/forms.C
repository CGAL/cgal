#include "forms.h"
#include "demo_tab.h"


const char* green_icon[]={
"16 16 2 1",
"g c green",
". c None",
"................",
"................",
"..gggggggggggg..",
"..gggggggggggg..",
"..gggggggggggg..",
"..ggg......ggg..",
"..ggg......ggg..",
"..ggg......ggg..",
"..ggg......ggg..",
"..ggg......ggg..",
"..ggg......ggg..",
"..gggggggggggg..",
"..gggggggggggg..",
"..gggggggggggg..",
"................",
"................"};



const char* white_icon[]={
"16 16 2 1",
"g c green",
". c None",
"................",
"................",
"................",
"................",
"................",
"................",
"................",
"................",
"................",
"................",
"................",
"................",
"................",
"................",
"................",
"................"};

////////////////////////////////////////////////////////////////////////
OptionsForm::OptionsForm(  QTabWidget * bar, QWidget* parent ,int number_of_tabs , const char* name, bool modal, WFlags f  ): 
	QDialog( parent, name, modal, f ),
	myBar(bar)
{
    setCaption( "Union -- Options" );
    resize( 320, 290 );

    optionsFormLayout = new QVBoxLayout( this, 11, 6 );

    arrLayout1 = new QHBoxLayout( 0, 0, 6 );

    textLabel1 = new QLabel( "Arr 1", this );
    arrLayout1->addWidget( textLabel1 );

    arrComboBox1 = new QComboBox( FALSE, this );

	for (int i=0; i < number_of_tabs; i++)
		arrComboBox1->insertItem( myBar->label(i) );
		
	arrLayout1->addWidget( arrComboBox1 );
    optionsFormLayout->addLayout( arrLayout1 );

	arrLayout2 = new QHBoxLayout( 0, 0, 6 );

    textLabel2 = new QLabel( "Arr 2", this );
    arrLayout2->addWidget( textLabel2 );

    arrComboBox2 = new QComboBox( FALSE, this );
    
	for (int i=0; i < number_of_tabs; i++)
		arrComboBox2->insertItem( myBar->label(i) );
	
	arrLayout2->addWidget( arrComboBox2 );
    optionsFormLayout->addLayout( arrLayout2 );

	buttonsLayout = new QHBoxLayout( 0, 0, 6 );

    okPushButton = new QPushButton( "OK", this );
    okPushButton->setDefault( TRUE );
    buttonsLayout->addWidget( okPushButton );

    cancelPushButton = new QPushButton( "Cancel", this );
    buttonsLayout->addWidget( cancelPushButton );
    optionsFormLayout->addLayout( buttonsLayout );

    connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept() ) );
    connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );

    textLabel1->setBuddy( arrComboBox1 );
	textLabel2->setBuddy( arrComboBox2 );
   
}
////////////////////////////////////////////////////////////////////////
PropertiesForm::PropertiesForm(  QTabWidget * bar, QWidget* parent ,int number_of_tabs , const char* name, bool modal, WFlags f  ): 
	QDialog( parent, name, modal, f ),
	myBar(bar)
{
    setCaption( "Properties -- Options" );
    resize( 320, 290 );

    optionsFormLayout = new QVBoxLayout( this, 11, 6 );

    arrLayout1 = new QHBoxLayout( 0, 0, 6 );
    textLabel1 = new QLabel( "Width", this );
    arrLayout1->addWidget( textLabel1 );
	box1 = new QSpinBox( 300, 1000, 50, this, "box1" );
	box1->setValue(700);
	arrLayout1->addWidget( box1 );
    optionsFormLayout->addLayout( arrLayout1 );

	arrLayout2 = new QHBoxLayout( 0, 0, 6 );
    textLabel2 = new QLabel( "Hight", this );
    arrLayout2->addWidget( textLabel2 );
	box2 = new QSpinBox( 300, 1000, 50, this, "box2" );
	box2->setValue(700);
	arrLayout2->addWidget( box2 );
	optionsFormLayout->addLayout( arrLayout2 );

	arrLayout3 = new QHBoxLayout( 0, 0, 6 );
    textLabel3 = new QLabel( "Line Width", this );
    arrLayout3->addWidget( textLabel3 );
	box3 = new QSpinBox( 1, 5, 1, this, "box3" );
	box3->setValue(2);
	arrLayout3->addWidget( box3 );
	optionsFormLayout->addLayout( arrLayout3 );

	arrLayout4 = new QHBoxLayout( 0, 0, 6 );
    textLabel4 = new QLabel( "Scaling Factor", this );
    arrLayout4->addWidget( textLabel4 );
	box4 = new MySpinBox( 10, 100, 1, this, "box4" );
	box4->setValue(20);
	arrLayout4->addWidget( box4 );
	optionsFormLayout->addLayout( arrLayout4 );
    
	buttonsLayout = new QHBoxLayout( 0, 0, 6 );
    okPushButton = new QPushButton( "OK", this );
    okPushButton->setDefault( TRUE );
    buttonsLayout->addWidget( okPushButton );
    cancelPushButton = new QPushButton( "Cancel", this );
    buttonsLayout->addWidget( cancelPushButton );
    optionsFormLayout->addLayout( buttonsLayout );
    connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept() ) );
    connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );

    textLabel1->setBuddy( box1 );
	textLabel2->setBuddy( box2 );
	textLabel3->setBuddy( box3 );
	textLabel4->setBuddy( box4 );
   
}
////////////////////////////////////////////////////////////////////////
OverLay::OverLay(  QTabWidget * bar, QWidget* parent ,int number_of_tabs , const char* name, bool modal, WFlags f  ): 
	QDialog( parent, name, modal, f ),
	myBar(bar)
{
    setCaption( "Union -- Options" );
    resize( 320, 290 );

    optionsFormLayout = new QVBoxLayout( this, 11, 6 );

    arrLayout1 = new QHBoxLayout( 0, 0, 6 );

    textLabel1 = new QLabel( "Arr 1", this );
    arrLayout1->addWidget( textLabel1 );

    arrComboBox1 = new QComboBox( FALSE, this );

	for (int i=0; i < number_of_tabs; i++)
		arrComboBox1->insertItem( myBar->label(i) );
		
	arrLayout1->addWidget( arrComboBox1 );
    optionsFormLayout->addLayout( arrLayout1 );

	arrLayout2 = new QHBoxLayout( 0, 0, 6 );

    textLabel2 = new QLabel( "Arr 2", this );
    arrLayout2->addWidget( textLabel2 );

    arrComboBox2 = new QComboBox( FALSE, this );
    
	for (int i=0; i < number_of_tabs; i++)
		arrComboBox2->insertItem( myBar->label(i) );
	
	arrLayout2->addWidget( arrComboBox2 );
    optionsFormLayout->addLayout( arrLayout2 );

	arrLayout3 = new QHBoxLayout( 0, 0, 6 );

    textLabel3 = new QLabel( "Arr 3", this );
    arrLayout3->addWidget( textLabel3 );

    arrComboBox3 = new QComboBox( FALSE, this );

	for (int i=0; i < number_of_tabs; i++)
		arrComboBox3->insertItem( myBar->label(i) );
		
	arrLayout3->addWidget( arrComboBox3 );
    optionsFormLayout->addLayout( arrLayout3 );

	buttonsLayout = new QHBoxLayout( 0, 0, 6 );

    okPushButton = new QPushButton( "OK", this );
    okPushButton->setDefault( TRUE );
    buttonsLayout->addWidget( okPushButton );

    cancelPushButton = new QPushButton( "Cancel", this );
    buttonsLayout->addWidget( cancelPushButton );
    optionsFormLayout->addLayout( buttonsLayout );

    connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept() ) );
    connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );

    textLabel1->setBuddy( arrComboBox1 );
	textLabel2->setBuddy( arrComboBox2 );
	textLabel3->setBuddy( arrComboBox3 );
   
}

////////////////////////////////////////////////////////////////////////
OverlayForm::OverlayForm(  QTabWidget * bar, QWidget* parent ,int tab_number , const char* name, bool modal, WFlags f  ): 
	QDialog( parent, name, modal, f ),
	myBar(bar)
{
    setCaption( "Planar Maps  --  Overlay" );
    resize( 590, 390 );

    optionsFormLayout = new QVBoxLayout( this, 11, 6 );

    split = new QSplitter(this);
    listBox1 = new DDListBox( split );
	listBox2 = new DDListBox( split );
	QString traits;
	Qt_widget_base_tab	*w_demo_p;
	for (int i=0; i < tab_number; i++)
	{
		if ( myBar->isTabEnabled( myBar->page(i) ) )
		{
			// We peform downcasting from QWigdet* to Qt_widget_base_tab*, as we know that only
			// Qt_widget_base_tab objects are stored in the tab pages.
			w_demo_p = static_cast<Qt_widget_base_tab *> (myBar->page(i));
			switch ( w_demo_p->traits_type ) {
			case SEGMENT_TRAITS:
				traits = " ( segment traits )";
				break;
			case POLYLINE_TRAITS:
				traits = " ( polyline traits )";
				break;
			case CONIC_TRAITS:
				traits = " ( conic traits )";
				break;
			}
			listBox1->insertItem( QPixmap( green_icon ) , myBar->label(i) + traits );
		}
	}
	listBox1->set_max_items(listBox1->count());

	arrLayout = new QHBoxLayout();

    textLabel1 = new QLabel( "Possible Planar Maps", this );
    arrLayout->addWidget( textLabel1 );

	textLabel2 = new QLabel( "Chosen Planar Maps", this );
	arrLayout->addWidget( textLabel2 );

	buttonsLayout = new QHBoxLayout( 0, 0, 6 );
    okPushButton = new QPushButton( "OK", this );
    okPushButton->setDefault( TRUE );
    buttonsLayout->addWidget( okPushButton );

    cancelPushButton = new QPushButton( "Cancel", this );
    buttonsLayout->addWidget( cancelPushButton );

    optionsFormLayout->addLayout( arrLayout );
	optionsFormLayout->addWidget( split );
	optionsFormLayout->addLayout( buttonsLayout );

    connect( okPushButton, SIGNAL( clicked() ), this, SLOT( accept() ) );
    connect( cancelPushButton, SIGNAL( clicked() ), this, SLOT( reject() ) );
   
	setAcceptDrops(TRUE);
}
//////////////////////////////////////////////////////////////////////////////
DDListBox::DDListBox( QWidget * parent, const char * name, WFlags f ) :
    QListBox( parent, name, f ),
	max_items(0),
	flag(false)
{
    setAcceptDrops( TRUE );
    dragging = FALSE;
}


void DDListBox::dragEnterEvent( QDragEnterEvent *evt )
{
    if (  QTextDrag::canDecode( evt ) ) 
	evt->accept();
}


void DDListBox::dropEvent( QDropEvent *evt )
{
    QString text;
	if (  QTextDrag::decode( evt, text ) ) 
		insertItem( QPixmap( green_icon ) , text );
	if (count() == max_items && max_items != 0)
	{
		flag = true;
		for (unsigned int i = 0; i < count(); i++)
			item(i)->setSelectable( true );
	}
}


void DDListBox::mousePressEvent( QMouseEvent *evt )
{	
	QListBox::mousePressEvent( evt );
    dragging = TRUE;
}


void DDListBox::mouseMoveEvent( QMouseEvent * )
{	
	if (count() == max_items && max_items != 0 && flag)
	{
		for (unsigned int i = 0; i < count(); i++)
			changeItem( QPixmap( green_icon ) , text(i) , i);
		flag = false;
	}

    if ( dragging && item(currentItem())->isSelectable() ) 
	{
		QDragObject *d = new  QTextDrag( currentText() , this );
		d->dragCopy(); // do NOT delete d.
		dragging = FALSE;
		unsigned int current = currentItem();
		if (count() == max_items && max_items != 0)
		{
			char s[100];
			strcpy(s, currentText()); 
			char * traits;
			traits = strtok(s," ");
			traits = strtok(NULL, " ");
			traits = strtok(NULL, " ");
			traits = strtok(NULL, " ");
			
			for (unsigned int i = 0; i < max_items; i++)
			{
				char s_i[100];
				strcpy(s_i, text(i)); 
				char * traits_i;
				traits_i = strtok(s_i," ");
				traits_i = strtok(NULL, " ");
				traits_i = strtok(NULL, " ");
				traits_i = strtok(NULL, " ");
				bool b = (strcmp(traits,traits_i) == 0);
				if (!b && i != current )
				{
					changeItem( QPixmap( white_icon ) , text(i) , i);
					item(i)->setSelectable( b ); 
				}
			}
		}
				
		removeItem ( current );
    }
}

void DDListBox::set_max_items(int num)
{
	max_items = num;
}


