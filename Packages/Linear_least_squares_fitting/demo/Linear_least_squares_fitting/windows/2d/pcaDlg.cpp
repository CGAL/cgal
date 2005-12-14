// pcaDlg.cpp : implementation file
//

#include "stdafx.h"
#include "pca.h"
#include "pcaDlg.h"
#include ".\pcadlg.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()


// CpcaDlg dialog



CpcaDlg::CpcaDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CpcaDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);

	// init opengl
	m_hGLContext = NULL;
	m_GLPixelIndex = 0;

	// init mouse
	m_RightButtonDown = FALSE;
	m_LeftButtonDown = FALSE;
}

void CpcaDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CpcaDlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	ON_WM_DESTROY()
	ON_WM_CREATE()
	ON_WM_SIZE()
	ON_COMMAND(ID_SET_CLEAR, OnSetClear)
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_MOUSEMOVE()
  ON_COMMAND(ID_FIT_LINE, OnFitLine)
END_MESSAGE_MAP()


// CpcaDlg message handlers

BOOL CpcaDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here
	
	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CpcaDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CpcaDlg::OnPaint() 
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CPaintDC dc(this); // device context for painting
		render();
		SwapBuffers(dc.m_ps.hdc);
		CDialog::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CpcaDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

void CpcaDlg::OnDestroy()
{
	CDialog::OnDestroy();

	if(wglGetCurrentContext() != NULL)
		wglMakeCurrent(NULL,NULL);

	if(m_hGLContext != NULL)
	{
		wglDeleteContext(m_hGLContext);
		m_hGLContext = NULL;
	}
}

BOOL CpcaDlg::SetWindowPixelFormat(HDC hDC)
{
	PIXELFORMATDESCRIPTOR pixelDesc;

	pixelDesc.nSize = sizeof(PIXELFORMATDESCRIPTOR);
	pixelDesc.nVersion = 1;
	pixelDesc.dwFlags = PFD_DRAW_TO_WINDOW | 
						PFD_SUPPORT_OPENGL |
						PFD_DOUBLEBUFFER |
						PFD_STEREO_DONTCARE;
	pixelDesc.iPixelType = PFD_TYPE_RGBA;
	pixelDesc.cColorBits = 32;
	pixelDesc.cRedBits = 8;
	pixelDesc.cRedShift = 16;
	pixelDesc.cGreenBits = 8;
	pixelDesc.cGreenShift = 8;
	pixelDesc.cBlueBits = 8;
	pixelDesc.cBlueShift = 0;
	pixelDesc.cAlphaBits = 0;
	pixelDesc.cAlphaShift = 0;
	pixelDesc.cAccumBits = 64;
	pixelDesc.cAccumRedBits = 16;
	pixelDesc.cAccumGreenBits = 16;
	pixelDesc.cAccumBlueBits = 16;
	pixelDesc.cAccumAlphaBits = 0;
	pixelDesc.cDepthBits = 32;
	pixelDesc.cStencilBits = 8;
	pixelDesc.cAuxBuffers = 0;
	pixelDesc.iLayerType = PFD_MAIN_PLANE;
	pixelDesc.bReserved = 0;
	pixelDesc.dwLayerMask = 0;
	pixelDesc.dwVisibleMask = 0;
	pixelDesc.dwDamageMask = 0;

	m_GLPixelIndex = ChoosePixelFormat(hDC,&pixelDesc);
	if(m_GLPixelIndex == 0) // Choose default
	{
		m_GLPixelIndex = 1;
		if(DescribePixelFormat(hDC,m_GLPixelIndex,
			sizeof(PIXELFORMATDESCRIPTOR),&pixelDesc)==0)
			return FALSE;
	}

	if(SetPixelFormat(hDC,m_GLPixelIndex,&pixelDesc)==FALSE)
		return FALSE;

	return TRUE;
}


BOOL CpcaDlg::CreateViewGLContext(HDC hDC)
{
	m_hGLContext = wglCreateContext(hDC);
	if(m_hGLContext==NULL)
		return FALSE;
	if(wglMakeCurrent(hDC,m_hGLContext)==FALSE)
		return FALSE;
	return TRUE;
}


int CpcaDlg::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if(CDialog::OnCreate(lpCreateStruct) == -1)
		return -1;

	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);

	if(SetWindowPixelFormat(hDC)==FALSE)
		return 0;
	if(CreateViewGLContext(hDC)==FALSE)
		return 0;

	/*
	COLORREF color = ::GetSysColor(COLOR_3DFACE);
	glClearColor((float)GetRValue(color)/255.0f,
				       (float)GetGValue(color)/255.0f,
				       (float)GetBValue(color)/255.0f,
				        1.0);
	*/
	glClearColor(1.0f,1.0f,1.0f,1.0f);
	glPolygonMode(GL_FRONT,GL_FILL);

	return 0;
}

void CpcaDlg::OnSize(UINT nType, int cx, int cy)
{
	CDialog::OnSize(nType, cx, cy);
	glViewport(0,0,(GLsizei)cx,(GLsizei)cy);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0,1.0,0.0,1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glDrawBuffer(GL_BACK);
	InvalidateRect(NULL,FALSE);
}

void CpcaDlg::render()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// draw point set
	glColor3ub(0,0,0);
	glPointSize(5.0f);
	std::list<Point_2>::iterator it;
	for(it = m_points.begin();
		  it != m_points.end();
			it++)
	{
		const Point_2& p = *it;
    draw_disc(p,0.003,0,0,0);
	}

  // draw fitting line
	glColor3ub(150,0,0);
  const Point_2& a = m_fitting_line.point(-100);
  const Point_2& b = m_fitting_line.point(100);
	glBegin(GL_LINES);
		glVertex2d((double)a.x(),(double)a.y());
		glVertex2d((double)b.x(),(double)b.y());
	glEnd();

  draw_disc(m_centroid,0.01,150,0,0);

	glFlush();
}

void CpcaDlg::draw_disc(const Point_2& c,
                        const FT radius,
                        unsigned char r,
                        unsigned char g,
                        unsigned char b)
{
  GLUquadricObj *pQuadric = gluNewQuadric();
  glPushMatrix();
  glTranslated(c.x(),c.y(),0.0);
  gluDisk(pQuadric,0,radius,30,1);
  glPopMatrix();
  gluDeleteQuadric(pQuadric);
}

void CpcaDlg::OnSetClear()
{
	m_points.clear();
	InvalidateRect(NULL,FALSE);
}

void CpcaDlg::OnLButtonDown(UINT nFlags, CPoint point)
{
	m_LeftButtonDown = TRUE;
	SetCapture();
	CDialog::OnLButtonDown(nFlags, point);
}

void CpcaDlg::OnLButtonUp(UINT nFlags, CPoint point)
{
	m_LeftButtonDown = FALSE;
	ReleaseCapture();

	// add point to point set
	m_points.push_back(convert(point));
  OnFitLine();
	InvalidateRect(NULL,FALSE);

	CDialog::OnLButtonUp(nFlags, point);
}

void CpcaDlg::OnMouseMove(UINT nFlags, CPoint point)
{
  if(!inside(point))
    return;

	if(m_LeftButtonDown)
	{
		m_points.push_back(convert(point));
    //OnFitLine();
		InvalidateRect(NULL,FALSE);
	}

	CDialog::OnMouseMove(nFlags, point);
}

Point_2 CpcaDlg::convert(const CPoint& point)
{
	CRect rect;
	GetClientRect(rect);
	FT cx = (FT)point.x;
	FT cy = (FT)point.y;
	FT w = (FT)rect.Width();
	FT h = (FT)rect.Height();
	FT x = cx/w;
	FT y = 1.0f-cy/h;
	return Point_2(x,y);
}

bool CpcaDlg::inside(const CPoint& point)
{
	CRect rect;
	GetClientRect(rect);
	FT cx = (FT)point.x;
	FT cy = (FT)point.y;
	FT w = (FT)rect.Width();
	FT h = (FT)rect.Height();
	FT x = cx/w;
	FT y = 1.0f-cy/h;
  if(x > 0.0f && x < 1.0f && 
     y > 0.0f && y < 1.0f)
     return true;
	return false;
}


void CpcaDlg::OnFitLine()
{
  FT quality = linear_least_squares_fitting_2(m_points.begin(),
                                              m_points.end(),
                                              m_fitting_line,
                                              m_centroid);
}
