// StaticView.cpp : implementation file
//

#include "stdafx.h"
#include "Mesh.h"
#include "StaticView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CStaticView

CStaticView::CStaticView()
{
	m_xRotation = 0.0f;
	m_yRotation = 0.0f;
	m_zRotation = 0.0f;
	m_xTranslation = 0.0f;
	m_yTranslation = 0.0f;
	m_zTranslation = -5.0f;
	m_xScaling = 1.0f;
	m_yScaling = 1.0f;
	m_zScaling = 1.0f;
	m_SpeedRotation = 1.0f / 3.0f;
	m_SpeedTranslation = 1.0f / 50.0f;
	m_xyRotation = 1;
}

CStaticView::~CStaticView()
{
}

BEGIN_MESSAGE_MAP(CStaticView, CStatic)
	//{{AFX_MSG_MAP(CStaticView)
	ON_WM_PAINT()
	ON_WM_DESTROY()
	ON_WM_SIZE()
	ON_WM_TIMER()
	ON_WM_CREATE()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CStaticView message handlers


void CStaticView::OnPaint() 
{
	CPaintDC dc(this); 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	glPushMatrix();
 
		glTranslated(m_xTranslation,m_yTranslation,m_zTranslation);
		glRotatef(m_xRotation, 1.0, 0.0, 0.0);
		glRotatef(m_yRotation, 0.0, 1.0, 0.0);
		glRotatef(m_zRotation, 0.0, 0.0, 1.0);
		glScalef(m_xScaling,m_yScaling,m_zScaling);
		//glCallList(1);
		BuildListCube();

	glPopMatrix();

	// Double buffer
	SwapBuffers(dc.m_ps.hdc);
	glFlush();
}

//********************************************
// SetWindowPixelFormat
//********************************************
BOOL CStaticView::SetWindowPixelFormat(HDC hDC)
{
	PIXELFORMATDESCRIPTOR pixelDesc;
	
	pixelDesc.nSize = sizeof(PIXELFORMATDESCRIPTOR);
	pixelDesc.nVersion = 1;
	
	pixelDesc.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL |
		PFD_DOUBLEBUFFER | PFD_STEREO_DONTCARE;
	
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
	if(m_GLPixelIndex == 0) 
	{
		m_GLPixelIndex = 1;
		if(DescribePixelFormat(hDC,m_GLPixelIndex,
			sizeof(PIXELFORMATDESCRIPTOR),&pixelDesc)==0)
			return FALSE;
	}
	if(!SetPixelFormat(hDC,m_GLPixelIndex,&pixelDesc))
		return FALSE;
	return TRUE;
}



//********************************************
// CreateViewGLContext
// Create an OpenGL rendering context
//********************************************
BOOL CStaticView::CreateViewGLContext(HDC hDC)
{
	m_hGLContext = wglCreateContext(hDC);
	
	if(m_hGLContext==NULL)
		return FALSE;
	
	if(wglMakeCurrent(hDC,m_hGLContext)==FALSE)
		return FALSE;
	
	return TRUE;
}

void CStaticView::OnDestroy() 
{
	if(wglGetCurrentContext() != NULL)
		wglMakeCurrent(NULL,NULL);
	
	if(m_hGLContext != NULL)
	{
		wglDeleteContext(m_hGLContext);
		m_hGLContext = NULL;
	}
	CStatic::OnDestroy();
}

void CStaticView::OnSize(UINT nType, int cx, int cy) 
{
	CStatic::OnSize(nType, cx, cy);
	
	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	wglMakeCurrent(hDC,m_hGLContext);

	double aspect = (cy == 0) ? 
	  cx : (double)cx/(double)cy;
	
	glViewport(0,0,cx,cy);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45,aspect,0.1,1000.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	::ReleaseDC(hWnd,hDC);
}

//***********************************************
// BuildList
//***********************************************
void CStaticView::BuildListCube()
{
	GLUquadricObj* pQuadric = gluNewQuadric();

	float x = 1.0f;
  glColor3ub(0,0,0);

	glBegin(GL_POLYGON);
		glVertex3d( x,  x, x);
		glVertex3d( x, -x, x);
		glVertex3d(-x, -x, x);
		glVertex3d(-x,  x, x);
	glEnd();

	glBegin(GL_POLYGON);
		glVertex3d(-x,  x, -x);
		glVertex3d(-x, -x, -x);
		glVertex3d( x, -x, -x);
		glVertex3d( x,  x, -x);
	glEnd();

	glBegin(GL_POLYGON);
		glVertex3d( x,  x,  x);
		glVertex3d( x,  x, -x);
		glVertex3d( x, -x, -x);
		glVertex3d( x, -x,  x);
	glEnd();

	glBegin(GL_POLYGON);
		glVertex3d(-x,  x,  x);
		glVertex3d(-x,  x, -x);
		glVertex3d(-x, -x, -x);
		glVertex3d(-x, -x,  x);
	glEnd();

	glBegin(GL_POLYGON);
		glVertex3d(-x, -x,  x);
		glVertex3d( x, -x,  x);
		glVertex3d( x, -x, -x);
		glVertex3d(-x, -x, -x);
	glEnd();


	glBegin(GL_POLYGON);
		glVertex3d(-x,  x,  x);
		glVertex3d( x,  x,  x);
		glVertex3d( x,  x, -x);
		glVertex3d(-x,  x, -x);
	glEnd();

	gluDeleteQuadric(pQuadric);
}

void CStaticView::OnTimer(UINT nIDEvent) 
{
	m_xRotation += 0.3f;	
	m_yRotation += 3.0f;	
	InvalidateRect(NULL,FALSE);
	CStatic::OnTimer(nIDEvent);
}

int CStaticView::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CStatic::OnCreate(lpCreateStruct) == -1)
		return -1;
	return 0;
}

//***************************************
// Init
//***************************************
int CStaticView::Init() 
{ 
	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);

	if(SetWindowPixelFormat(hDC)==FALSE)
		return FALSE;
	
	if(CreateViewGLContext(hDC)==FALSE)
		return FALSE;

  glShadeModel(GL_FLAT);
	glEnable(GL_NORMALIZE);

	// Lights properties
  float	ambientProperties[]  = {0.7f, 0.7f, 0.7f, 1.0f};
	float	diffuseProperties[]  = {0.8f, 0.8f, 0.8f, 1.0f};
  float	specularProperties[] = {1.0f, 1.0f, 1.0f, 1.0f};
	
  glLightfv( GL_LIGHT0, GL_AMBIENT, ambientProperties);
  glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuseProperties);
  glLightfv( GL_LIGHT0, GL_SPECULAR, specularProperties);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);

	COLORREF color = GetSysColor(COLOR_3DFACE);
	float r = (float)GetRValue(color)/255.0f;
	float g = (float)GetGValue(color)/255.0f;
	float b = (float)GetBValue(color)/255.0f;
	glClearColor(r,g,b,1.0f);
  glClearDepth(1.0f);
	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

	// Perspective
	CRect rect;
	GetClientRect(&rect);
	double aspect = (rect.Height() == 0) ? 
		rect.Width() : (double)rect.Width()/(double)rect.Height();
	gluPerspective(45,aspect,0.1,1000.0);

	glPolygonMode(GL_FRONT,GL_LINE);
	glPolygonMode(GL_BACK,GL_LINE);

	// Antialiasing
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
	glLineWidth(1.1f);

	SetTimer(1,1,NULL);
	return 1;
}
