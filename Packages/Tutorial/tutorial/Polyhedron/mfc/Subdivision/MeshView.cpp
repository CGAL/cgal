// MeshView.cpp : implementation of the CMeshView class
//

#include "stdafx.h"
#include "Mesh.h"
#include "MeshDoc.h"
#include "MeshView.h"
#include "Lib/dump_eps.h"
#include "Lib/texture.h"
#include "DialogRenderingMode.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CMeshView

IMPLEMENT_DYNCREATE(CMeshView, CView)

BEGIN_MESSAGE_MAP(CMeshView, CView)
	// Standard printing commands
	ON_COMMAND(ID_FILE_PRINT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, CView::OnFilePrintPreview)
	ON_WM_CREATE()
	ON_WM_PAINT()
	ON_WM_SIZE()
	ON_WM_ACTIVATE()
	ON_WM_ERASEBKGND()
	ON_WM_DESTROY()
	ON_WM_MOUSEMOVE()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONUP()
	ON_COMMAND(ID_RENDER_LIGHT, OnRenderLight)
	ON_UPDATE_COMMAND_UI(ID_RENDER_LIGHT, OnUpdateRenderLight)
	ON_COMMAND(ID_RENDER_CULLING, OnRenderCulling)
	ON_UPDATE_COMMAND_UI(ID_RENDER_CULLING, OnUpdateRenderCulling)
	ON_COMMAND(ID_RENDER_BACKCOLOR, OnRenderBackcolor)
	ON_COMMAND(ID_EDIT_COPY, OnEditCopy)
	ON_COMMAND(ID_MODE_FILL, OnModeFill)
	ON_UPDATE_COMMAND_UI(ID_MODE_FILL, OnUpdateModeFill)
	ON_COMMAND(ID_MODE_WIREFRAME, OnModeWireframe)
	ON_UPDATE_COMMAND_UI(ID_MODE_WIREFRAME, OnUpdateModeWireframe)
	ON_COMMAND(ID_MODE_POINT, OnModePoint)
	ON_UPDATE_COMMAND_UI(ID_MODE_POINT, OnUpdateModePoint)
	ON_COMMAND(ID_RENDER_SMOOTH, OnRenderSmooth)
	ON_UPDATE_COMMAND_UI(ID_RENDER_SMOOTH, OnUpdateRenderSmooth)
	ON_COMMAND(ID_VIEW_ALL, OnViewAll)
	ON_COMMAND(ID_RENDER_SUPERIMPOSE_EDGES, OnRenderSuperimposeEdges)
	ON_UPDATE_COMMAND_UI(ID_RENDER_SUPERIMPOSE_EDGES, OnUpdateRenderSuperimposeEdges)
	ON_COMMAND(ID_COLORS_MESH, OnColorsMesh)
	ON_COMMAND(ID_COLORS_SUPERIMPOSEDEDGES, OnColorsSuperimposededges)
	ON_COMMAND(ID_RENDER_ANTIALIASING, OnRenderAntialiasing)
	ON_UPDATE_COMMAND_UI(ID_RENDER_ANTIALIASING, OnUpdateRenderAntialiasing)
	ON_COMMAND(ID_PREDEFINED_MESH, OnPredefinedmodesMesh)
	ON_UPDATE_COMMAND_UI(ID_PREDEFINED_MESH, OnUpdatePredefinedmodesMesh)
	ON_COMMAND(ID_PREDEFINED_CONTROLMESH, OnPredefinedControlmesh)
	ON_UPDATE_COMMAND_UI(ID_PREDEFINED_CONTROLMESH, OnUpdatePredefinedControlmesh)
	ON_COMMAND(ID_PREDEFINED_THICK, OnPredefinedMesh)
	ON_UPDATE_COMMAND_UI(ID_PREDEFINED_THICK, OnUpdatePredefinedMesh)
	ON_COMMAND(ID_COLORS_THICKEDGES, OnColorsThickedges)
	ON_COMMAND(ID_RENDER_VISUALCHOOSER, OnRenderVisualchooser)
	ON_COMMAND(ID_SUPERIMPOSE_VERTICES, OnSuperimposeVertices)
	ON_UPDATE_COMMAND_UI(ID_SUPERIMPOSE_VERTICES, OnUpdateSuperimposeVertices)
	ON_COMMAND(ID_COLORS_VERTICES, OnColorsVertices)
	ON_COMMAND(ID_MATERIAL_PEARL, OnMaterialPearl)
	ON_COMMAND(ID_MATERIAL_BRASS, OnMaterialBrass)
	ON_COMMAND(ID_MATERIAL_BLACKPLASTIC, OnMaterialBlackplastic)
	ON_COMMAND(ID_MATERIAL_GOLD, OnMaterialGold)
	ON_COMMAND(ID_MATERIAL_SILVER, OnMaterialSilver)
	ON_COMMAND(ID_MATERIAL_JADE, OnMaterialJade)
	ON_COMMAND(ID_MATERIAL_RUBY, OnMaterialRuby)
	ON_COMMAND(ID_FILE_DUMPTOPS, OnFileDumptops)
	ON_UPDATE_COMMAND_UI(ID_FILE_DUMPTOPS, OnUpdateFileDumptops)
	ON_COMMAND(ID_BOUNDIXBOX_SHOW, OnBoundixboxShow)
	ON_UPDATE_COMMAND_UI(ID_BOUNDIXBOX_SHOW, OnUpdateBoundixboxShow)
	ON_COMMAND(ID_BOUNDIXBOX_SHOWWHENMOVING, OnBoundixboxShowwhenmoving)
	ON_UPDATE_COMMAND_UI(ID_BOUNDIXBOX_SHOWWHENMOVING, OnUpdateBoundixboxShowwhenmoving)
	ON_COMMAND(ID_SPECIAL_DRAWVORONOIEDGES, OnSpecialDrawvoronoiedges)
	ON_UPDATE_COMMAND_UI(ID_SPECIAL_DRAWVORONOIEDGES, OnUpdateSpecialDrawvoronoiedges)
END_MESSAGE_MAP()

// CMeshView construction/destruction

CMeshView::CMeshView()
{
	// OpenGL
	m_hGLContext = NULL;
	m_GLPixelIndex = 0;

	// rendering options
	m_Lighting = true;
	m_Culling = false;
	m_FirstView = true;
	m_UseNormals = true;
	m_Antialiasing = false;
	m_SmoothShading = true;
	m_PolygonMode = GL_FILL;
	m_SuperimposeEdges = true;
	m_DrawVoronoiEdges = false;
	m_DrawBoundingBox = false;
	m_SuperimposeVertices = false;
	m_ThickerControlEdges = false;
	m_ThicknessControlEdges = 3.0f;
	m_DrawBoundingBoxWhenMoving = false;
	m_SuperimposeOnlyControlMesh = true;

	// other options
	m_PointSize = 3.0f;
	m_BackColor[0] = 	m_BackColor[1] = m_BackColor[2] = 0.0f;
	m_MeshColor[0] = 	m_MeshColor[1] = m_MeshColor[2] = 1.0f;
	m_EdgeColor[0] = 	m_EdgeColor[1] = m_EdgeColor[2] = 0.0f;
	m_VertexColor[0] = m_VertexColor[1] = m_VertexColor[2] = 0.0;
	m_ControlEdgeColor[0] = m_ControlEdgeColor[1] = m_ControlEdgeColor[2] = 0.0f;

	// mouse
	m_LeftButtonDown = false;
	m_RightButtonDown = false;
	m_Moving = false;
}

CMeshView::~CMeshView()
{
}

BOOL CMeshView::PreCreateWindow(CREATESTRUCT& cs)
{
	return CView::PreCreateWindow(cs);
}

// CMeshView drawing
// OnPaint message is preferred
void CMeshView::OnDraw(CDC* /*pDC*/)
{
}


// CMeshView printing
BOOL CMeshView::OnPreparePrinting(CPrintInfo* pInfo)
{
	// default preparation
	return DoPreparePrinting(pInfo);
}

void CMeshView::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
}
void CMeshView::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
}

// CMeshView diagnostics

#ifdef _DEBUG
void CMeshView::AssertValid() const
{
	CView::AssertValid();
}

void CMeshView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CMeshDoc* CMeshView::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CMeshDoc)));
	return (CMeshDoc*)m_pDocument;
}
#endif //_DEBUG

//////////////////////////////////////////////
// OPENGL
//////////////////////////////////////////////

//********************************************
// SetWindowPixelFormat
//********************************************
BOOL CMeshView::SetWindowPixelFormat(HDC hDC)
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
	if(m_GLPixelIndex == 0) // Choose default
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
BOOL CMeshView::CreateViewGLContext(HDC hDC)
{
	m_hGLContext = wglCreateContext(hDC);

	if(m_hGLContext==NULL)
		return FALSE;

	if(wglMakeCurrent(hDC,m_hGLContext)==FALSE)
		return FALSE;

	return TRUE;	
}

//********************************************
// Create an OpenGL framework
//********************************************
int CMeshView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if(CView::OnCreate(lpCreateStruct) == -1)
		return -1;

	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);

	if(SetWindowPixelFormat(hDC)==FALSE)
		return -1;
	
	if(CreateViewGLContext(hDC)==FALSE)
		return -1;

	// Default mode
	glEnable(GL_DEPTH_TEST);
	glClearColor(m_BackColor[0],
		           m_BackColor[1],
							 m_BackColor[2],
							 1.0f);

	// lighting
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	float lpos[4] = { -.2f, .2f, .9797958971f, 0.0f };
	glLightfv(GL_LIGHT0,GL_POSITION,lpos);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0f);

	// back material
  float MatAmbientBack[]  = {0.0f, 1.0f, 0.0f, 1.0f};
  glMaterialfv(GL_BACK,GL_AMBIENT,MatAmbientBack);
	ChangeMaterial(CString("Silver"),false);

	// init camera
	InitCamera();

	return 0;
}

void CMeshView::OnPaint()
{
	// device context for painting
	CPaintDC dc(this); 

	// model is stored in Document
	CMeshDoc *pDoc = (CMeshDoc *)GetDocument();
	if(pDoc->m_pMesh == NULL)
		return;

	// setup camera (only once)
	ViewAll();

	// useful in multidoc templates
	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	//wglMakeCurrent(hDC,m_hGLContext);
	wglMakeCurrent(dc.m_ps.hdc,m_hGLContext);
	

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// shading option
	if(m_SmoothShading)
		glShadeModel(GL_SMOOTH);
	else
		glShadeModel(GL_FLAT);

	// culling option
	if(m_Culling)
		glEnable(GL_CULL_FACE);
	else
		glDisable(GL_CULL_FACE);

	// polygon mode (point, line or fill)
	glPolygonMode(GL_FRONT_AND_BACK,m_PolygonMode);

	// set mesh color
	glColor3f(m_MeshColor[0],m_MeshColor[1],m_MeshColor[2]);

	// antialiasing
	if(m_Antialiasing)
	{
		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
		glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
		glLineWidth(1.5f);
	}
	else
	{
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);
		glLineWidth(1.0f);
	}

	// drawing
	glPushMatrix();

	  // setup viewpoint from current arcball
		m_Arcball.glDraw();

		if(m_SuperimposeEdges || m_SuperimposeVertices)
		{
			// enable polygon offset
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(3.0f,1.0f);
		}

		// draw the mesh 
		if((m_Moving && m_DrawBoundingBoxWhenMoving) || 
			 m_DrawBoundingBox)
		{
			glColor3f(1.0f,0.0f,0.0f);
			glDisable(GL_LIGHTING);
		  pDoc->m_pMesh->gl_draw_bounding_box();
		}

		// lighting option
		if(m_Lighting)
		{
			m_UseNormals = true;
			glEnable(GL_LIGHTING);
		}
		else
		{
			m_UseNormals = false;
			glDisable(GL_LIGHTING);
		}

		if(!m_Moving || !m_DrawBoundingBoxWhenMoving)
			pDoc->m_pMesh->gl_draw(m_SmoothShading,
			                       m_UseNormals);

		// disable lighting
		if(m_SuperimposeEdges || m_SuperimposeVertices)
			glDisable(GL_LIGHTING);

		// draw the mesh once again with a few
		// options desactivated
		if(m_SuperimposeEdges && 
			 !(m_Moving && m_DrawBoundingBoxWhenMoving))
		{
			// set line mode
			glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);

			// edge color
			glColor3f(m_EdgeColor[0],m_EdgeColor[1],m_EdgeColor[2]);

			// superimpose edges on the mesh
			bool skip_ordinary_edges = m_SuperimposeOnlyControlMesh;
			bool skip_control_edges = m_ThickerControlEdges;
			pDoc->m_pMesh->superimpose_edges(skip_ordinary_edges,
			                                 skip_control_edges,
																			 m_DrawVoronoiEdges);

			if(m_ThickerControlEdges)
			{
				glLineWidth(m_ThicknessControlEdges);
				// thick edge color
				glColor3f(m_ControlEdgeColor[0],
					        m_ControlEdgeColor[1],
									m_ControlEdgeColor[2]);
				skip_ordinary_edges = true;
				skip_control_edges = false;
				pDoc->m_pMesh->superimpose_edges(skip_ordinary_edges,
			                                   skip_control_edges,
																				 m_DrawVoronoiEdges);
			}
		} // end superimpose edges

		// superimpose vertices
		if(m_SuperimposeVertices &&
			 !(m_Moving && m_DrawBoundingBoxWhenMoving))
		{
			glColor3f(m_VertexColor[0],m_VertexColor[1],m_VertexColor[2]);
			pDoc->m_pMesh->superimpose_spheres(0.1);
		} // end superimpose vertices

	// disable polygon offset
	if(m_SuperimposeEdges || m_SuperimposeVertices)
		glDisable(GL_POLYGON_OFFSET_FILL);

	glPopMatrix();

	// Double buffer
	SwapBuffers(dc.m_ps.hdc);
	glFlush();

	// release
	::ReleaseDC(hWnd,hDC);
}

// resize client
void CMeshView::OnSize(UINT nType, int cx, int cy)
{
	CView::OnSize(nType, cx, cy);
	
	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	wglMakeCurrent(hDC,m_hGLContext);

	// set OpenGL viewport and perspective
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	m_Viewport.SetSize(cx,cy);
	m_Viewport.glDraw();
	m_Camera.glDraw(m_Viewport);
	glDrawBuffer(GL_BACK);
	::ReleaseDC(hWnd,hDC);
}

void CMeshView::OnActivateView(BOOL bActivate, 
															 CView* pActivateView, 
															 CView* pDeactiveView) 
{
	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	wglMakeCurrent(hDC,m_hGLContext);
	::ReleaseDC(hWnd,hDC);
	CView::OnActivateView(bActivate, 
		                    pActivateView, 
												pDeactiveView);
}

// avoid flickering
BOOL CMeshView::OnEraseBkgnd(CDC*)
{
	return TRUE;
}

// destroy rendering context
void CMeshView::OnDestroy()
{
	if(wglGetCurrentContext() != NULL)
		wglMakeCurrent(NULL,NULL);

	if(m_hGLContext != NULL)
	{
		wglDeleteContext(m_hGLContext);
		m_hGLContext = NULL;
	}
}

// handle mouse button
void CMeshView::HandleMouseButton(int x, 
																	int y)
{
  CVector3d vec = m_Arcball.Intersect(x,m_Viewport.yRes()-y,
		                                  m_Camera,m_Viewport);
  m_Arcball.EndDrag(vec);
	m_Arcball.SetMode(m_LeftButtonDown+2*m_RightButtonDown);
  vec = m_Arcball.Intersect(x,m_Viewport.yRes()-y,m_Camera,m_Viewport);
	m_Arcball.BeginDrag(vec);
}

// mouse move
void CMeshView::OnMouseMove(UINT nFlags, CPoint point)
{
	if(m_LeftButtonDown || m_RightButtonDown)
	{
		CVector3d vec = m_Arcball.Intersect(point.x,m_Viewport.yRes()-point.y,
										                    m_Camera,m_Viewport);
		m_Arcball.Motion(vec);
		m_Moving = true;
		InvalidateRect(NULL,FALSE);
	}
	else
		m_Moving = false;
}

void CMeshView::OnLButtonDown(UINT nFlags, CPoint point)
{
	SetCapture();
	m_LeftButtonDown = TRUE;
	HandleMouseButton(point.x,point.y);
	CView::OnLButtonDown(nFlags, point);
}

void CMeshView::OnLButtonUp(UINT nFlags, CPoint point)
{
	ReleaseCapture();
	m_LeftButtonDown = FALSE;
	m_Moving = false;
	HandleMouseButton(point.x,point.y);
	CView::OnLButtonUp(nFlags, point);
	InvalidateRect(NULL,FALSE);
}

void CMeshView::OnRButtonDown(UINT nFlags, CPoint point)
{
	SetCapture();
	m_RightButtonDown = TRUE;
	HandleMouseButton(point.x,point.y);
	CView::OnRButtonDown(nFlags, point);
}

void CMeshView::OnRButtonUp(UINT nFlags, CPoint point)
{
	ReleaseCapture();
	m_RightButtonDown = FALSE;
	m_Moving = false;
	HandleMouseButton(point.x,point.y);
	CView::OnRButtonUp(nFlags, point);
	InvalidateRect(NULL,FALSE);
}


void CMeshView::OnEditCopy()
{
  // Clean clipboard of contents, and copy the DIB.
  if(OpenClipboard())
   {
		CMeshDoc *pDoc = (CMeshDoc *)GetDocument();
		pDoc->StatusMessage("Bmp output...");
    BeginWaitCursor();

		// Snap
		CRect rect;
		GetClientRect(&rect);
		CSize size(rect.Width(),rect.Height());
		pDoc->StatusMessage("Bmp output...snap...");

		unsigned char *pixel = new unsigned char[3*size.cx*size.cy];
		if(pixel == NULL)
		{
			AfxMessageBox("insufficient memory");
			return;
		}
		glPixelStorei(GL_PACK_ALIGNMENT,1) ; 
    glPixelStorei(GL_PACK_ROW_LENGTH, size.cx) ; 
    glReadPixels(0,0,size.cx,size.cy,GL_RGB,GL_UNSIGNED_BYTE,pixel); 
    glPixelStorei(GL_PACK_ROW_LENGTH, 0) ; 

		// Link image - buffer
		CTexture image;
		pDoc->StatusMessage("Bmp output...snap...read buffer...");
		VERIFY(image.ReadBuffer(pixel,size.cx,size.cy,24));
		image.BGRtoRGB();

		// Cleanup memory
		delete [] pixel;

		pDoc->StatusMessage("Bmp output...snap...read buffer...fill clipboard...");
    VERIFY(EmptyClipboard());
    VERIFY(SetClipboardData(CF_DIB,image.ExportHandle()));
    VERIFY(CloseClipboard());
    EndWaitCursor();

		pDoc->StatusMessage("Ready");
   }
}

void CMeshView::OnInitialUpdate()
{
	CView::OnInitialUpdate();
}


//////////////////////////////////////////////
// OPENGL OPTIONS
//////////////////////////////////////////////

// ligthing 
void CMeshView::OnRenderLight()
{
	m_Lighting = !m_Lighting;
	CMeshDoc *pDoc = (CMeshDoc *)GetDocument();
	InvalidateRect(NULL,FALSE);
}
void CMeshView::OnUpdateRenderLight(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_Lighting);
}

// culling
void CMeshView::OnRenderCulling()
{
	m_Culling = !m_Culling;
	InvalidateRect(NULL,FALSE);
}
void CMeshView::OnUpdateRenderCulling(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_Culling);
}

// fill rendering mode
void CMeshView::OnModeFill()
{
	m_PolygonMode = GL_FILL;
	InvalidateRect(NULL,FALSE);
}
void CMeshView::OnUpdateModeFill(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_PolygonMode == GL_FILL);
}

// line rendering mode
void CMeshView::OnModeWireframe()
{
	m_PolygonMode = GL_LINE;
	InvalidateRect(NULL,FALSE);
}
void CMeshView::OnUpdateModeWireframe(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_PolygonMode == GL_LINE);
}

// point rendering mode
void CMeshView::OnModePoint()
{
	m_PolygonMode = GL_POINT;
	InvalidateRect(NULL,FALSE);
}
void CMeshView::OnUpdateModePoint(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_PolygonMode == GL_POINT);
}

// smooth vs flat shading
void CMeshView::OnRenderSmooth()
{
	m_SmoothShading = !m_SmoothShading;
	CMeshDoc *pDoc = (CMeshDoc *)GetDocument();
	InvalidateRect(NULL,FALSE);
}
void CMeshView::OnUpdateRenderSmooth(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_SmoothShading);
}

// superimpose edges
void CMeshView::OnRenderSuperimposeEdges()
{
	m_SuperimposeEdges = !m_SuperimposeEdges;
	InvalidateRect(NULL,FALSE);
}
void CMeshView::OnUpdateRenderSuperimposeEdges(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_SuperimposeEdges);
}

// superimpose vertices
void CMeshView::OnSuperimposeVertices()
{
	m_SuperimposeVertices = !m_SuperimposeVertices;
	InvalidateRect(NULL,FALSE);
}
void CMeshView::OnUpdateSuperimposeVertices(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_SuperimposeVertices);
}


void CMeshView::OnPredefinedmodesMesh()
{
	m_SuperimposeEdges = true;
	m_ThickerControlEdges = false;
	m_SuperimposeOnlyControlMesh = false;
	InvalidateRect(NULL,FALSE);
}
void CMeshView::OnUpdatePredefinedmodesMesh(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_SuperimposeEdges && 
		               !m_ThickerControlEdges &&
									 !m_SuperimposeOnlyControlMesh);
}

void CMeshView::OnPredefinedControlmesh()
{
	m_SuperimposeEdges = true;
	m_ThickerControlEdges = false;
	m_SuperimposeOnlyControlMesh = true;
	InvalidateRect(NULL,FALSE);
}

void CMeshView::OnUpdatePredefinedControlmesh(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_SuperimposeEdges && 
		               !m_ThickerControlEdges &&
									 m_SuperimposeOnlyControlMesh);
}

void CMeshView::OnPredefinedMesh()
{
	m_SuperimposeEdges = true;
	m_ThickerControlEdges = true;
	m_SuperimposeOnlyControlMesh = false;
	InvalidateRect(NULL,FALSE);
}

void CMeshView::OnUpdatePredefinedMesh(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_SuperimposeEdges && 
		               m_ThickerControlEdges &&
									 !m_SuperimposeOnlyControlMesh);
}


void CMeshView::OnRenderAntialiasing()
{
	m_Antialiasing = !m_Antialiasing;
	InvalidateRect(NULL,FALSE);
}
void CMeshView::OnUpdateRenderAntialiasing(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_Antialiasing);
}




//////////////////////////////////////////////
// CAMERA
//////////////////////////////////////////////

// view all the current object
void CMeshView::ViewAll(bool check_first)
{
	// make it only once
	if(!m_FirstView && check_first)
		  return;

	CMeshDoc *pDoc = (CMeshDoc *)GetDocument();
	ASSERT(pDoc != NULL);
	if(pDoc->m_pMesh == NULL)
		return;
	m_FirstView = false;

	// set up the camera to visualize the whole object
	TRACE("setup camera.");
	pDoc->m_pMesh->compute_bounding_box();
  CMatrix44 ArcballMatrix = m_Arcball.GetMatrix();

  CVector3d minBound, maxBound;
	minBound.Set(pDoc->m_pMesh->xmin(),
		           pDoc->m_pMesh->ymin(),
							 pDoc->m_pMesh->zmin());
	maxBound.Set(pDoc->m_pMesh->xmax(),
		           pDoc->m_pMesh->ymax(),
							 pDoc->m_pMesh->zmax());
  minBound = ArcballMatrix * minBound;
  maxBound = ArcballMatrix * maxBound;
	m_Camera.ViewAll(minBound[0],
		               maxBound[0],
		               minBound[1],
		               maxBound[1],
		               minBound[2],
		               maxBound[2],
									 m_Viewport);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  m_Camera.glDraw(m_Viewport);
  glMatrixMode(GL_MODELVIEW);
	TRACE("..ok\n");
}

// init camera parameters
void CMeshView::InitCamera()
{
	// perspective
	CRect rect;
	GetClientRect(&rect);

	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);
	wglMakeCurrent(hDC,m_hGLContext);

	// set viewport and camera
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	m_Viewport.SetOrigin(0,0);
	m_Viewport.SetSize(rect.Width(),rect.Height());
	m_Camera.SetHeightAngle(45.);
	m_Camera.SetPosition(0.,0.,5.);
	m_Camera.SetOrientation(0.,1.,0.,0);
	m_Camera.SetNearDistance(.1);
	m_Camera.SetFarDistance(1000.);
	m_Viewport.glDraw();
	m_Camera.glDraw(m_Viewport);
	glDrawBuffer(GL_BACK);

	::ReleaseDC(hWnd,hDC);
}

void CMeshView::OnViewAll()
{
	ViewAll(false);
	InvalidateRect(NULL,FALSE);
}


bool CMeshView::EditColor(float pColor[3])
{
	unsigned char r = (unsigned char)(pColor[0]*255);
	unsigned char g = (unsigned char)(pColor[1]*255);
	unsigned char b = (unsigned char)(pColor[2]*255);
	COLORREF color = RGB(r,g,b);
	CColorDialog dlg(color);
	if(dlg.DoModal() == IDOK)
	{
		color = dlg.GetColor();
		pColor[0] = (float)GetRValue(color)/255.0f;
		pColor[1] = (float)GetGValue(color)/255.0f;
		pColor[2] = (float)GetBValue(color)/255.0f;
		return true;
	}
	return false;
}


//////////////////////////////////////////////
// COLORS
//////////////////////////////////////////////

void CMeshView::OnRenderBackcolor()
{
	if(EditColor(m_BackColor))
	{
		glClearColor(m_BackColor[0],
								 m_BackColor[1],
								 m_BackColor[2],
								 1.0f);
	  InvalidateRect(NULL,FALSE);
	}
}

void CMeshView::OnColorsMesh()
{
	if(EditColor(m_MeshColor))
	  InvalidateRect(NULL,FALSE);
}

void CMeshView::OnColorsSuperimposededges()
{
	if(EditColor(m_EdgeColor))
	  InvalidateRect(NULL,FALSE);
}

void CMeshView::OnColorsThickedges()
{
	if(EditColor(m_ControlEdgeColor))
	  InvalidateRect(NULL,FALSE);
}

void CMeshView::OnColorsVertices()
{
	if(EditColor(m_VertexColor))
	  InvalidateRect(NULL,FALSE);
}


void CMeshView::OnRenderVisualchooser()
{
	CDialogRenderingMode dlg;
	if(dlg.DoModal())
	{
		int mode = dlg.m_Mode;
		switch(mode)
		{
		case 0:
			m_Lighting = false;
			m_SuperimposeOnlyControlMesh = false;
			m_SuperimposeEdges = false;
			m_ThickerControlEdges = false;
			m_PolygonMode = GL_LINE;
			m_Culling = false;
			break;
		case 1:
			m_Lighting = false;
			m_SuperimposeOnlyControlMesh = false;
			m_SuperimposeEdges = false;
			m_ThickerControlEdges = false;
			m_PolygonMode = GL_LINE;
			m_Culling = true;
			break;
		case 2:
			m_Lighting = true;
			m_SuperimposeOnlyControlMesh = false;
			m_SuperimposeEdges = false;
			m_ThickerControlEdges = false;
			m_PolygonMode = GL_FILL;
			m_Culling = false;
			break;
		case 3:
			m_Lighting = true;
			m_SuperimposeOnlyControlMesh = false;
			m_SuperimposeEdges = true;
			m_ThickerControlEdges = false;
			m_PolygonMode = GL_FILL;
			m_Culling = false;
			break;
		case 4:
			m_Lighting = true;
			m_SuperimposeOnlyControlMesh = false;
			m_SuperimposeEdges = true;
			m_ThickerControlEdges = true;
			m_PolygonMode = GL_FILL;
			m_Culling = false;
			break;
		case 5:
			m_Lighting = true;
			m_SuperimposeOnlyControlMesh = true;
			m_SuperimposeEdges = true;
			m_ThickerControlEdges = false;
			m_PolygonMode = GL_FILL;
			m_Culling = false;
			break;
		case 6:
			m_Lighting = false;
			m_SuperimposeOnlyControlMesh = false;
			m_SuperimposeEdges = true;
			m_ThickerControlEdges = false;
			m_PolygonMode = GL_FILL;
			m_Culling = false;
			break;
		case 7:
			m_Lighting = false;
			m_SuperimposeOnlyControlMesh = false;
			m_SuperimposeEdges = true;
			m_ThickerControlEdges = true;
			m_PolygonMode = GL_FILL;
			m_Culling = false;
			break;
		case 8:
			m_Lighting = false;
			m_SuperimposeOnlyControlMesh = true;
			m_SuperimposeEdges = true;
			m_ThickerControlEdges = false;
			m_PolygonMode = GL_FILL;
			m_Culling = false;
		}
		CMeshDoc *pDoc = (CMeshDoc *)GetDocument();
		InvalidateRect(NULL,FALSE);
	}
}

// change material
void CMeshView::ChangeMaterial(CString &string,
															 bool update)
{
	float	ambient[]  = {0.0f,0.0f,0.0f,1.0f};
	float	diffuse[]  = {0.0f,0.0f,0.0f,1.0f};
	float	specular[]  = {0.0f,0.0f,0.0f,1.0f};
	float	emission[]  = {0.3f,0.3f,0.3f,1.0f};
	float shininess[] = {0.0f};

	// Change
	if(string == "Silver")
	{
		// Ambient
		ambient[0] = 0.19225f;
		ambient[1] = 0.19225f;
		ambient[2] = 0.19225f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.50754f;
		diffuse[1] = 0.50754f;
		diffuse[2] = 0.50754f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.508273f;
		specular[1] = 0.508273f;
		specular[2] = 0.508273f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 51.2f;
	}

	else

	// Change
	if(string == "Gold")
	{
		// Ambient
		ambient[0] = 0.24725f;
		ambient[1] = 0.1995f;
		ambient[2] = 0.0745f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.75164f;
		diffuse[1] = 0.60648f;
		diffuse[2] = 0.22648f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.928281f;
		specular[1] = 0.855802f;
		specular[2] = 0.666065f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 51.2f;
	}

	else

	// Change
	if(string == "Jade")
	{
		// Ambient
		ambient[0] = 0.135f;
		ambient[1] = 0.2225f;
		ambient[2] = 0.1575f;
		ambient[3] = 0.95f;
		// Diffuse
		diffuse[0] = 0.54f;
		diffuse[1] = 0.89f;
		diffuse[2] = 0.63f;
		diffuse[3] = 0.95f;
		// Specular
		specular[0] = 0.316228f;
		specular[1] = 0.316228f;
		specular[2] = 0.316228f;
		specular[3] = 0.95f;
		// Shininess
		shininess[0] = 12.8f;
	}

	else

	// Change
	if(string == "Light blue")
	{
		// Ambient
		ambient[0] = 0.0f;
		ambient[1] = 0.5f;
		ambient[2] = 0.75f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.0f;
		diffuse[1] = 0.5f;
		diffuse[2] = 1.0f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.75f;
		specular[1] = 0.75f;
		specular[2] = 0.75f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 64.0f;
	}

	else

	// Change
	if(string == "Emerald")
	{
		// Ambient
		ambient[0] = 0.0215f;
		ambient[1] = 0.1745f;
		ambient[2] = 0.0215f;
		ambient[3] = 0.55f;
		// Diffuse
		diffuse[0] = 0.07568f;
		diffuse[1] = 0.61424f;
		diffuse[2] = 0.07568f;
		diffuse[3] = 0.55f;
		// Specular
		specular[0] = 0.633f;
		specular[1] = 0.727811f;
		specular[2] = 0.633f;
		specular[3] = 0.55f;
		// Shininess
		shininess[0] = 76.8f;
	}

	else

	// Change
	if(string == "Polished silver")
	{
		// Ambient
		ambient[0] = 0.23125f;
		ambient[1] = 0.23125f;
		ambient[2] = 0.23125f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.2775f;
		diffuse[1] = 0.2775f;
		diffuse[2] = 0.2775f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.773911f;
		specular[1] = 0.773911f;
		specular[2] = 0.773911f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 89.6f;
	}

	else

	// Change
	if(string == "Chrome")
	{
		// Ambient
		ambient[0] = 0.25f;
		ambient[1] = 0.25f;
		ambient[2] = 0.25f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.4f;
		diffuse[1] = 0.4f;
		diffuse[2] = 0.4f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.774597f;
		specular[1] = 0.774597f;
		specular[2] = 0.774597f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 76.8f;
	}

	else

	// Change
	if(string == "Copper")
	{
		// Ambient
		ambient[0] = 0.19125f;
		ambient[1] = 0.0735f;
		ambient[2] = 0.0225f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.7038f;
		diffuse[1] = 0.27048f;
		diffuse[2] = 0.0828f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.256777f;
		specular[1] = 0.137622f;
		specular[2] = 0.086014f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 12.8f;
	}

	else

	// Change
	if(string == "Polished gold")
	{
		// Ambient
		ambient[0] = 0.24725f;
		ambient[1] = 0.2245f;
		ambient[2] = 0.0645f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.34615f;
		diffuse[1] = 0.3143f;
		diffuse[2] = 0.0903f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.797357f;
		specular[1] = 0.723991f;
		specular[2] = 0.208006f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 83.2f;
	}

	else

	// Change
	if(string == "Pewter")
	{
		// Ambient
		ambient[0] = 0.105882f;
		ambient[1] = 0.058824f;
		ambient[2] = 0.113725f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.427451f;
		diffuse[1] = 0.470588f;
		diffuse[2] = 0.541176f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.333333f;
		specular[1] = 0.333333f;
		specular[2] = 0.521569f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 9.84615f;
	}

	else

	// Change
	if(string == "Obsidian")
	{
		// Ambient
		ambient[0] = 0.05375f;
		ambient[1] = 0.05f;
		ambient[2] = 0.06625f;
		ambient[3] = 0.82f;
		// Diffuse
		diffuse[0] = 0.18275f;
		diffuse[1] = 0.17f;
		diffuse[2] = 0.22525f;
		diffuse[3] = 0.82f;
		// Specular
		specular[0] = 0.332741f;
		specular[1] = 0.328634f;
		specular[2] = 0.346435f;
		specular[3] = 0.82f;
		// Shininess
		shininess[0] = 38.4f;
	}

	else

	// Change
	if(string == "Black plastic")
	{
		// Ambient
		ambient[0] = 0.0f;
		ambient[1] = 0.0f;
		ambient[2] = 0.0f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.01f;
		diffuse[1] = 0.01f;
		diffuse[2] = 0.01f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.5f;
		specular[1] = 0.5f;
		specular[2] = 0.5f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 32.0f;
	}

	else

	// Change
	if(string == "Polished bronze")
	{
		// Ambient
		ambient[0] = 0.25f;
		ambient[1] = 0.148f;
		ambient[2] = 0.006475f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.4f;
		diffuse[1] = 0.2368f;
		diffuse[2] = 0.1036f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.774597f;
		specular[1] = 0.458561f;
		specular[2] = 0.200621f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 76.8f;
	}

	
	else

	// Change
	if(string == "Polished copper")
	{
		// Ambient
		ambient[0] = 0.2295f;
		ambient[1] = 0.08825f;
		ambient[2] = 0.0275f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.5508f;
		diffuse[1] = 0.2118f;
		diffuse[2] = 0.066f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.580594f;
		specular[1] = 0.223257f;
		specular[2] = 0.0695701f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 51.2f;
	}

	else

	// Change
	if(string == "Pearl")
	{
		// Ambient
		ambient[0] = 0.25f;
		ambient[1] = 0.20725f;
		ambient[2] = 0.20725f;
		ambient[3] = 0.922f;
		// Diffuse
		diffuse[0] = 1.0f;
		diffuse[1] = 0.829f;
		diffuse[2] = 0.829f;
		diffuse[3] = 0.922f;
		// Specular
		specular[0] = 0.296648f;
		specular[1] = 0.296648f;
		specular[2] = 0.296648f;
		specular[3] = 0.922f;
		// Shininess
		shininess[0] = 11.264f;
	}

	else

	// Change
	if(string == "Ruby")
	{
		// Ambient
		ambient[0] = 0.1745f;
		ambient[1] = 0.01175f;
		ambient[2] = 0.01175f;
		ambient[3] = 0.55f;
		// Diffuse
		diffuse[0] = 0.61424f;
		diffuse[1] = 0.04136f;
		diffuse[2] = 0.04136f;
		diffuse[3] = 0.55f;
		// Specular
		specular[0] = 0.727811f;
		specular[1] = 0.626959f;
		specular[2] = 0.626959f;
		specular[3] = 0.55f;
		// Shininess
		shininess[0] = 76.8f;
	}

	else

	// Change
	if(string == "Turquoise")
	{
		// Ambient
		ambient[0] = 0.1f;
		ambient[1] = 0.18725f;
		ambient[2] = 0.1745f;
		ambient[3] = 0.8f;
		// Diffuse
		diffuse[0] = 0.396f;
		diffuse[1] = 0.74151f;
		diffuse[2] = 0.69102f;
		diffuse[3] = 0.8f;
		// Specular
		specular[0] = 0.297254f;
		specular[1] = 0.30829f;
		specular[2] = 0.306678f;
		specular[3] = 0.8f;
		// Shininess
		shininess[0] = 12.8f;
	}

	else

	// Change
	if(string == "Brass")
	{
		// Ambient
		ambient[0] = 0.329412f;
		ambient[1] = 0.223529f;
		ambient[2] = 0.027451f;
		ambient[3] = 1.0f;
		// Diffuse
		diffuse[0] = 0.780392f;
		diffuse[1] = 0.268627f;
		diffuse[2] = 0.113725f;
		diffuse[3] = 1.0f;
		// Specular
		specular[0] = 0.992157f;
		specular[1] = 0.741176f;
		specular[2] = 0.807843f;
		specular[3] = 1.0f;
		// Shininess
		shininess[0] = 27.8974f;
	}

	// apply
	glMaterialfv( GL_FRONT, GL_AMBIENT,   ambient);
	glMaterialfv( GL_FRONT, GL_DIFFUSE,   diffuse);
	glMaterialfv( GL_FRONT, GL_SPECULAR,  specular);
	glMaterialfv( GL_FRONT, GL_SHININESS, shininess);
  glMaterialfv( GL_FRONT, GL_EMISSION,  emission);

	if(update)
		InvalidateRect(NULL,FALSE);
}

void CMeshView::OnMaterialPearl()
{
	ChangeMaterial(CString("Pearl"),false);
	InvalidateRect(NULL,FALSE);
}

void CMeshView::OnMaterialBrass()
{
	ChangeMaterial(CString("Brass"),false);
	InvalidateRect(NULL,FALSE);
}

void CMeshView::OnMaterialBlackplastic()
{
	ChangeMaterial(CString("Black plastic"),false);
	InvalidateRect(NULL,FALSE);
}

void CMeshView::OnMaterialGold()
{
	ChangeMaterial(CString("Gold"),false);
	InvalidateRect(NULL,FALSE);
}

void CMeshView::OnMaterialSilver()
{
	ChangeMaterial(CString("Silver"),false);
	InvalidateRect(NULL,FALSE);
}

void CMeshView::OnMaterialJade()
{
	ChangeMaterial(CString("Jade"),false);
	InvalidateRect(NULL,FALSE);
}

void CMeshView::OnMaterialRuby()
{
	ChangeMaterial(CString("Ruby"),false);
	InvalidateRect(NULL,FALSE);
}


void CMeshView::OnFileDumptops()
{
  static char BASED_CODE filter[] = "EPS Files (*.eps)|*.eps";
  CFileDialog SaveDlg(FALSE,"*.eps","mesh.eps",
		OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,filter);
  if(SaveDlg.DoModal() == IDOK)
  {
    CString string = SaveDlg.GetPathName();
    char *pFilename = string.GetBuffer(MAX_PATH);

		// drawing
		glPushMatrix();

		// setup viewpoint from current arcball
		m_Arcball.glDraw();

		// get camera
		GLint viewport[4];
		GLdouble modelMatrix[16];
		GLdouble projMatrix[16];
		glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
		glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
		glGetIntegerv(GL_VIEWPORT,viewport);

		glPopMatrix();

		CMeshDoc* pDoc = GetDocument();
		ASSERT_VALID(pDoc);
		Dumper_eps<Enriched_polyhedron<Enriched_kernel,Enriched_items>,Enriched_kernel> dumper;
		dumper.dump(pDoc->m_pMesh,
								pFilename,
								modelMatrix,
								projMatrix,
								viewport,
								false);
    string.ReleaseBuffer();
  }

}
void CMeshView::OnUpdateFileDumptops(CCmdUI *pCmdUI)
{
}

// bounding box
void CMeshView::OnBoundixboxShow()
{
	m_DrawBoundingBox = !m_DrawBoundingBox;
	InvalidateRect(NULL,FALSE);
}
void CMeshView::OnUpdateBoundixboxShow(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_DrawBoundingBox);
}

// show bounding box
void CMeshView::OnBoundixboxShowwhenmoving()
{
	m_DrawBoundingBoxWhenMoving = !m_DrawBoundingBoxWhenMoving;
	InvalidateRect(NULL,FALSE);
}
void CMeshView::OnUpdateBoundixboxShowwhenmoving(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_DrawBoundingBoxWhenMoving);
}

// draw Voronoi edges
void CMeshView::OnSpecialDrawvoronoiedges()
{
	m_DrawVoronoiEdges = !m_DrawVoronoiEdges;
	InvalidateRect(NULL,FALSE);
}
void CMeshView::OnUpdateSpecialDrawvoronoiedges(CCmdUI *pCmdUI)
{
	pCmdUI->SetCheck(m_DrawVoronoiEdges);
}
