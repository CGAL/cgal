// GyrovizDoc_4.cpp : implementation of the CGyrovizDoc_4 class
//

#include "stdafx.h"

#include "GyrovizDoc_4.h"
#include "MainFrm.h"
#include "DialogOptions_4.h"

using namespace cimg_library;

#include "Gyroviz_segmentation2.h"

#include <string>
#include <sstream>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CGyrovizDoc_4

IMPLEMENT_DYNCREATE(CGyrovizDoc_4, CDocument)

BEGIN_MESSAGE_MAP(CGyrovizDoc_4, CDocument)
END_MESSAGE_MAP()


// CGyrovizDoc_4 construction/destruction

CGyrovizDoc_4::CGyrovizDoc_4()
{
}

CGyrovizDoc_4::~CGyrovizDoc_4()
{
}

BOOL CGyrovizDoc_4::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;
	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)
	return TRUE;
}

// CGyrovizDoc_4 serialization

void CGyrovizDoc_4::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}


// CGyrovizDoc_4 diagnostics

#ifdef _DEBUG
void CGyrovizDoc_4::AssertValid() const
{
	CDocument::AssertValid();
}

void CGyrovizDoc_4::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

// CGyrovizDoc_4 commands


// Update the number of vertices and faces in the status bar
void CGyrovizDoc_4::update_status()
{   
	CWinApp *pApp = AfxGetApp();
	if(pApp->m_pMainWnd != NULL) 
	{ 
		CStatusBar* pStatus = 
			(CStatusBar*)AfxGetApp()->m_pMainWnd->GetDescendantWindow(
			AFX_IDW_STATUS_BAR);

		if(pStatus != NULL) 
		{
			//CString vertices;
			//vertices.Format("%d vertices",m_gyroviz_dt.number_of_vertices());

			//CString faces;
			//faces.Format("%d faces",m_gyroviz_dt.number_of_faces());

			//// Update status bar
			//pStatus->SetPaneText(1,vertices);
			//pStatus->SetPaneText(2,faces);
			//pStatus->UpdateWindow(); 
		}
	}
}

// User message in status bar
void CGyrovizDoc_4::status_message(char* fmt,...)
{   
	CWinApp *pApp = AfxGetApp();
	if(pApp->m_pMainWnd != NULL) 
	{ 
		char buffer[256];
		CStatusBar* pStatus = 
			(CStatusBar*)AfxGetApp()->m_pMainWnd->GetDescendantWindow(
			AFX_IDW_STATUS_BAR);

		// fill buffer
		va_list argptr;      
		va_start(argptr,fmt);
		vsprintf(buffer,fmt,argptr);
		va_end(argptr);

		if(pStatus != NULL) 
		{
			pStatus->SetPaneText(0,buffer);
			pStatus->UpdateWindow(); 
		}
	}
	return;
}


BOOL CGyrovizDoc_4::OnSaveDocument(LPCTSTR lpszPathName)
{
	// save pslg to a file
	//return m_pslg.save((char *)lpszPathName);
	return CDocument::OnSaveDocument(lpszPathName);
}


bool CGyrovizDoc_4::read_sequence(const std::string& first_image, 
								  const std::string& first_feature_file)
{
	const char keys[] = "1234567890";

	//
	// Compute images mask
	//

	// get span until first number in string out in image_base
	int image_base_length;
	image_base_length = first_image.find_last_of (keys);
	image_base_length = first_image.find_last_not_of (keys, image_base_length) + 1;
	std::string image_base = first_image.substr(0, image_base_length);

	// get extension
	std::string image_extension = first_image.substr(first_image.length()-4);

	// create mask from first_image
	int image_number_of_digits = first_image.length()-image_base.length()-image_extension.length();
	std::string image_mask = image_base + std::string(image_number_of_digits,'?') + image_extension;

	// get path
	int image_path_length;
	image_path_length = first_image.find_last_of ("\\") + 1;
	std::string image_path = first_image.substr(0, image_path_length);

	//
	// Compute feature files (.pnt) mask
	//

	// get span until first number in string out in feature_file_base
	int feature_file_base_length;
	feature_file_base_length = first_feature_file.find_last_of (keys);
	feature_file_base_length = first_feature_file.find_last_not_of (keys, feature_file_base_length) + 1;
	std::string feature_file_base = first_feature_file.substr(0, feature_file_base_length);

	// get extension
	std::string feature_file_extension = first_feature_file.substr(first_feature_file.length()-4);

	// create mask from first_feature_file
	int feature_file_number_of_digits = first_feature_file.length()-feature_file_base.length()-feature_file_extension.length();
	//std::string feature_file_mask = feature_file_base + std::string(feature_file_number_of_digits,'?') + feature_file_extension;

	//
	// Get all image files
	//

	std::vector<std::string> images;

	// search for the first image file
	WIN32_FIND_DATA FindFileData; 
	HANDLE hFind;				  
	DWORD a = 0;				  
	hFind = FindFirstFile(image_mask.c_str(), &FindFileData);
	if (hFind==INVALID_HANDLE_VALUE)
		a = ERROR_NO_MORE_FILES;
	while (a != ERROR_NO_MORE_FILES )
	{
		// Note file names
		if ( strcmp(FindFileData.cFileName,"." )!=0 && 
			strcmp(FindFileData.cFileName,"..")!=0 &&
			(FindFileData.dwFileAttributes != FILE_ATTRIBUTE_DIRECTORY) )
		{
			images.push_back(image_path + std::string(FindFileData.cFileName));
		}

		// Move to next entry
		if (!FindNextFile(hFind, &FindFileData))
			a = GetLastError();
	}
	FindClose(hFind);

	// sort by image index
	sort(images.begin(), images.end());

	//
	// Loop over all image files
	//

	for (std::vector<std::string>::iterator image_it = images.begin(); 
		image_it != images.end(); 
		image_it++)
	{
		string image_name = *image_it;

		//
		// Read and segment image
		//

		CImg <unsigned char> cimg_interm_image;
		CImg<unsigned char> cimg_gray_image;
		CImg<unsigned char> cimg_filt_image;
		CImg<unsigned char> cimg_seg_image;

		unsigned char*  original_image;  
		unsigned char* grayscaled_image;
		unsigned char* filtered_image;
		unsigned char* segmented_image;

		try {
			cimg_interm_image = CImg<unsigned char>(image_name.c_str());
		}
		catch(...) {
			AfxMessageBox("Unable to open image file");
			return false;
		}
		original_image = cimg_image_multiplexer_char(cimg_interm_image);

		cimg_gray_image = to_grayscale(cimg_interm_image);
		grayscaled_image = cimg_image_multiplexer_char(cimg_gray_image);

		cimg_filt_image = gauss3(cimg_gray_image);
		filtered_image = cimg_image_multiplexer_char(cimg_filt_image);

		cimg_seg_image = grad_freiChen(cimg_filt_image,15); // Edge detection on filtered image def : 30
		segmented_image = cimg_image_multiplexer_char(cimg_seg_image);

		//
		// Read corresponding Voodoo 2D feature file (.pnt)
		//

		CDt2 gyroviz_dt;

		int image_index = atoi(image_name.substr(image_base_length, image_number_of_digits).c_str());
		std::ostringstream out;
		out << setfill('0') << std::setw(feature_file_number_of_digits) << image_index;
		std::string feature_index = out.str();
		std::string feature_file_name = feature_file_base + feature_index + feature_file_extension;

		if(!gyroviz_dt.read_pnt(feature_file_name.c_str()))
		{
			AfxMessageBox("Unable to open features file");
			return false;
		}
		
		gyroviz_dt.nw_add_constraints(cimg_seg_image, 5);

		m_list_cdt2.push_back(gyroviz_dt);
		m_vector_triangle3_wc = CDt2::list_cdt2_to_vector_twc(m_list_cdt2);


		//
		// Clean up images
		//

		delete [] original_image;  
		delete [] grayscaled_image;
		delete [] filtered_image;
		delete [] segmented_image;


	}

	// Debug
	int size_m_list_cdt2 = m_list_cdt2.size();
	int size_m_vector_triangle3_wc = m_vector_triangle3_wc.size();

	return true;
}


BOOL CGyrovizDoc_4::OnOpenDocument(LPCTSTR lpszPathName)
{
	//
	// Read bitmap image
	//

	if (!CDocument::OnOpenDocument(lpszPathName))
		return FALSE;

	// get extension
	CString file = lpszPathName;
	CString extension = lpszPathName;
	extension = extension.Right(4);
	extension.MakeLower();

	// set current path
	// path "c:\path\file.wrl" -> c:\path
	CString path = lpszPathName;
	path = path.Left(path.ReverseFind('\\'));
	SetCurrentDirectory(path);

	if(extension != ".bmp")
	{
		AfxMessageBox("File format not supported");
		return FALSE;
	}

	//
	// Read corresponding Voodoo 2D feature file (.pnt)
	//

	// file filters
	static char szFilter[] = "Voodoo 2D Feature Files (*.pnt)|*.pnt|All Files (*.*)|*.*||";

	// create the Open dialog
	CFileDialog dlgOpen(true, "pnt", NULL, 
		OFN_DONTADDTORECENT | OFN_FILEMUSTEXIST, szFilter, AfxGetMainWnd());

	// dialog title
	dlgOpen.m_ofn.lpstrTitle = "Select corresponding Voodoo 2D feature file";

	// show the dialog
	if (dlgOpen.DoModal() == IDOK)
	{
		// get extension
		CString file = dlgOpen.m_ofn.lpstrFile;
		CString extension = dlgOpen.m_ofn.lpstrFile;
		extension = extension.Right(4);
		extension.MakeLower();

		// set current path
		CString path = dlgOpen.m_ofn.lpstrFile;
		path = path.Left(path.ReverseFind('\\'));
		SetCurrentDirectory(path);

		// if .pnt extension
		if(extension.CompareNoCase(".pnt") == 0)
		{
			double init = clock();
			if(!read_sequence(lpszPathName, dlgOpen.m_ofn.lpstrFile))
			{
				AfxMessageBox("Unable to open sequence");
				return FALSE;
			}
			status_message("Constrained Delaunay triangle soup (%lf s)",duration(init));
		}
		else
		{
			AfxMessageBox("File format not supported");
			return FALSE;
		}
	}
	else
	{
		return FALSE;
	}

	update_status();
	UpdateAllViews(NULL);
	return TRUE;
}

double CGyrovizDoc_4::duration(const double time_init)
{
	return (clock() - time_init)/CLOCKS_PER_SEC;
}

