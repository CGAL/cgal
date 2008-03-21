#pragma once

#include "resource.h"


// Boîte de dialogue CDialogOptions_3

class CDialogOptions_3 : public CDialog
{
	DECLARE_DYNAMIC(CDialogOptions_3)

public:
	CDialogOptions_3(CWnd* pParent = NULL);   // constructeur standard
	virtual ~CDialogOptions_3();

// Données de boîte de dialogue
	enum { IDD = IDD_DIALOG_OPTIONS_3 };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // Prise en charge DDX/DDV

	DECLARE_MESSAGE_MAP()

public:
	double m_point_size; // OpenGL point size
};
