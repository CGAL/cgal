#pragma once


// Boîte de dialogue CDialogOptions

class CDialogOptions : public CDialog
{
	DECLARE_DYNAMIC(CDialogOptions)

public:
	CDialogOptions(CWnd* pParent = NULL);   // constructeur standard
	virtual ~CDialogOptions();

// Données de boîte de dialogue
	enum { IDD = IDD_DIALOG_OPTIONS };

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // Prise en charge DDX/DDV

	DECLARE_MESSAGE_MAP()
public:
	double m_criterion_uniform_size;
	unsigned int m_max_nb_vertices;
	unsigned int m_refresh_each;
	unsigned int m_init_nb_vertices;
};
