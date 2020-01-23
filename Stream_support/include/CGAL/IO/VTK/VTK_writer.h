// Copyright (c) 2018
// GeometryFactory (France) All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb

#ifndef CGAL_IO_VTK_WRITER_H
#define CGAL_IO_VTK_WRITER_H

namespace CGAL {
namespace Stream_support {
namespace internal {

template <class FT>
void write_vector(std::ostream& os,
                  const std::vector<FT>& vect)
{
  const char* buffer = reinterpret_cast<const char*>(&(vect[0]));
  std::size_t size = vect.size()*sizeof(FT);

  os.write(reinterpret_cast<const char *>(&size), sizeof(std::size_t)); // number of bytes encoded
  os.write(buffer, vect.size()*sizeof(FT)); // encoded data
}

class ErrorObserverVtk
  : public vtkCommand
{
public:
  ErrorObserverVtk()
    : Error(false), Warning(false), ErrorMessage(""), WarningMessage("")
  { }

  static ErrorObserverVtk *New() { return new ErrorObserverVtk; }

  bool GetError() const          { return this->Error; }
  bool GetWarning() const        { return this->Warning; }
  std::string GetErrorMessage()   { return ErrorMessage; }
  std::string GetWarningMessage() { return WarningMessage; }

  void Clear()
  {
    this->Error = false;
    this->Warning = false;
    this->ErrorMessage = "";
    this->WarningMessage = "";
  }

  virtual void Execute(vtkObject *vtkNotUsed(caller),
                       unsigned long event,
                       void *calldata)
  {
    switch (event)
    {
      case vtkCommand::ErrorEvent:
        ErrorMessage = static_cast<char *>(calldata);
        this->Error = true;
        break;
      case vtkCommand::WarningEvent:
        WarningMessage = static_cast<char *>(calldata);
        this->Warning = true;
        break;
    }
  }

private:
  bool Error;
  bool Warning;
  std::string ErrorMessage;
  std::string WarningMessage;
};

template <class vtkReader>
vtkSmartPointer<vtkReader>
read_vtk_file(const std::string& input_filename,
              vtkSmartPointer<CGAL::ErrorObserverVtk> errorObserver)
{
  vtkSmartPointer<vtkReader> reader = vtkSmartPointer<vtkReader>::New();
  reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
  reader->AddObserver(vtkCommand::WarningEvent, errorObserver);
  reader->SetFileName(input_filename.data());
  reader->Update();
  return reader;
}

} // namespace internal
} // namespace Stream_support
} // namespace CGAL

#endif // CGAL_IO_VTK_WRITER_H
