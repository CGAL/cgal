#
# QtChoice.cmake needs to be included into CMakeList when a macros compatibility between Qt4 and Qt5 
# necessary.
#
# The purpose of this file is to substitute Qt version notification into Qt macro.
#
# For instance qt4_wrap_ui or qt5_wrap_ui could be replace by just qt_wrap_ui into CMakeList.
#
# The Qt macros supported are the followings with their Qt4 and Qt5 equivalents :
#
#  Name of the macro                        Qt4                              Qt5
#    qt_automoc                         qt4_automoc                    Qt5 doesn't need automoc
#    qt_wrap_ui                         qt4_wrap_ui                    qt5_wrap_ui
#    qt_wrap_cpp                        qt4_wrap_cpp                   qt5_wrap_cpp
#    qt_generate_moc                    qt4_generate_moc               qt5_generate_moc
#    qt_add_dbus_adaptor                qt4_add_dbus_adaptor           qt5_add_dbus_adaptor
#    qt_add_dbus_interfaces             qt4_add_dbus_interfaces        qt5_add_dbus_interfaces
#    qt_add_dbus_interface              qt4_add_dbus_interface         qt5_add_dbus_interface
#    qt_generate_dbus_interface         qt_generate_dbus_interface     qt5_generate_dbus_interface
#    qt_add_resources                   qt_add_resources               qt5_add_resources
#
#
#

macro(qt_automoc)
  if(QT5)
  	message(STATUS "Qt5 doesn't need automoc.")
  elseif(QT4)
   	qt4_automoc(${ARGN})
  else()
        message("No Qt selected for add_resources.")
  endif()
endmacro()

macro(qt_wrap_ui)
  if(QT5)
    qt5_wrap_ui(${ARGN})
  elseif(QT4)
   	qt4_wrap_ui(${ARGN})
  else()
   	message("No Qt selected for wrap_ui.")
  endif()
endmacro()

macro(qt_wrap_cpp)
  if(QT5)
     qt5_wrap_cpp(${ARGN})
  elseif(QT4)
     qt4_wrap_cpp(${ARGN})
  else()
   message("No Qt selected for wrap_cpp.")
  endif()
endmacro()

macro(qt_generate_moc)
  if(QT5)
     qt5_generate_moc(${ARGN})
  elseif(QT4)
     qt4_generate_moc(${ARGN})
  else()
   message("No Qt selected for generate_moc.")
  endif()
endmacro()

macro(qt_add_dbus_adaptor)
  if(QT5)
     qt5_add_dbus_adaptor(${ARGN})
  elseif(QT4)
     qt4_add_dbus_adaptor(${ARGN})
  else()
   message("No Qt selected for add_dbus_adaptor.")
  endif()
endmacro()

macro(qt_add_dbus_interfaces)
  if(QT5)
      qt5_add_dbus_interfaces(${ARGN})
  elseif(QT4)
      qt4_add_dbus_interfaces(${ARGN})
  else()
   message("No Qt selected for add_dbus_interfaces.")
  endif()
endmacro()

macro(qt_add_dbus_interface)
  if(QT5)
  	qt5_add_dbus_interface(${ARGN})
  elseif(QT4)
  	qt4_add_dbus_interface(${ARGN})
  else()
   message("No Qt selected for add_dbus_interface.")
  endif()
endmacro()

macro(qt_generate_dbus_interface)
  if(QT5)
  	qt5_generate_dbus_interface(${ARGN})
  elseif(QT4)
  	qt4_generate_dbus_interface(${ARGN})
  else()
   message("No Qt selected for generate_dbus_interface.")
  endif()
endmacro()

macro(qt_add_resources)
  if(QT5)
  	qt5_add_resources(${ARGN})
  elseif(QT4)
   	qt4_add_resources(${ARGN})
  else()
   message("No Qt selected for add_resources.")
  endif()
endmacro()

