@echo off

@where /q python
if not errorlevel 1 ( set python=python )

@where /q python2
if not errorlevel 1 ( set python=python2 )

@where /q python2.6
if not errorlevel 1 ( set python=python2.6 )

@where /q python2.7
if not errorlevel 1 ( set python=python2.7 )

:go
%python% ${CMAKE_BINARY_DIR}/pkglist_filter.py %1

@echo on
