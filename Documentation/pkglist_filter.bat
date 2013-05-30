@echo off

@where python
if not errorlevel 1 ( set python=python )

@where python2
if not errorlevel 1 ( set python=python2 )

@where python2.6
if not errorlevel 1 ( set python=python2.6 )

@where python2.7
if not errorlevel 1 ( set python=python2.7 )

@echo on

:go
%python% ${CMAKE_BINARY_DIR}/pkglist_filter.py %1
