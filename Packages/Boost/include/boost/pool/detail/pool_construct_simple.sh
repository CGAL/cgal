#!/bin/sh
#
# Copyright (C) 2001 Stephen Cleary
#
# This file can be redistributed and/or modified under the terms found
#  in "copyright.html"
# This software and its documentation is provided "as is" without express or
#  implied warranty, and with no claim as to its suitability for any purpose.
#
# See http://www.boost.org for updates, documentation, and revision history.
#

m4 -P -E -DNumberOfArguments=$1 pool_construct_simple.m4 > pool_construct_simple.inc
