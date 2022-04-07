/**************************************************************************
// Copyright (c) 2004  Max-Planck-Institut Saarbruecken (Germany)
// All rights reserved.
//
// This file is part of BenchmarkParser
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s) : Lutz Kettner
**************************************************************************/

#include <benchmark_visitor.h>
#include <benchmark_parser.h>
#include <iostream>
#include <sstream>
#include <cstdlib>

void Benchmark_visitor::error_handler( std::string s) {
    std::ostringstream message;
    message << "ERROR: Benchmark_parser: " << benchmark_filename << " line "
            << benchmark_linenumber << ": " << s;
    if ( ! m_error) {
        m_error = true;
        m_error_message = message.str();
    }
    if ( m_mode & PRINT_ERRORS) {
        std::cerr <<  message.str() << std::endl;
    }
    switch ( m_mode & 0x06) {
    case THROW_ERRORS:
        throw Benchmark_parser_exception( message.str());
        break;
    case EXIT_ERRORS:
        exit( 1);
        break;
    default:
        break;
    }
}

void Benchmark_visitor::parse_error( std::string s) {
    //error_handler( std::string("parse error: ") + s);
    error_handler( s);
}

void Benchmark_visitor::unknown_token( std::string s) {
    error_handler( std::string("unknown token '") + s
                   + std::string("' found."));
}

void Benchmark_visitor::token_not_handled( std::string s) {
    error_handler( std::string("known token '") + s
                   + std::string("' found but not handled by visitor."));
}

void Benchmark_visitor::reset_header() {
    m_major = -1;
    m_minor = -1;
    m_format_options = std::string("");
    m_benchmark_name = std::string("");
}

void Benchmark_visitor::accept_file_format(
    std::string s, int major, int minor, std::string options)
{
    reset_header();
    if ( s != BENCHMARK_FORMAT_NAME) {
        error_handler( std::string( "wrong file format '") + s +
                       std::string("', must be '") + BENCHMARK_FORMAT_NAME
                       + std::string("'."));
        return;
    }
    if ( major > BENCHMARK_FORMAT_MAJOR) {
        std::ostringstream message;
        message << "version major " << major << " is larger than "
                << BENCHMARK_FORMAT_MAJOR << " what the parser understands.";
        error_handler( message.str());
        return;
    }
    m_major = major;
    m_minor = minor;
    // Options are comma separated. We add one left and one right to
    // have simple contains-an-option ",OPT," queries.
    m_format_options = std::string(",") + options + std::string(",");
}

void Benchmark_visitor::accept_benchmark_name( std::string s) {
    m_benchmark_name = s;
}

void Benchmark_visitor::accept_classification( std::string problem,
                                               std::string geom,
                                               std::string clas,
                                               std::string family,
                                               std::string instance,
                                               std::string release) {
    std::cout<<"inclassifi visitor" <<std::endl;
    error_handler ("classification not supported");
}
