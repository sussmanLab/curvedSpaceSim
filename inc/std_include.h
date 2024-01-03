#ifndef STDINCLUDE
#define STDINCLUDE

#include <cmath>
#include <algorithm>
#include <memory>
#include <ctype.h>
#include <random>
#include <stdio.h>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <stdexcept>
#include <cassert>
#include <cpuid.h>

#include "dataTypes.h" //defines simple, non-CGAL data types that will be used

using namespace std; // bad form
//!A utility function for checking if a file exists
inline bool fileExists(const std::string& name)
    {
    std::ifstream f(name.c_str());
    return f.good();
    }
//!Report somewhere that code needs to be written
static void unwrittenCode(const char *message, const char *file, int line)
    {
    printf("\nCode unwritten (file %s; line %d)\nMessage: %s\n",file,line,message);
    throw std::exception();
    }

//spot-checking of code for debugging
#define DEBUGCODEHELPER printf("\nReached: file %s at line %d\n",__FILE__,__LINE__);
//A macro to say code needs to be written
#define UNWRITTENCODE(message) (unwrittenCode(message,__FILE__,__LINE__))
//a macro to say something is wrong!
#define ERRORERROR(message) (unwrittenCode(message,__FILE__,__LINE__))
#endif
