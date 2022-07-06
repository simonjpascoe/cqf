// Taken and modified from http://unixwiz.net/techtips/outputdebugstring.html
#include "stdafx.h"
#include "Debugging.h"

#include <sstream>
#include <string>
#include <vector>
#include <Windows.h>

using namespace std;

namespace Common {

void odprintf(string format, ...)
{
    char    buf[4096];
    char *p = buf;
    va_list args;

    ostringstream ss;
    ss << "[CQF]/[" << GetCurrentThreadId() << "]" << format;
    string fmt2 = ss.str();

    va_start(args, format);
    auto n = _vsnprintf_s(buf, sizeof buf - 3, fmt2.c_str(), args); // buf-3 is room for CR/LF/NUL
    va_end(args);

    p += (n < 0) ? sizeof buf - 3 : n;

    while ( p > buf  &&  isspace(p[-1]) )
            *--p = '\0';

    *p++ = '\r';
    *p++ = '\n';
    *p   = '\0';

    OutputDebugStringA(buf);
}

string vector_to_string(vector<double> &doubles)
{   
    ostringstream stm; 
    for (auto d: doubles) {
        stm << d << " ";
    }
    auto result = stm.str();
    return result;
}

}