#pragma once
#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include <vector>

// Only include debugging statements if we are in debug build.
#ifdef _DEBUG
#define dbg(format, ...) odprintf(format, __VA_ARGS__)
#else
#define dbg(format, ...)
#endif

using namespace std;

namespace Common {
	void odprintf(string format, ...);
	string vector_to_string(vector<double> &doubles);
}