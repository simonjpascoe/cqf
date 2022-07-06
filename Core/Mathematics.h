#pragma once
#include <math.h>
#include <vector>
#include <functional>

#include "Matrix.h"

using namespace std;

namespace Common {
    const double M_SQRT2PI = 2.50662827463100050242;
    const double M_1_SQRTPI = M_2_SQRTPI/2.0;

    double interpolate(double xp, const vector<double>& x, const vector<double>& y);

    vector<double> grad(function<double(vector<double>)> f, vector<double> x);
    vector<double> grad(function<double(vector<double>)> f, vector<double> x, double fx);

    pair<double, vector<double>> minimize_bfgs(
        const vector<double>& initialGuess,
        const function<double(vector<double>)>& objective);

    double StandardCumulativeNormal(double u);
    double InverseStandardCumulativeNormal(double p);
}