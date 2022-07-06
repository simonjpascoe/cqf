#pragma once

#include <functional>
#include <memory>
#include "Matrix.h"
#include "HeathJarrowMortonPriceable.h"

using namespace Common;

namespace HJM {

class CalibrationResult
{
public:
    Matrix Coefficients;
    PCAResult AnalysisResult;
};

class PricingResult
{
public:
    Matrix FairValues;
    Matrix Convergence;
};

class HeathJarrowMortonPricer 
{
public:
    HeathJarrowMortonPricer(void) { }
    virtual ~HeathJarrowMortonPricer(void) {}

    CalibrationResult CalibrateFactors(
        const Matrix& maturities, 
        const Matrix& timeSeriesData, 
        double coverageFactor) const;
    PricingResult Simulate(
        const vector<HeathJarrowMortonPriceable*>& products,
        const Matrix& coefficients, 
        const Matrix& tenors, 
        const Matrix& seeds, 
        double maturity, 
        double dt,
        long simulationCount) const;

private:
    double calculateCovariance(const Matrix& primary, const Matrix& secondary) const;
    double average(const Matrix& vector) const;
    function<vector<double>(double)> createBasisFunctions(double y1, double y2, double y3, double y4) const;
    double integrate(const function<double(double)>& f, double lower, double upper) const;
};

}