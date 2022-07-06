#include <xlw/xlw.h>
#include <xlw/XlfServices.h>

#include "ExcelHelpers.h"
#include "XDictionary.h"
#include "KeyStore.h"
#include "Matrix.h"
#include "ObjectCache.h"
#include "ObjectOwner.h"
#include "HeathJarrowMortonPricer.h"
#include "Instruments.h"

#include "ZeroCouponBond.h"
#include "Cap.h"

using namespace xlw;
using namespace std;
using namespace Common;
using namespace HJM;
using namespace Instruments;

extern "C" {
    LPXLFOPER EXCEL_EXPORT xlCalibrateFactors(LPXLFOPER dictionary)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;

        XlfOper xDictionary(dictionary);
        auto dict = ObjectCache<XDictionary>::Instance().Lookup(xDictionary.AsString());

        auto xMats = dict.Lookup("tenors").StringValue();
        auto xSource = dict.Lookup("timeseries").StringValue();
        auto xCoverage = dict.Lookup("coverage").NumericValue();

        auto cMaturities = ObjectCache<Matrix>::Instance().Lookup(xMats);
        auto cTimeSeriesData = ObjectCache<Matrix>::Instance().Lookup(xSource);
        HeathJarrowMortonPricer hjm;
        auto calibrationResult = hjm.CalibrateFactors(cMaturities, cTimeSeriesData, xCoverage);

        XDictionary output;
        vector<pair<string, freeF>> children;
        auto coefficients = ObjectCache<Matrix>::Instance().Store(calibrationResult.Coefficients);
        children.push_back(make_pair(coefficients, ObjectCache<Matrix>::Instance().freeF()));
        output.Add("coefficients", CellMatrix(coefficients));

        auto eigenvalues = ObjectCache<Matrix>::Instance().Store(Matrix::FromVector(calibrationResult.AnalysisResult.Eigenvalues));
        children.push_back(make_pair(eigenvalues, ObjectCache<Matrix>::Instance().freeF()));
        output.Add("eigenvalues", CellMatrix(eigenvalues));

        auto rows = calibrationResult.AnalysisResult.Eigenvectors[0].Rows();
        auto cols = calibrationResult.AnalysisResult.Eigenvectors.size();
        Matrix M(rows,cols);
        for (int i=0; i<rows; i++)
        {
            for (int j=0; j<cols; j++)
            {
                M(i,j) = calibrationResult.AnalysisResult.Eigenvectors[j](i,0);
            }
        }

        auto eigenvectors = ObjectCache<Matrix>::Instance().Store(M);
        children.push_back(make_pair(eigenvectors, ObjectCache<Matrix>::Instance().freeF()));
        output.Add("eigenvectors", CellMatrix(eigenvectors));

        output.Add("coverage", CellMatrix(calibrationResult.AnalysisResult.Coverage));

        // ownership management
        string key = ObjectCache<XDictionary>::Instance().Store(output);
        auto owner = XlfServices.Information.GetCallingCell().AsRef();
        auto old = ObjectOwner::Instance().Mark(owner, make_pair(key, children));
        if (old.first!="") {
            ObjectCache<XDictionary>::Instance().Remove(old.first);
            for (auto cached : old.second)
            {
                cached.second(cached.first);
            }
        } 

        return XlfOper(key);
        EXCEL_END;
    }

    LPXLFOPER EXCEL_EXPORT xlSimulate(
        LPXLFOPER dictionary
        )
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;

        XlfOper xDictionary(dictionary);
        auto dict = ObjectCache<XDictionary>::Instance().Lookup(xDictionary.AsString());

        auto xCoefficients = dict.Lookup("coefficients").StringValue();
        auto xTenors = dict.Lookup("tenors").StringValue();
        auto xSeeds = dict.Lookup("seed").StringValue();
        auto xTimestep = dict.Lookup("timestep").NumericValue();
        auto xSimulationCount = (unsigned long)dict.Lookup("simulations").NumericValue();

        auto cCoefficients = ObjectCache<Matrix>::Instance().Lookup(xCoefficients);
        auto cTenors = ObjectCache<Matrix>::Instance().Lookup(xTenors);
        auto cSeeds = ObjectCache<Matrix>::Instance().Lookup(xSeeds);
        
        double maturity = 0;
        vector<HeathJarrowMortonPriceable*> instruments;
        auto keyStore = dict.Lookup("instruments").StringValue();
        auto instrumentKeys = ObjectCache<KeyStore>::Instance().Lookup(keyStore);
        for (auto k : instrumentKeys.Values()) {
            auto instrument = ObjectCache<shared_ptr<Instruments::Instrument>>::Instance().Lookup(k);
            maturity = max(maturity, instrument->Maturity());
            auto two = dynamic_cast<HeathJarrowMortonPriceable*>(instrument.get());
            instruments.push_back(two);
        }

        HeathJarrowMortonPricer hjm;
        auto result = hjm.Simulate(instruments, cCoefficients, cTenors, cSeeds, maturity, xTimestep, xSimulationCount);
        
        XDictionary output;
        vector<pair<string, freeF>> children;
        auto fairvalues = ObjectCache<Matrix>::Instance().Store(result.FairValues);
        children.push_back(make_pair(fairvalues, ObjectCache<Matrix>::Instance().freeF()));
        output.Add("fairvalues", CellMatrix(fairvalues));

        auto convergence = ObjectCache<Matrix>::Instance().Store(result.Convergence);
        children.push_back(make_pair(convergence, ObjectCache<Matrix>::Instance().freeF()));
        output.Add("convergence", CellMatrix(convergence));

        // ownership management
        string key = ObjectCache<XDictionary>::Instance().Store(output);
        auto owner = XlfServices.Information.GetCallingCell().AsRef();
        auto old = ObjectOwner::Instance().Mark(owner, make_pair(key, children));
        if (old.first!="") {
            ObjectCache<XDictionary>::Instance().Remove(old.first);
            for (auto cached : old.second)
            {
                cached.second(cached.first);
            }
        } 
    
        return XlfOper(key);
        EXCEL_END;
    }    

    
}

namespace {
    XLRegistration::Arg argsXlXDictionary[] = {
        { "dictionary", "Dictionary", "XLF_OPER" }
    };

    XLRegistration::XLFunctionRegistrationHelper register__calibrate_factors(
        "xlCalibrateFactors", "hjm_calibrate_factors", "Perform PCA on the input timeseries data",
        "CQF", argsXlXDictionary, 1);

    XLRegistration::XLFunctionRegistrationHelper register__simulate(
        "xlSimulate", "hjm_simulate", "Simulation and price products",
        "CQF", argsXlXDictionary, 1);



}
