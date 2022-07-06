#include <xlw/xlw.h>
#include <xlw/XlfServices.h>
#include <xlw/XlOpenClose.h>

#include "ExcelHelpers.h"
#include "XDictionary.h"
#include "KeyStore.h"
#include "Matrix.h"
#include "Mathematics.h"

using namespace xlw;
using namespace std;
using namespace Common;

extern "C" {
    LPXLFOPER EXCEL_EXPORT xlDependency(LPXLFOPER value, LPXLFOPER dummy)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;

        return XlfOper(value);
        EXCEL_END;
    }

    LPXLFOPER EXCEL_EXPORT xlXDictionaryCreate(LPXLFOPER source)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;

        XlfOper xSource(source);

        if (xSource.columns() != 2) { THROW_XLW("An XDictionary must have two columns"); }
        
        XDictionary dict;
        for (auto i = 0; i<xSource.rows(); i++)
        {
            dict.Add(xSource(i,0).AsString(), xSource(i,1).AsCellMatrix());
        }

        auto key = cache_item(dict);
        
        return XlfOper(key);
        EXCEL_END;
    }

    LPXLFOPER EXCEL_EXPORT xlXDictionaryLookup(LPXLFOPER source, LPXLFOPER key)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;

        XlfOper xSource(source);
        XlfOper xKey(key);

        auto dict = ObjectCache<XDictionary>::Instance().Lookup(xSource.AsString());
        return XlfOper(dict.LookupAsMatrix(xKey.AsString()));        
        
        EXCEL_END;
    }

    LPXLFOPER EXCEL_EXPORT xlKeyStoreCreate(LPXLFOPER source)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;

        XlfOper xSource(source);

        auto cells = xSource.AsCellMatrix();
        if (cells.ColumnsInStructure() > 1 && cells.RowsInStructure() >1)
            { THROW_XLW("At least one dimension for key store should 1"); }

        KeyStore store;
        for (size_t i = 0; i<cells.RowsInStructure(); i++) 
        {
            for (size_t j = 0; j<cells.ColumnsInStructure(); j++)
            {
                store.Add(cells(i,j).StringValue());
            }
        }

        auto key = cache_item(store);
        
        return XlfOper(key);
        EXCEL_END;
    }

    LPXLFOPER EXCEL_EXPORT xlInterpolate(LPXLFOPER dictionary)
    {
        EXCEL_BEGIN;
        NO_FUNC_WIZARD;

        XlfOper xDictionary(dictionary);
        auto dict = ObjectCache<XDictionary>::Instance().Lookup(xDictionary.AsString());

        auto xX = dict.Lookup("x").StringValue();
        auto xY = dict.Lookup("y").StringValue();
        auto xP = dict.Lookup("xp").NumericValue();

        auto x = ObjectCache<Matrix>::Instance().Lookup(xX);
        auto y = ObjectCache<Matrix>::Instance().Lookup(xY);

        auto result = interpolate(xP, x.ToVector(0), y.ToVector(0));

        return XlfOper(result);
        EXCEL_END;
    }
}

namespace {
    XLRegistration::Arg argsXlDependency[] = {
        { "value", "Value to return", "XLF_OPER" },
        { "dummy", "Dummy link cell", "XLF_OPER" }
    };
    
    XLRegistration::XLFunctionRegistrationHelper register__dependency(
        "xlDependency", "dependency", "Add a cell dependency into the excel calc tree",
        "CQF", argsXlDependency, 2);

    XLRegistration::Arg argsXlXDictionaryCreate[] = {
        { "content", "Row based dictionary values", "XLF_OPER" }
    };
    
    XLRegistration::XLFunctionRegistrationHelper register__dictionary_create(
        "xlXDictionaryCreate", "dictionary_create", "Create a dictionary of values",
        "CQF", argsXlXDictionaryCreate, 1);

    XLRegistration::Arg argsXlXDictionaryLookup[] = {
        { "dictionary", "Source dictionary", "XLF_OPER" },
        { "key", "Case sensitive key to lookup", "XLF_OPER" }
    };
    
    XLRegistration::XLFunctionRegistrationHelper register__dictionary_lookup(
        "xlXDictionaryLookup", "dictionary_lookup", "Lookup a value from a dictionary",
        "CQF", argsXlXDictionaryLookup, 2);

    XLRegistration::Arg argsXlKeyStoreCreate[] = {
        { "content", "Vector of keys", "XLF_OPER" }
    };

    XLRegistration::XLFunctionRegistrationHelper register__keystore_create(
        "xlKeyStoreCreate", "keystore_create", "Create a list of keys",
        "CQF", argsXlKeyStoreCreate, 1);
    
    XLRegistration::Arg argsXlInterpolate[] = {
        { "dictionary", "dictionary", "XLF_OPER" }
    };

    XLRegistration::XLFunctionRegistrationHelper register__interpolate(
        "xlInterpolate", "interpolate", "Interpolate vectors",
        "CQF", argsXlInterpolate, 1);

}