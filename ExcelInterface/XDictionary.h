#pragma once

#include <memory>
#include <map>
#include <xlw/xlw.h>
#include <xlw/XlfServices.h>

using namespace std;
using namespace xlw;

namespace Common {

typedef CellValue XValue;

class XDictionary 
{
public:
    virtual ~XDictionary()
    {
        m_cache.clear();
    }

    const CellMatrix& LookupAsMatrix(string identity) const
    {
        auto it = m_cache.find(identity);
        if (it!=m_cache.end())
        {
            return it->second;
        }
        else
        {
            throw invalid_argument("Identity `" + identity + "`not found");
        }
    }

    const XValue& Lookup(string identity) const
    {
        return LookupAsMatrix(identity)(0,0);
    }

    void Add(const string key, const CellMatrix& item)
    {
        m_cache.insert(make_pair(key, item));
    }

private:
    map<string, CellMatrix> m_cache;

};
}