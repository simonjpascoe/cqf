#pragma once

#include <memory>
#include <vector>
#include <xlw/xlw.h>
#include <xlw/XlfServices.h>

using namespace std;
using namespace xlw;

namespace Common {

class KeyStore 
{
public:
    virtual ~KeyStore()
    {
        m_keys.clear();
    }

    const vector<string>& Values() const
    {
        return m_keys;
    }

    void Add(const string key)
    {

        m_keys.push_back(key);
    }

private:
    vector<string> m_keys;

};

}