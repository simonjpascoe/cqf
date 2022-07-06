#pragma once

#include <memory>
#include <vector>
#include <xlw/xlw.h>
#include <xlw/XlfServices.h>

using namespace std;
using namespace xlw;

namespace Common {

class Portfolio 
{
public:
    virtual ~Portfolio()
    {
        _values.clear();
    }

    const vector<pair<string,double>>& Values() const
    {
        return _values;
    }

    void Add(const string instrument, const double quantity)
    {

        _values.push_back(make_pair(instrument, quantity));
    }

private:
    vector<pair<string,double>> _values;

};

}