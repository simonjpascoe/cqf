#pragma once

#include <vector>
#include <string>

using namespace std;

namespace Financials {

template <typename T>
class Quote
{
public: 
    T Instrument;
    double Bid;
    double Offer;
};

template <typename T>
class Quotes
{
public:
    void Add(const Quote<T>& quote)
    {
        _quotes.push_back(quote);
    }

    const vector<Quote<T>>& All() const
    {
        return _quotes;
    }

private:
    vector<Quote<T>> _quotes;
};

class MarketData
{
public:
    double HighVolatility;
    double LowVolatility;
    double InterestRate;
    double Spot;
};

}