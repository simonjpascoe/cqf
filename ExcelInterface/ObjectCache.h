#pragma once
#include <functional>
#include <map>
#include "Singleton.h"

using namespace std;

template<typename T>
class ObjectCache : public Singleton<ObjectCache<T>>
{
public:
    const T& Lookup(string identity) const
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

    string Store(const T& item, string overrideName = "")
    {
        // TODO: correct any threading issue here
        m_key++;
        char ident[199];
        sprintf_s(ident, 199, "@@%d:%s",m_key, overrideName == "" ? typeid(T).name() : overrideName.c_str());
        auto p = make_pair(ident, item);
        m_cache.insert(p);
        
        return ident;
    }

    void Remove(string identity)
    {
        m_cache.erase(identity);
    }

    function<void(string)> freeF()
    {
        auto f = [this] (string identity) { this->Remove(identity);};
        return f;
    }

private:
    int m_key;
    map<string, T> m_cache;

};