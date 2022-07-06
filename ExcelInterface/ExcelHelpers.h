#pragma once

#include <map>
#include <string>

#include <xlw/xlw.h>
#include <xlw/XlfServices.h>

#include "ObjectCache.h"
#include "ObjectOwner.h"
#include "XDictionary.h"

using namespace std;
using namespace xlw;
using namespace Common;

#define NO_FUNC_WIZARD if (XlfExcel::Instance().IsCalledByFuncWiz()) return XlfOper("No execution in wizard.");

template<typename T>
inline string cache_item(const T& item, string overrideName = "") { 
    string key = ObjectCache<T>::Instance().Store(item, overrideName);
    auto owner = XlfServices.Information.GetCallingCell().AsRef();
    auto old = ObjectOwner::Instance().Mark(owner, key);
    if (old.first!="") {
        ObjectCache<T>::Instance().Remove(old.first);
        for (auto cached : old.second)
        {
            cached.second(cached.first);
        }
    } 

    return key; 
}

inline string cache_dictionary(const XDictionary& output, const vector<pair<string, freeF>>& children, string overrideName = "") { 
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

    return key; 
}


template<typename T>
inline string cache_child(const T&item, vector<pair<string, freeF>>& children) { 
    auto key = ObjectCache<T>::Instance().Store(item);
    children.push_back(make_pair(key, ObjectCache<T>::Instance().freeF()));
    return key;
}