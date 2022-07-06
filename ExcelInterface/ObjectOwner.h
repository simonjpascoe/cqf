#pragma once

#include <functional>

#include <xlw/xlw.h>
#include <xlw/XlfServices.h>
#include "Singleton.h"

using namespace std;
using namespace xlw;

typedef function<void(string)> freeF;

typedef pair<string, vector<pair<string,freeF>>> OOEntry;

class ObjectOwner : public Singleton<ObjectOwner>
{
public:
    ObjectOwner(void);
    virtual ~ObjectOwner(void);

    OOEntry Find(XlfRef reference) const;
    OOEntry Mark(XlfRef reference, string handle);
    OOEntry Mark(XlfRef reference, OOEntry handles);
    void Free(XlfRef reference);

private:
    string toKey(XlfRef reference) const;
    map<string, OOEntry> m_owners;
};
