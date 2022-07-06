#include "ObjectOwner.h"

using namespace std;

ObjectOwner::ObjectOwner(void)
{
}

ObjectOwner::~ObjectOwner(void)
{
}

OOEntry ObjectOwner::Find(XlfRef reference) const
{
    auto it = m_owners.find(toKey(reference));
    if (it != m_owners.end())
    {
        return it->second;
    }
    else
    {
        return make_pair("", vector<pair<string, freeF>>());
    }
}

OOEntry ObjectOwner::Mark(XlfRef reference, string handle) 
 {
     return Mark(reference, make_pair(handle, vector<pair<string, freeF>>()));
}


OOEntry ObjectOwner::Mark(XlfRef reference, OOEntry handles) 
{
    auto oldHandle = Find(reference);
    m_owners[toKey(reference)] = handles;
    return oldHandle;
}

void ObjectOwner::Free(XlfRef reference)
{
    m_owners.erase(reference.GetTextA1());
}

string ObjectOwner::toKey(XlfRef reference) const
{
    std::ostringstream stream;
    stream << reference.GetSheetId() << reference.GetTextA1();
    return stream.str();
}
    