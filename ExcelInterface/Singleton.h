#pragma once

template<typename T>
class Singleton
{
public:
	static T& Instance()
	{
		static T object;
		return object;
	}
	virtual ~Singleton() { }
protected:
	Singleton() { }
private:
	Singleton(const Singleton& other);
	Singleton& operator=(const Singleton& other);
};