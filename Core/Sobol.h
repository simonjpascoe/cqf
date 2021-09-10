#pragma once
#include <vector>

namespace Common {

class Sobol 
{
public:
	Sobol(unsigned long dimensions);
	~Sobol();

	std::vector<double> Draw();
	void Skip(unsigned long n);

private:
	unsigned long m_dimensions;

	static const int MAXDIM = 6;
	static const int MAXBIT = 30;
	unsigned long long m_in;

	std::vector<unsigned> m_ix;
	unsigned *m_iu[MAXBIT];

	int m_mdeg[MAXDIM];
	unsigned m_ip[MAXDIM];

	unsigned m_iv[MAXDIM*MAXBIT];
	double m_fac;

};

}