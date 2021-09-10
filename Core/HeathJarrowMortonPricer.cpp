#include "stdafx.h"
#include <vector>
#include <functional>
#include <ppl.h>

#include "HeathJarrowMortonPricer.h"
#include "Sobol.h"
#include "Mathematics.h"

using namespace std;

namespace HJM {

CalibrationResult HeathJarrowMortonPricer::CalibrateFactors(
	const Matrix& maturities_,
	const Matrix& timeSeriesData, 
	double coverageFactor) const
{
	// 0. check inputs
	// maturities should be column vector
	auto maturities = maturities_;
	if (maturities_.Cols() != 1) { maturities = maturities_.Transpose(); }

	// 1. calculate differences matrix on the series data
	Matrix differences(timeSeriesData.Rows()-1, timeSeriesData.Cols());
	for (int c = 0; c<timeSeriesData.Cols(); ++c) {
		double columnSum = 0;
		for (int r = 0; r<timeSeriesData.Rows()-1; ++r) {
			differences(r,c) = timeSeriesData(r+1,c) - timeSeriesData(r,c);
		}
	}

	// 2. calculate the covariance matrix
	// transpose because the Matrix class is Row major 
    // so lookups are more efficient in this direction
	// then do the covar calc on rows rather than columns for same result
	auto differencesT = differences.Transpose();
	Matrix covariance(differencesT.Rows(),differencesT.Rows());
	for (int i = 0; i<differencesT.Rows(); ++i)
	{
		auto primary = differencesT.Row(i);
		for (int j = 0; j<differencesT.Rows(); ++j)
		{
			if (i<=j) {
				auto secondary = differencesT.Row(j);
				covariance(i,j) = calculateCovariance(primary, secondary) * 252;
			}
			else
			{
				covariance(i,j) = covariance(j,i);
			}
		}
	}

	// 3. Calculate the PCA. covar matrix is SPD so no need to transpose from the above result
	auto pca = covariance.PrincipalComponentAnalysis(coverageFactor);
	
	// 4. Fit the volatility functions
	auto project = createBasisFunctions(1,1,1,1);

	Matrix projections(maturities.Rows(),4);
	for (int i = 0; i < maturities.Rows(); i++)
	{
		auto p = project(maturities(i,0));
		projections(i,0) = p[0];
		projections(i,1) = p[1];
		projections(i,2) = p[2];
		projections(i,3) = p[3];
	}
	
	// 4a. Fit the curves, return the basis function coefficients
	Matrix parameters(4,pca.Eigenvectors.size());
	for (int i=0; i<pca.Eigenvectors.size(); i++)
	{
		// (reuse memory)
		pca.Eigenvectors[i] = pca.Eigenvectors[i] * sqrt(pca.Eigenvalues[i]);
		auto result = projections.fit(pca.Eigenvectors[i]);
		for (int j=0; j<4; j++) {
			parameters(j, i) = result(j,0);
		}
	}

	// setup return struct
	CalibrationResult cr;
	cr.Coefficients = parameters;
	cr.AnalysisResult = pca;
	return cr;
}

PricingResult HeathJarrowMortonPricer::Simulate(
	const vector<HeathJarrowMortonPriceable*>& products,
	const Matrix& coefficients, 
	const Matrix& tenors, 
	const Matrix& seeds, 
	double maturity, 
	double dt,
	long simulationCount) const
{
	// build the volatility functions
	vector<function<double(double)>> vol;
	for (int i=0; i<coefficients.Cols(); i++) {
		// fast(er) but fixes the functional form to up to quartic level only.
		auto volf = [i,&coefficients] (double x)
			{
				return coefficients(0,i) + x * (coefficients(1,i) + x * (coefficients(2,i) + coefficients(3,i)*x));
			};
		vol.push_back(volf);
	}

	// thus the volatility given by the first vol function at t = 0.4 is 
	// auto v = vol[0](0.4);
	// here we are assuming that in the standard model, v(t) = v(0,t) c.f. slide notes

	// in the Musiela parameterisation, define the drift function m
	// m(t) = \sum_{i=0}^N v_i (tau) \int_0^tau v_i(s) ds
	auto m = [this,&vol] (double tau) 
		{
			double m = 0;
			for (auto vi : vol) {
				m += vi(tau) * integrate(vi,0,tau);
			}
			return m;
		};

	size_t resolution = (size_t)ceil(maturity/dt);
	// define the dF function which takes a pre-calculated m value
	auto dF = [&vol,dt] (double m_tau, double tau, const vector<double>& dX) 
		{
			double volsum = 0;
			for (int i =0; i<vol.size(); i++) { 
				volsum+= vol[i](tau)*dX[i];
			}
			return m_tau*dt +volsum * sqrt(dt);
		};

	// perform pricing and simulation
	
	
	// Optional random skip for Sobol generator
	//int r = rand();
	//for (int j=0; j<r; j++) { sobolGen.Draw(vol.size());}
	
	Matrix simulation(resolution+1, tenors.Cols());

	// 1. seed
	// 2. cache the m(tau) function values as they remain the same for each maturity, rather than keep calculating them
	if (tenors.Cols() != seeds.Cols()) { throw "Tenors and Seed should have the same number of columns"; }
	vector<double> mcache;
	for (auto i=0; i<seeds.Cols(); i++) { 
		simulation(0,i) = seeds(0,i); 
		mcache.push_back(m(tenors(0,i)));
	}
	
	// for the convergence plot, find all powers of 2 < simulationCount
	// and the convergence results follow including the final count as well
	auto columns = (long)floor(sqrt(simulationCount)) + 1;
	if (sqrt(simulationCount) == (long)floor(sqrt(simulationCount))) { columns--;} // if simcount is square, no need for extra column
	Matrix convergence(products.size()+1,columns);
	Matrix fairvalues(products.size(), simulationCount+1);

	// precalcuate skipped Sobol generators
	Sobol sobolCore((unsigned long)vol.size());
	sobolCore.Skip(4096); // skip the first part of the series for all cases

	vector<Sobol> sobolGens;
	for (long sim=1; sim<=simulationCount; sim++)
	{
		Sobol sobolChild(sobolCore);
		sobolGens.push_back(sobolChild);
		sobolCore.Skip(resolution+1);
	}

	// lambda to execution one simulation result
	auto simulate = 
		[&] (long sim) 
		{
			// 3. evolve and price one simulation
			auto sobolGen = sobolGens[sim-1];
			for (auto i=1; i<=resolution; i++) {
				// 3a. draw 'random' deviates, uncorrelated
				vector<double> dX;
				for (auto u: sobolGen.Draw()) { dX.push_back(InverseStandardCumulativeNormal(u));}

				// 3b. most maturities fowards diff for dFbar/dtau
				for (auto j=0; j<simulation.Cols()-1; j++)
				{
					simulation(i,j) = simulation(i-1,j) + dF(mcache[j],tenors(0,j),dX) + (simulation(i-1,j+1)- simulation(i-1,j))/(tenors(0,j+1)-tenors(0,j)) * dt;
				}
				// 3c. final maturity - backwards diff for dFbar/dtau
				auto j = simulation.Cols()-1;
				simulation(i,j) = simulation(i-1,j) + dF(mcache[j],tenors(0,j),dX) + (simulation(i-1,j)- simulation(i-1,j-1))/(tenors(0,j)-tenors(0,j-1)) * dt;
			}

			// 4. update product prices
			for (auto k=0; k<products.size(); k++) {
				fairvalues(k,sim) = products[k]->HeathJarrowMortonEvaluate(dt, simulation, tenors);
			}
		};

	concurrency::parallel_for((long)1, simulationCount+1, simulate);
			
	Matrix finals(products.size(), 1);
	// construct simulation convergence matrix
	for (auto sim = 1; sim<simulationCount+1; sim++)
	{
		for (auto k=0; k<products.size(); k++) {
			finals(k,0) = ((sim-1) * finals(k,0) + fairvalues(k,sim))/sim;
			auto n = (long)floor(sqrt(sim));
			if (sqrt(sim) == n) { // i.e. sim = n*n for some n 
				convergence(0,n-1) = (double)sim;
				convergence(k+1,n-1) = finals(k,0);
			}
		}
	}

	// 5. complete convergence result by the copying of the results to final column of matrix
	for (auto k=0; k<products.size(); k++) {
		convergence(k+1,columns-1) = finals(k,0);
	}
	convergence(0,columns-1)=simulationCount;

	PricingResult pr;
	pr.Convergence = convergence;
	pr.FairValues = finals;

	return pr;
}

double HeathJarrowMortonPricer::calculateCovariance(const Matrix& primary, const Matrix& secondary) const
{
	if (primary.Rows() != 1) { throw "Provided primary is not a row vector";}
	if (secondary.Rows() != 1) { throw "Provided secondary is not a row vector";}
	if (primary.Cols() != secondary.Cols()) { throw "Vectors for covariance calculation are not same size";}

	double pAverage = average(primary);
	double sAverage = average(secondary);
	double covariance = 0;

	for (int i=0; i<primary.Cols(); ++i)
	{
		covariance += (primary(0,i) - pAverage) * (secondary(0,i) - sAverage);
	}

	return covariance / primary.Cols();
}

double HeathJarrowMortonPricer::average(const Matrix& vector) const
{
	if (vector.Rows() != 1) {throw "Provided vector is not a row vector";}

	double sum = 0;
	for (int i=0; i<vector.Cols(); ++i)
	{
		sum+=vector(0,i);
	}
	return sum / vector.Cols();
}

function<vector<double>(double)> HeathJarrowMortonPricer::createBasisFunctions(double y1, double y2, double y3, double y4) const
{
	return [y1,y2,y3,y4] (double x) 
	{
		vector<double> projections;
		projections.push_back(y1);
		projections.push_back(y2*x);
		projections.push_back(y3*x*x);
		projections.push_back(y4*x*x*x);
		return projections;
	};
}

double HeathJarrowMortonPricer::integrate(const function<double(double)>& f, double lower, double upper) const
{
	// simple algorithm for now : Simpson's rule
	// divide up the region into 50 points
	auto stepsize = (upper - lower) / 50.0;
	auto stepsizeh = stepsize/2.0;
	auto coeff = stepsize / 6.0;

	auto integral = 0.0;
	for (auto i = 0; i<50; i++)
	{
		auto left = i*stepsize;
		integral += coeff * (f(left) + 4*f(left + stepsizeh) + f(left + stepsize));
	}
	return integral;
}

}