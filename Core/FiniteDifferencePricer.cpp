#include "stdafx.h"
#include <vector>
#include <functional>
#include <ppl.h>

#include "Matrix.h"
#include "Mathematics.h"
#include "FiniteDifferencePricer.h"
#include "Debugging.h"

using namespace Common;
using namespace concurrency;

namespace FD {

FiniteDifferenceResult FiniteDifferencePricer::Evaluate(
    const Enums::EvaluationMode evaluationMode,
    const vector<pair<FiniteDifferencePriceable*,double>>& exotic, 
    const Quotes<FiniteDifferencePriceable*>& hedgeInstruments,
    const FiniteDifferenceParameters& parameters,
    const MarketData& marketData) const
{
    // 1. Calculate the maximum maturity in order to initialise the grid. This isn't used
    // here except for returning to the user interface though.
    auto exoticMaturity = calculateMaturity(exotic);
    auto hedgeMaturity  = calculateMaturity(hedgeInstruments);

    auto space = CalculateSpaceVector(max(exoticMaturity, hedgeMaturity), marketData, parameters);
    pair<double, FiniteDifferenceSolution> offerSolution;
    pair<double, FiniteDifferenceSolution> bidSolution;

    // 2. Perform the main calculations. Optimization performs in parallel for the bid/offer to return
    // answers quicker (provided machine hardware is available)
    if (evaluationMode == Enums::Optimize) {
        parallel_invoke(
            [this, &offerSolution, &parameters, &exotic, &marketData, &hedgeInstruments] {
                dbg("offer");
                offerSolution = OptimizeWeights(parameters, exotic, marketData, hedgeInstruments, false);
                dbg("offer finished with quality %f", offerSolution.first);
                },
            [this, &bidSolution, &parameters, &exotic, &marketData, &hedgeInstruments] {
                dbg("bid");
                bidSolution = OptimizeWeights(parameters, exotic, marketData, hedgeInstruments, true);
                dbg("bid finished with quality %f", bidSolution.first);
                }
            );
    } 
    // 2b. Standard pricing evaluates via the same function but without any hedging instruments
    else if (evaluationMode == Enums::Pricing) {
        vector<double> emptyWeights;
        Quotes<FiniteDifferencePriceable*> emptyHedges;
        bidSolution.first = 0;
        offerSolution.first = 0;
        bidSolution.second = EvaluateHedge(parameters, exotic, marketData, emptyHedges, emptyWeights);
        offerSolution.second = EvaluateHedge(parameters, invertPosition(exotic), marketData, emptyHedges, emptyWeights);
        offerSolution.second.Value *= -1.0;
    } else throw ("Unsupported evaluation mode");

    // 3. Return full solution to user
    FiniteDifferenceResult fdr;
    fdr.Offer = offerSolution;
    fdr.Bid = bidSolution;
    fdr.SpotLevels = Matrix::FromVector(space);
    return fdr;
}

// The core top level evaluation function that is used for all pricing purposes
FiniteDifferenceSolution FiniteDifferencePricer::EvaluateHedge(
    const FiniteDifferenceParameters& parameters,
    const vector<pair<FiniteDifferencePriceable*,double>>& exotic, 
    const MarketData& marketData,
    const Quotes<FiniteDifferencePriceable*>& hedgeInstruments,
    const vector<double>& hedgeWeights) const
{
    // 1. Setup. Calculate the maturity and merge the exotic and hedging instruments
    // determine the maturity payoff of the exotic
    auto maturity = calculateMaturity(exotic);
    auto package = exotic;
    auto hedges = hedgeInstruments.All();

    for (auto i = 0; i < hedges.size(); i++) {
        package.push_back(make_pair(hedges[i].Instrument, hedgeWeights[i]));
    }

    // 2. Call the pricing Kernel, maybe using Richardson Extrapolation
    const double extrapolationFactor = 1.6;
    vector<double> result;
    switch (parameters.Scheme) {
        case Enums::Explicit :            { result = Kernel(parameters, package, marketData); break;}
        case Enums::ExplicitRichardson :{ result = RichardsonKernel(parameters, package, marketData, extrapolationFactor); break;}
        case Enums::Implicit :            { result = KernelImplicit(parameters, package, marketData); break;} 
    }

    // 3. Calculate the cost(/gain) of the hedge
    auto hedge =  0.0;
    for (auto i = 0; i < hedges.size(); i++) {
        hedge += hedges[i].Bid * hedgeWeights[i];
    }
    Matrix hedgeM(result.size(), 1);
    hedgeM.Fill(hedge);

    // 4. Calculate the space discretizaton for returning to user and to calculate maturity payoffs
    auto space = CalculateSpaceVector(calculateMaturity(package), marketData, parameters);

    // 5. Calculate the maturity payout for the exotic part. This only returns the sum of instruments maturing at the end.
    vector<double> maturityExoticPayoff(parameters.SpaceSteps+1,0.0);
    for (auto fdp : exotic)
    {
        if (fdp.first->Maturity() == maturity)
        {
            for (auto i=0; i<maturityExoticPayoff.size(); i++)
            {
                maturityExoticPayoff[i] += fdp.first->Payoff(space[i]) * fdp.second;
            }
        }
    }
    // 6. Determine the residual payoff of the entire package (any component before/after exotic maturity will be ignored)
    vector<double> residualMaturityPayoff(parameters.SpaceSteps+1,0.0);
    for (auto fdp : package)
    {
        if (fdp.first->Maturity() == maturity)
        {
            for (auto i=0; i<residualMaturityPayoff.size(); i++)
            {
                residualMaturityPayoff[i] += fdp.first->Payoff(space[i]) * fdp.second;
            }
        }
    }

    // 7. Group and return to user
    FiniteDifferenceSolution solution;
    solution.Payoff =  Matrix::FromVector(maturityExoticPayoff);
    solution.ResidualPayoff = Matrix::FromVector(residualMaturityPayoff);
    solution.Curve = Matrix::FromVector(result);
    solution.Value = interpolate(marketData.Spot, space, result) - hedge; 
    solution.HedgeWeights = Matrix::FromVector(hedgeWeights);

    return solution;
}

pair<double, FiniteDifferenceSolution> FiniteDifferencePricer::OptimizeWeights(
    const FiniteDifferenceParameters& parameters,
    const vector<pair<FiniteDifferencePriceable*,double>>& exotic, 
    const MarketData& marketData,
    const Quotes<FiniteDifferencePriceable*>& hedgeInstruments,
    const bool solveForBid) const
{
    // 0. Flip the exotic if we are solving for the offer price.
    auto exo = solveForBid ? exotic : invertPosition(exotic);

    // 1. Create the optimisation objective function. The lambda creates a function which is
    // f : R^N -> R, suitable for multidimensionl minimisation
    auto objective = [this, &parameters, &exo, &marketData, &hedgeInstruments] (const vector<double>& weights) 
        { 
            return -EvaluateHedge(parameters, exo, marketData, hedgeInstruments, weights).Value;
        };

    // 2. Perform optimisation. Start with initial guess of 1.0
    vector<double> initialWeights(hedgeInstruments.All().size(), 1.0 );

    try {
        auto optimal = minimize_bfgs(initialWeights, objective);

        // 3. Finally evaluate the optimal returned solution
        auto solution = EvaluateHedge(parameters, exo, marketData, hedgeInstruments, optimal.second);

        // adjust for inversion of the position for calculating offer
        if (!solveForBid)
        {
            solution.Value *= -1.0;
        }

        // 4. Return the quality of fit indicator and the solution
        return make_pair(optimal.first, solution);
    }
    catch (exception&)
    {
        throw "Optimisation has failed to find a suitable solution";
    }
}

//
// The PDE UVM Pricing kernel function (Explicit scheme)
//
vector<double> FiniteDifferencePricer::Kernel(
    const FiniteDifferenceParameters& parameters,
    const vector<pair<FiniteDifferencePriceable*, double>>& instruments,
    const MarketData& marketData) const
{
    // 1. Initialise space vector
    auto maturity = calculateMaturity(instruments);
    auto space = CalculateSpaceVector(maturity, marketData, parameters);
    auto spaceStep = space[1]-space[0];

    // 2. Initialise time step
    // For explicit finite difference, the timestep needs to be small to stay stable
    // as we *could* use high volatility all the time in pricing, we use the
    // high volatility for this calculation.
    // This takes into account that the grid may not start at S = 0.
    // Allow excessive override if needed...
    auto timeStep = min(parameters.MaxTimestepLength, 0.99 * spaceStep*spaceStep / space[parameters.SpaceSteps] / space[parameters.SpaceSteps]
                                                        /marketData.HighVolatility / marketData.HighVolatility); 
                         
    long timesteps = (long)ceil(maturity / timeStep + 0.5);
    auto timeStep2 = maturity/timesteps;

    dbg("Kernel timestep: %9f (total steps %d)", timeStep2, timesteps);

    // 3. Index all instruments by timestep to know when to add them into the pricing grid
    vector<pair<long, pair<FiniteDifferencePriceable*, double>>> indexedInstruments;
    for (auto instrument : instruments)
    {
        auto step = (long)floor(instrument.first->Maturity() / timeStep2);
        indexedInstruments.push_back(make_pair(step, instrument));
    }

    // 4. Rollback the pde solver. The entire solution is not kept in memory; only the timestep just solved. These
    // two vectors are then swapped to provide a good use of memory.
    auto_ptr<vector<double>> current(new vector<double>(parameters.SpaceSteps+1,0.0));
    auto_ptr<vector<double>> future(new vector<double>(parameters.SpaceSteps+1,0.0));
    
    // handle the possibly non-zero low bound for S.
    auto shift = space[0] / spaceStep;

    // 4a. Rollback. Start at timesteps+1, and then adding the maturity instruments and any early expiration
    // instruments can be handled by the same code without having to split the vector of instruments.
    for (auto i = timesteps+1; i>=0; i--)
    {
        // 4b. Update for options to be added to this timestep
        for (auto iis : indexedInstruments)
        {
            if (iis.first == i) 
            {
                auto instrument = iis.second;
                for (unsigned int k = 0; k <= parameters.SpaceSteps; k++)
                    {
                        future->at(k) += instrument.first->Payoff(space[k]) * instrument.second;
                    }
            }
        }

        // 4c. Lower boundary conditions
        if (shift == 0) {
            current->at(0) = future->at(0) * exp(-timeStep2 * marketData.InterestRate);
        } else {
            current->at(0) = 2*current->at(1) - current->at(2);
        }
        
        // 4d. Internal grid valuation
        for (unsigned long j = 1; j<parameters.SpaceSteps; j++)
        {
            auto s = shift + j;
            // 4e. Find out which vol to use by sign of gamma (we don't need the exact value)
            // (We treat gamma zero as high volatility)
            auto vol = (future->at(j+1) - 2 * future->at(j) + future->at(j-1)) > 0 ?
                            marketData.LowVolatility  : marketData.HighVolatility;
            current->at(j) = future->at(j) + timeStep2 * (
                + 0.5 * s * future->at(j+1) * (vol * vol * s + marketData.InterestRate)
                + 0.5 * s * future->at(j-1) * (vol * vol * s - marketData.InterestRate)
                - (vol * vol * s * s + marketData.InterestRate) * future->at(j));
        }
        
        // 4e. Upper boundary conditions
        current->at(parameters.SpaceSteps) = 2*current->at(parameters.SpaceSteps-1) - current->at(parameters.SpaceSteps-2);

        // 4f. Swap vectors and repeat until initial current time (t=0).
        swap(current,future);
    }

    vector<double> result(*current);
    return result;
}

//
// The PDE UVM Pricing kernel function (ExplicitRichardson scheme)
// Uses the normal Explicit solver internally
//
vector<double> FiniteDifferencePricer::RichardsonKernel(
    const FiniteDifferenceParameters& parameters,
    const vector<pair<FiniteDifferencePriceable*, double>>& instruments,
    const MarketData& marketData,
    const double factor) const
{
    dbg("RichardsonKernel");
    // 1. Calculate the specification for the bumped and original grids
    auto parameters2 = parameters; 
    parameters2.SpaceSteps = (long)ceil(parameters2.SpaceSteps * factor);

    auto maturity = calculateMaturity(instruments);

    auto space  = CalculateSpaceVector(maturity, marketData, parameters);
    auto space2 = CalculateSpaceVector(maturity, marketData, parameters2);

    auto dS1 = space[2]-space[1];
    auto dS2 = space2[2]-space2[1];
    auto dS1s = dS1*dS1;
    auto dS2s = dS2*dS2;
    dbg("dS1: %.9f, dS2: %.9f", dS1, dS2);

    // 2. Solve in parallel to speed up the solution (provided hardware is available). Further threads may be used
    // if combining optimisation with the extrapolation.
    vector<double> solution;
    vector<double> solution2;
    parallel_invoke(
        [&solution,  this, &parameters,  &instruments, &marketData] () { solution  = Kernel(parameters,  instruments, marketData);},
        [&solution2, this, &parameters2, &instruments, &marketData] () { solution2 = Kernel(parameters2, instruments, marketData);}
        );

    // 3. Interpolate the larger grid solution back to the smaller grid space vector
    // and merge the results to form the final solution
    for (auto i = 0; i<solution.size(); i++)
    {
        auto v1 = solution[i];
        auto v2 = interpolate(space[i], space2, solution2);
        auto n = (dS2s * v1 - dS1s * v2) / (dS2s - dS1s);
        solution[i] = n;
    }

    return solution;
}

//
// The PDE UVM Pricing kernel function (Implicit scheme)
//
vector<double> FiniteDifferencePricer::KernelImplicit(
        const FiniteDifferenceParameters& parameters,
        const vector<pair<FiniteDifferencePriceable*, double>>& instruments,
        const MarketData& marketData) const
{
    // 1. Initialise space vector
    auto maturity = calculateMaturity(instruments);
    auto space = CalculateSpaceVector(maturity, marketData, parameters);
    auto spaceStep = space[1]-space[0];

    // 2. Initialise time step (set to the maximum at the moment)
    auto timesteps = maturity/parameters.MaxTimestepLength;
    auto timeStep = parameters.MaxTimestepLength;

    dbg("Kernel timestep: %9f (total steps %d)", timeStep, timesteps);

    // 3. Index all instruments by timestep to know when to add them into the pricing grid
    vector<pair<long, pair<FiniteDifferencePriceable*, double>>> indexedInstruments;
    for (auto instrument : instruments)
    {
        auto step = (long)floor(instrument.first->Maturity() / timeStep);
        indexedInstruments.push_back(make_pair(step, instrument));
    }


    // 4. Rollback the pde solver. The entire solution is not kept in memory; only the timestep just solved. These
    // two vectors are then swapped to provide a good use of memory.
    auto M = Matrix::Identity(parameters.SpaceSteps+1);
    auto_ptr<Matrix> current(new Matrix(parameters.SpaceSteps+1,1));
    auto_ptr<Matrix> future(new Matrix(parameters.SpaceSteps+1,1));

    // handle the possibly non-zero low bound for S.
    auto shift = space[0] / spaceStep;

    // 4a. Rollback. Start at timesteps+1, and then adding the maturity instruments and any early expiration
    // instruments can be handled by the same code without having to split the vector of instruments.
    for (auto i = timesteps+1; i>=0; i--)
    {
        // 4b. Update for options to be added to this timestep
        for (auto iis : indexedInstruments)
        {
            if (iis.first == i) 
            {
                auto instrument = iis.second;
                for (unsigned int k = 0; k <= parameters.SpaceSteps; k++)
                    {
                        (*future)(k,0) += instrument.first->Payoff(space[k]) * instrument.second;
                    }
            }
        }

        auto vol = marketData.LowVolatility;

        // 4c. Set up the matrix at time i
        // boundary conditions
        M(0,0) = 2;
        M(0,1) = -1.0;

        M(parameters.SpaceSteps,parameters.SpaceSteps) = 2.0;
        M(parameters.SpaceSteps,parameters.SpaceSteps-1) = -1.0;

        // non - boundary conditions
        for (auto j=1; j<parameters.SpaceSteps; j++)
        {
            auto J = shift*vol*vol * (2*j + shift);
            auto a = 0.5 * timeStep * (marketData.InterestRate * j - vol*vol * j * j + marketData.InterestRate * shift - J);
            auto b = 1.0 + timeStep * (marketData.InterestRate + vol*vol * j * j  + J);
            auto c = 0.5 * timeStep * (- marketData.InterestRate * j - vol*vol * j * j - marketData.InterestRate * shift - J);
            M(j,j-1) = a;
            M(j,j) = b;
            M(j,j+1) = c;
        }

        *current = M.TriSolve(*future);


        // 4f. Swap vectors and repeat until initial current time (t=0).
        swap(current,future);
    }

    vector<double> result = current->ToVector(0);
    return result;
}

const vector<double> FiniteDifferencePricer::CalculateSpaceVector(
        const double maturity,
        const MarketData& marketData,
        const FiniteDifferenceParameters& parameters) const
{
    // Based on the volatility and maturity, expand a grid centered on the current spot (floored at S=0 as required)
    auto spaceWidth = marketData.HighVolatility * sqrt(maturity) * marketData.Spot;
    auto minSpaceVariable = max(0.0, marketData.Spot - (spaceWidth * parameters.GridWidth));
    auto maxSpaceVariable =  marketData.Spot + spaceWidth * parameters.GridWidth;
    auto spaceStep = (maxSpaceVariable - minSpaceVariable) / parameters.SpaceSteps;

    vector<double> spots;
    for (unsigned int i = 0; i <= parameters.SpaceSteps; i++)
    {
        auto level = i*spaceStep;
        spots.push_back(minSpaceVariable + i*spaceStep);
    }

    return spots;
}

double FiniteDifferencePricer::calculateMaturity(
    const vector<pair<FiniteDifferencePriceable*, double>>& instruments) const
{
    auto maturity = 0.0;
    for (auto fdp : instruments)
    {
        maturity = max(maturity, fdp.first->Maturity());
    }
    return maturity;
}

double FiniteDifferencePricer::calculateMaturity(
    const vector<FiniteDifferencePriceable*>& instruments) const
{
    auto maturity = 0.0;
    for (auto fdp : instruments)
    {
        maturity = max(maturity, fdp->Maturity());
    }
    return maturity;
}

double FiniteDifferencePricer::calculateMaturity(
    const Quotes<FiniteDifferencePriceable*>& quotes) const
{
    auto maturity = 0.0;
    for (auto fdp : quotes.All())
    {
        maturity = max(maturity, fdp.Instrument->Maturity());
    }
    return maturity;
}

vector<pair<FiniteDifferencePriceable*, double>> FiniteDifferencePricer::invertPosition( 
        const vector<pair<FiniteDifferencePriceable*, double>>& instruments) const
{
    auto flip = instruments;
    for (auto i = 0; i<flip.size(); i++)
    {
        flip[i].second *= -1;
    }
    return flip;
}
}