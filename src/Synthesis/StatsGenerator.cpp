#include "StatsGenerator.h"

#include "Synthesis/TextureFilterer.h"

using namespace TextureSynthesis;

void StatsGenerator::computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
        const std::vector<std::vector<Signal>>& modulationSignals,
        std::vector<double>& statistics)
{
    computeStatistics(cochlearEnvelopes, modulationSignals, statistics, nullptr);
}

void StatsGenerator::computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
        const std::vector<std::vector<Signal>>& modulationSignals,
        std::vector<double>& statistics,
        std::vector<std::vector<std::vector<double>>> jacobian)
{
    computeStatistics(cochlearEnvelopes, modulationSignals, statistics, &jacobian);
}

void StatsGenerator::computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
        const std::vector<std::vector<Signal>>& modulationSignals,
        std::vector<double>& statistics,
        std::vector<std::vector<std::vector<double>>>* jacobian)
{
    // TODO: optimize vector handling
    // TODO: finish implementing gradient functions

    // Cochlear statistics

    int signalLength = cochlearEnvelopes.at(0).size();
    double sampleRate = cochlearEnvelopes.at(0).sampleRate;

    int numEnvelopes = cochlearEnvelopes.size();

    std::vector<int> correlationDifs { 1, 2, 3, 5, 8, 11, 16, 21 };

    std::vector<double> cochlearMeans(numEnvelopes);
    std::vector<double> cochlearVariances(numEnvelopes);
    statistics.clear();
    if(jacobian) jacobian->clear();

    std::vector<double> zeros;
    for(int i = 0; i < cochlearEnvelopes[0].size(); i++)
        zeros.push_back(0);

    double mean, variance;
    std::vector<double> realPart;
    for(int i = 0; i < numEnvelopes; i++)
    {
        realPart = cochlearEnvelopes[i].realPart();

        std::vector<double> moments;
        mean = computeMean(realPart);
        cochlearMeans.push_back(mean);
        moments.push_back(mean);
        if(jacobian)
        {
            std::vector<std::vector<double>> firstMomentGrad;
            for(int env = 0; env < numEnvelopes; env++)
            {
                if(i == env)
                    firstMomentGrad.push_back(computeMeanGrad(realPart));
                else
                    firstMomentGrad.push_back(zeros);
            }
            jacobian->push_back(firstMomentGrad);
        }

        variance = computeVariance(realPart, mean);
        cochlearVariances.push_back(variance);
        moments.push_back(nthMoment(realPart, 2, mean, variance));
        moments.push_back(nthMoment(realPart, 3, mean, variance));
        moments.push_back(nthMoment(realPart, 4, mean, variance));
        for(int i = 0; i < 4; i++)
            statistics.push_back(moments[i]);   // find method for this when there's internet
        if(jacobian)
        {
            std::vector<std::vector<double>> secondMomentGrad;
            std::vector<std::vector<double>> thirdMomentGrad;
            std::vector<std::vector<double>> fourthMomentGrad;
            for(int env = 0; env < numEnvelopes; env++)
            {
                if(i == env)
                {
                    secondMomentGrad.push_back(nthMomentGrad(realPart, 2, variance, moments));
                    thirdMomentGrad.push_back(nthMomentGrad(realPart, 3, variance, moments));
                    fourthMomentGrad.push_back(nthMomentGrad(realPart, 4, variance, moments));
                }
                else
                {
                    secondMomentGrad.push_back(zeros);
                    thirdMomentGrad.push_back(zeros);
                    fourthMomentGrad.push_back(zeros);
                }
            }
            jacobian->push_back(secondMomentGrad);
            jacobian->push_back(thirdMomentGrad);
            jacobian->push_back(fourthMomentGrad);
        }

        for(int j : correlationDifs)
            if(i - j >= 0)
            {
                statistics.push_back(crossCorrelation(realPart,
                            cochlearEnvelopes[i - j].realPart(), mean, cochlearMeans[i - j],
                            variance, cochlearVariances[i - j]));
                if(jacobian)
                {
                    std::vector<std::vector<double>> ijCorrelationGrad;
                    for(int env = 0; env < numEnvelopes; env++)
                    {
                        if(env == i || env == i - j)
                            ijCorrelationGrad.push_back(crossCorrelationGrad(realPart,
                                    cochlearEnvelopes[i - j].realPart(), mean, cochlearMeans[i - j],
                                    variance, cochlearVariances[i - j], env == i));
                        else
                            ijCorrelationGrad.push_back(zeros);
                    }
                    jacobian->push_back(ijCorrelationGrad);
                }
            }
    }


    // Modulation statistics

    int modulationSize = modulationSignals[0].size();

    double defaultVariance = -1;
    std::vector<std::vector<double>> modulationVariances(numEnvelopes,
            std::vector<double>(modulationSize, defaultVariance));

    for(int i = 0; i < numEnvelopes; i++)
    {
        for(int j = 0; j < modulationSize; j++)
        {
            realPart = modulationSignals[i][j].realPart();
            mean = computeMean(realPart);
            variance = computeVariance(realPart, mean);
            modulationVariances[i][j] = variance;

            statistics.push_back(computePower(realPart, variance));
            if(jacobian)
            {
                std::vector<std::vector<double>> powerGrad;
                for(int env = 0; env < numEnvelopes; env++)
                {
                    if(env == i)
                        // TODO: pass appropriate filter to computePowerGrad
                        powerGrad.push_back(computePowerGrad(modulationSignals[i][j], mean, variance));
                    else
                        powerGrad.push_back(zeros);
                }
                jacobian->push_back(powerGrad);
            }
        }
    }

    // magic numbers here are presecribed in paper
    for(int i = 2; i < numEnvelopes; i++)
        for(int j = i - 1; j >= i - 2; j--)
            for(int n = 1; n < 7; n++)
            {
                statistics.push_back(c1ModulationCorrelation(
                            modulationSignals[i][n].realPart(),
                            modulationSignals[j][n].realPart(),
                            modulationVariances[i][n], modulationVariances[j][n]));
                if(jacobian)
                {
                    std::vector<std::vector<double>> ijnC1CorrelationGrad;
                    for(int env = 0; env < numEnvelopes; env++)
                    {
                        if(env == i || env == j)
                            ijnC1CorrelationGrad.push_back(c1ModulationCorrelationGrad(
                                    modulationSignals[i][n].realPart(),
                                    modulationSignals[j][n].realPart(),
                                    modulationVariances[i][n], modulationVariances[j][n],
                                    env == i));
                        else
                            ijnC1CorrelationGrad.push_back(zeros);
                    }
                    jacobian->push_back(ijnC1CorrelationGrad);
                }
            }

    int numC2s = 6;
    std::vector<Signal> analyticModSignals;
    for(int n = 0; n <= numC2s; n++)
        analyticModSignals.push_back(Signal(signalLength, sampleRate));

    std::complex<double> c2Corr;
    std::vector<std::complex<double>> c2CorrGrad;
    for(int i = 0; i < numEnvelopes; i++)
    {
        for(int n = 0; n <= numC2s; n++)
        {
            analyticModSignals[n].set(modulationSignals[i][n]);
            analyticModSignals[n].makeAnalytic();
        }
        for(int n = 0; n < numC2s; n++)
        {
            c2Corr = c2ModulationCorrelation(analyticModSignals[n],
                    analyticModSignals[n + 1], modulationVariances[i][n],
                    modulationVariances[i][n + 1]);
            statistics.push_back(std::real(c2Corr));
            statistics.push_back(std::imag(c2Corr));
            if(jacobian)
            {
                std::vector<std::vector<double>> inC2CorrelationRealGrad(numEnvelopes);
                std::vector<std::vector<double>> inC2CorrelationImagGrad(numEnvelopes);
                for(int env = 0; env < numEnvelopes; env++)
                {
                    if(env == n || env == n + 1)
                    {
                        c2CorrGrad = c2ModulationCorrelationGrad(analyticModSignals[n],
                            analyticModSignals[n + 1], modulationVariances[i][n],
                            modulationVariances[i][n + 1], env == n);
                        for(std::complex<double> c2Val : c2CorrGrad)
                        {
                            inC2CorrelationRealGrad[env].push_back(std::real(c2Val));
                            inC2CorrelationImagGrad[env].push_back(std::imag(c2Val));
                        }
                    }
                    else
                    {
                        inC2CorrelationRealGrad[env] = zeros;
                        inC2CorrelationImagGrad[env] = zeros;
                    }
                }
                jacobian->push_back(inC2CorrelationRealGrad);
                jacobian->push_back(inC2CorrelationImagGrad);
            }
        }
    }
}

double StatsGenerator::computeMean(const std::vector<double>& data)
{
    int size = data.size();
    double m = 0;
    for(int i = 0; i < size; i++)
        m += data[i];
    return m / size;
}

std::vector<double> StatsGenerator::computeMeanGrad(const std::vector<double>& data)
{
    int size = data.size();
    double dm = 1.0 / size;
    std::vector<double> grad;
    for(int i = 0; i < data.size(); i++)
        grad.push_back(dm);
    return grad;
}

double StatsGenerator::computeVariance(const std::vector<double>& data, double mean)
{
    int size = data.size();
    double v = 0;
    for(int i = 0; i < size; i++)
        v += (data[i] - mean) * (data[i] - mean);
    return v / size;
}

double StatsGenerator::nthMoment(const std::vector<double>& data,
        int n, double mean, double variance)
{
    switch(n)
    {
        case 1:
            return mean;

        case 2:
            if(mean < 1e-100)
                return 0;
            else
                return variance * variance / (mean * mean);
    }

    // default
    int size = data.size();
    double m = 0;
    for(int i = 0; i < size; i++)
        m += pow(data[i] - mean, n);
    double powVar = pow(variance, n);
    if(powVar < 1e-100)
        return 0;
    else
        return m / (size * powVar);
}

std::vector<double> StatsGenerator::nthMomentGrad(const std::vector<double>& data,
        int n, double variance, const std::vector<double>& moments)
{
    int size = data.size();
    std::vector<double> grad(size);
    double sizeInv = 1.0 / size;

    switch(n)
    {
        case 1:
            return computeMeanGrad(data);

        case 2:
        {
            double sqrtM = moments[1];
            for(int i = 0; i < size; i++)
                grad[i] = sizeInv * (1 - sizeInv) * (data[i] - moments[0]) / (moments[0] * sqrtM) -
                    sizeInv * sqrtM / pow(moments[0], 2);
            return grad;
        }

        case 3:
            for(int i = 0; i < size; i++)
                grad[i] = 3 / pow(moments[1], 1.5) * sizeInv * (pow(data[i] - moments[0], 2)
                        - moments[1] - moments[2] / moments[1] * (data[1] - moments[0]));
            return grad;

        case 4:
            for(int i = 0; i < size; i++)
                grad[i] = 4 * sizeInv * (pow(data[i] - moments[0], 3) - moments[2] -
                        moments[3] / moments[1] * (data[i] - moments[0])) / pow(moments[1], 2);
            return grad;
    }

    // default
    std::cerr << "Cannot compute gradient above 4th moment\n";
    exit(1);
}

double StatsGenerator::crossCorrelation(const std::vector<double>& data1,
        const std::vector<double>& data2, double mean1, double mean2,
        double variance1, double variance2)
{
    int size = std::min(data1.size(), data2.size());

    double c = 0;
    for(int i = 0; i < size; i++)
        c += (data1[i] - mean1) * (data2[i] - mean2);

    if(variance1 < 1e-10 || variance2 < 1e-10)
        return 0;
    else
        return c / (size * variance1 * variance2);
}

std::vector<double> StatsGenerator::crossCorrelationGrad(const std::vector<double>& data1,
        const std::vector<double>& data2, double mean1, double mean2, 
        double variance1, double variance2, bool varyingData1)
{
    int size = std::min(data1.size(), data2.size());
    std::vector<double> grad(size);
    double sizeInv = 1.0 / size;

    const std::vector<double> *varying, *fixed;
    double varyingMean, fixedMean, varyingVar, fixedVar;
    if(varyingData1)
    {
        varying = &data1; fixed = &data2;
        varyingMean = mean1; fixedMean = mean2;
        varyingVar = variance1; fixedVar = variance2;
    }
    else
    {
        varying = &data2; fixed = &data1;
        varyingMean = mean2; fixedMean = mean1;
        varyingVar = variance2; fixedVar = variance1;
    }

    double varyingStdDev = sqrt(varyingVar);
    double fixedStdDev = sqrt(fixedVar);

    double cbc = 0;
    for(int i = 0; i < size; i++)
        cbc += (varying->at(i) - varyingMean) * (fixed->at(i) - fixedMean);
    cbc /= size * varyingStdDev * fixedStdDev;

    for(int i = 0; i < size; i++)
        grad.push_back(
            (fixed->at(i) - fixedMean) * sizeInv / (varyingStdDev * fixedStdDev)
            - (varying->at(i) - varyingMean) * sizeInv / pow(varyingStdDev, 2)
            * cbc
        );

    return grad;
}

double StatsGenerator::computePower(const std::vector<double>& data, double variance)
{
    int size = data.size();

    double p = 0;
    for(int i = 0; i < size; i++)
        p += data[i] * data[i];

    if(variance < 1e-100)
        return 0;
    else
        return p;
}

std::vector<double> StatsGenerator::computePowerGrad(const Signal& data, const Filter& filter,
        double mean, double variance)
{
    int size = data.size();
    std::vector<double> grad;
    double sizeInv = 1.0 / size;

    Signal fData(data);
    filter.filter(fData);
    double sqFilteredMean = 0;
    for(int i = 0; i < size; i++)
    {
        sqFilteredMean += pow(std::real(fData[i]), 2);
        fData[i] *= sizeInv;
    }
    sqFilteredMean *= sizeInv;
    
    Signal ffData(fdata);
    filter.filter(ffdata);

    for(int i = 0; i < size; i++)
    {
        grad[i] = (2 * variance * ffData - 2 * sqFilteredMean
            * (std::real(data[i]) - mean) * sizeInv) / pow(variance, 2);
    }

    return grad;
}

double StatsGenerator::c1ModulationCorrelation(const std::vector<double>& data1,
        const std::vector<double>& data2, double variance1, double variance2)
{
    int size = std::min(data1.size(), data2.size());

    double c = 0;
    for(int i = 0; i < size; i++)
        c += data1[i] * data2[i];
    if(variance1 < 1e-10 || variance2 < 1e-10)
        return 0;
    else
        return c / (size * variance1 * variance2);
}

std::vector<double> StatsGenerator::c1ModulationCorrelationGrad(
        const std::vector<double>& data1, const std::vector<double>& data2,
        double variance1, double variance2, bool varyingData1)
{
    std::vector<double> grad(data1.size());
    return grad;
}

std::complex<double> StatsGenerator::c2ModulationCorrelation(const Signal& signal1,
        const Signal& signal2, double variance1, double variance2)
{
    int size = std::min(signal1.size(), signal2.size());

    std::complex<double> c = 0;
    for(int i = 0; i < size; i++)
    {
        double denom = std::norm(signal1[i]);
        if(denom > 1e-100)
            c += std::conj(std::pow(signal1[i], 2) /
                    denom) * signal2[i];
    }
    if(variance1 < 1e-10 || variance2 < 1e-10)
        return 0;
    else
        return c / (size * variance1 * variance2);
}

std::vector<std::complex<double>> StatsGenerator::c2ModulationCorrelationGrad(const Signal& signal1,
        const Signal& signal2, double variance1, double variance2, bool varyingData1)
{
    std::vector<std::complex<double>> grad(signal1.size());
    return grad;
}
