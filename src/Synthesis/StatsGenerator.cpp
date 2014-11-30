#include "StatsGenerator.h"

#include "Synthesis/TextureFilterer.h"

using namespace TextureSynthesis;

void StatsGenerator::computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
        const std::vector<std::vector<Signal>>& modulationSignals,
        std::vector<std::complex<double>>& statistics)
{
    computeStatistics(cochlearEnvelopes, modulationSignals, statistics, nullptr);
}

void StatsGenerator::computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
        const std::vector<std::vector<Signal>>& modulationSignals,
        std::vector<std::pair<double, int>>& statistics,
        std::vector<std::vector<std::complex<double>>> jacobian)
{
    computeStatistics(cochlearEnvelopes, modulationSignals, statistics, &jacobian);
}

void StatsGenerator::computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
        const std::vector<std::vector<Signal>>& modulationSignals,
        std::vector<std::pair<double, int>>& statistics,
        std::vector<std::vector<std::complex<double>>>* jacobian)
{
    // Cochlear statistics

    int signalLength = cochlearEnvelopes.at(0).size();
    double sampleRate = cochlearEnvelopes.at(0).sampleRate;

    int numEnvelopes = cochlearEnvelopes.size();

    std::vector<int> correlationDifs { 1, 2, 3, 5, 8, 11, 16, 21 };

    std::vector<double> cochlearMeans(numEnvelopes);
    std::vector<double> cochlearVariances(numEnvelopes);
    statistics.clear();
    if(jacobian) jacobian.clear();

    std::vector<double> zeros;
    for(int i = 0; i < cochlearEnvelopes[0].size(); i++)
        zeros.push_back(0);

    double mean, variance;
    std::vector<double> realPart;
    for(int i = 0; i < numEnvelopes; i++)
    {
        realPart = cochlearEnvelopes[i].realPart();

        mean = computeMean(realPart);
        cochlearMeans.push_back(mean);
        statistics.push_back(mean);
        if(jacobian)
        {
            std::vector<std::vector<double>> 1stMomentGrads;
            for(int env = 0; env < numEnvelopes; env++)
            {
                if(i == env)
                    1stMomentGrads.push_back(computeMeanGrads(realPart));
                else
                    1stMomentGrads.push_back(zeros);
            }
            jacobian.push_back(1stMomentGrads);
        }

        variance = computeVariance(realPart, mean);
        cochlearVariances.push_back(variance);
        statistics.push_back(nthMoment(realPart, 2, mean, variance));
        statistics.push_back(nthMoment(realPart, 3, mean, variance));
        statistics.push_back(nthMoment(realPart, 4, mean, variance));
        if(jacobian)
        {
            std::vector<std::vector<double>> 2ndMomentGrads;
            std::vector<std::vector<double>> 3rdMomentGrads;
            std::vector<std::vector<double>> 4thMomentGrads;
            for(int env = 0; env < numEnvelopes; env++)
            {
                if(i == env)
                {
                    2ndMomentGrads.push_back(nthMomentGrads(realPart, 2, mean, variance));
                    3rdMomentGrads.push_back(nthMomentGrads(realPart, 3, mean, variance));
                    4thMomentGrads.push_back(nthMomentGrads(realPart, 4, mean, variance));
                }
                else
                {
                    2ndMomentGrads.push_back(zeros);
                    3rdMomentGrads.push_back(zeros);
                    4thMomentGrads.push_back(zeros);
                }
            }
            jacobian.push_back(2ndMomentGrads);
            jacobian.push_back(3rdMomentGrads);
            jacobian.push_back(4thMomentGrads);
        }

        for(int j : correlationDifs)
            if(i - j >= 0)
            {
                statistics.push_back(crossCorrelation(realPart,
                            cochlearEnvelopes[i - j].realPart(), mean, cochlearMeans[i - j],
                            variance, cochlearVariances[i - j]));
                if(jacobian)
                {
                    std::vector<std::vector<double>> ijCorrelationGrads;
                    for(int env = 0; env < numEnvelopes; env++)
                    {
                        if(env == i || env == i - j)
                            ijCorrelationGrads.push_back(crossCorrelationGrads(realPart,
                                    cochlearEnvelopes[i - j].realPart(), mean, cochlearMeans[i - j],
                                    variance, cochlearVariances[i - j]), env == i - j);
                        else
                            ijCorrelationGrads.push_back(zeros);
                    }
                    jacobian.push_back(ijCorrelationGrads);
                }
            }
    }


    // TODO: continue correction jacobian-filling

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
            variance = computeVariance(realPart, computeMean(realPart));
            modulationVariances[i][j] = variance;

            statistics.push_back(computePower(realPart, variance));
            if(jacobian) jacobian.push_back(computePowerGrads(realPart, variance));
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
                if(jacobian) jacobian.push_back(c1ModulationCorrelationGrads(
                            modulationSignals[i][n].realPart(),
                            modulationSignals[j][n].realPart(),
                            modulationVariances[i][n], modulationVariances[j][n]));
            }

    int numC2s = 6;
    std::vector<Signal> analyticModSignals;
    for(int n = 0; n <= numC2s; n++)
        analyticModSignals.push_back(Signal(signalLength, sampleRate));

    for(int i = 0; i < numEnvelopes; i++)
    {
        for(int n = 0; n <= numC2s; n++)
        {
            analyticModSignals[n].set(modulationSignals[i][n]);
            analyticModSignals[n].makeAnalytic();
        }
        for(int n = 0; n < numC2s; n++)
        {
            statistics.push_back(c2ModulationCorrelation(analyticModSignals[n],
                    analyticModSignals[n + 1], modulationVariances[i][n],
                    modulationVariances[i][n + 1]));
            if(jacobian) jacobian.push_back(c2ModulationCorrelationGrads(analyticModSignals[n],
                    analyticModSignals[n + 1], modulationVariances[i][n],
                    modulationVariances[i][n + 1]));
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

std::vector<double> computeMeanGrads(const std::vector<double>& data)
{
    int size = data.size();
    double dm = 1.0 / size;
    std::vector<double> grads;
    for(int i = 0; i < data.size(); i++)
        grads.push_back(dm);
    return grads;
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

        default:
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

