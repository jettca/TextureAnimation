#include "StatsGenerator.h"

#include "Synthesis/TextureFilterer.h"

using namespace TextureSynthesis;

void StatsGenerator::computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
        const std::vector<std::vector<Signal>>& modulationSignals,
        std::vector<std::complex<double>>& statistics)
{
    // Cochlear statistics

    int cochlearSize = cochlearEnvelopes.size();

    std::vector<int> correlationDifs { 1, 2, 3, 5, 8, 11, 16, 21 };

    std::vector<double> cochlearMeans(cochlearSize);
    std::vector<double> cochlearVariances(cochlearSize);
    statistics.clear();

    double mean, variance;
    std::vector<double> realPart;
    for(int i = 0; i < cochlearSize; i++)
    {
        realPart = cochlearEnvelopes[i].realPart();

        mean = computeMean(realPart);
        cochlearMeans.push_back(mean);
        statistics.push_back(mean);

        variance = computeVariance(realPart, mean);
        cochlearVariances.push_back(variance);
        statistics.push_back(nthMoment(realPart, 2, mean, variance));
        statistics.push_back(nthMoment(realPart, 3, mean, variance));
        statistics.push_back(nthMoment(realPart, 4, mean, variance));

        for(int j : correlationDifs)
            if(i - j >= 0)
                statistics.push_back(crossCorrelation(realPart,
                            cochlearEnvelopes[i - j].realPart(), mean, cochlearMeans[i - j],
                            variance, cochlearVariances[i - j]));
    }


    // Modulation statistics

    int modulationSize = modulationSignals[0].size();

    double defaultVariance = -1;
    std::vector<std::vector<double>> modulationVariances(cochlearSize,
            std::vector<double>(modulationSize, defaultVariance));

    for(int i = 0; i < cochlearSize; i++)
    {
        for(int j = 0; j < modulationSize; j++)
        {
            realPart = modulationSignals[i][j].realPart();
            variance = computeVariance(realPart, computeMean(realPart));
            modulationVariances[i][j] = variance;

            statistics.push_back(computePower(realPart, variance));
        }
    }

    // magic numbers here are presecribed in paper
    for(int i = 2; i < cochlearSize; i++)
        for(int j = i - 1; j >= i - 2; j--)
            for(int n = 1; n < 7; n++)
                statistics.push_back(c1ModulationCorrelation(
                            modulationSignals[i][n].realPart(),
                            modulationSignals[j][n].realPart(),
                            modulationVariances[i][n], modulationVariances[j][n]));
    for(int i = 0; i < cochlearSize; i++)
        for(int n = 0; n < 6; n++)
        {
            statistics.push_back(c2ModulationCorrelation(modulationSignals[i][n],
                    modulationSignals[i][n + 1], modulationVariances[i][n],
                    modulationVariances[i][n + 1]));
        }
}

double StatsGenerator::computeMean(const std::vector<double>& data)
{
    int size = data.size();
    double m;
    for(int i = 0; i < size; i++)
        m += data[i];
    return m / size;
}

double StatsGenerator::computeVariance(const std::vector<double>& data, double mean)
{
    int size = data.size();
    double v;
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

std::complex<double> StatsGenerator::c2ModulationCorrelation(Signal signal1, Signal signal2,
        double variance1, double variance2)
{
    int size = std::min(signal1._signal.size(), signal2._signal.size());

    signal1.makeAnalytic();
    signal2.makeAnalytic();

    std::complex<double> c = 0;
    for(int i = 0; i < size; i++)
    {
        double denom = std::norm(signal1._signal[i]);
        if(denom > 1e-100)
            c += std::conj(std::pow(signal1._signal[i], 2) /
                    denom) * signal2._signal[i];
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

