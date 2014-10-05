#include "AuditoryModel.h"
#include "CochlearFilter.h"

#include <cmath>
#include <complex>

using namespace TextureSynthesis;

double computeMean(const std::vector<double>& data)
{
    int size = data.size();
    double m;
    for(int i = 0; i < size; i++)
        m += data[i];
    return m / size;
}

double computeVariance(const std::vector<double>& data, double mean)
{
    int size = data.size();
    double v;
    for(int i = 0; i < size; i++)
        v += (data[i] - mean) * (data[i] - mean);
    return v / size;
}

double nthMoment(const std::vector<double>& data, int n, double mean, double variance)
{
    switch(n)
    {
        case 1:
            return mean;

        case 2:
            return variance * variance / (mean * mean);

        default:
            int size = data.size();
            double m = 0;
            for(int i = 0; i < size; i++)
                m += pow(data[i] - mean, n);
            return m / (size * pow(variance, n));
    }
}

double crossCorrelation(const std::vector<double>& data1,
        const std::vector<double>& data2, double mean1, double mean2,
        double variance1, double variance2)
{
    int size = std::min(data1.size(), data2.size());

    double c = 0;
    for(int i = 0; i < size; i++)
        c += (data1[i] - mean1) * (data2[i] - mean2);
    return c / (size * variance1 * variance2);
}

double c1ModulationCorrelation(const std::vector<double>& data1,
        const std::vector<double>& data2, double variance1, double variance2)
{
    int size = std::min(data1.size(), data2.size());

    double c = 0;
    for(int i = 0; i < size; i++)
        c += data1[i] * data2[i];
    return c / (size * variance1 * variance2);
}

std::complex<double> c2ModulationCorrelation(Signal signal1, Signal signal2,
        double variance1, double variance2)
{
    int size = std::min(signal1._signal.size(), signal2._signal.size());

    signal1.makeAnalytic();
    signal2.makeAnalytic();

    std::complex<double> c = 0;
    for(int i = 0; i < size; i++)
        c += std::conj(std::pow(signal1._signal[i], 2) /
                std::norm(signal1._signal[i])) * signal2._signal[i];
    return c / (size * variance1 * variance2);
}

double computePower(const std::vector<double>& data, double variance)
{
    int size = data.size();

    double p = 0;
    for(int i = 0; i < size; i++)
        p += data[i] * data[i];
    return p / (size * variance * variance);
}


void TextureSynthesis::computeStatistics(const Signal& signal,
        std::vector<std::complex<double>>& statistics)
{
    // Compute filtered signals
    std::vector<Signal> cochlearSignals;
    std::vector<std::vector<Signal>> modulationSignals;
    auditoryFilter(signal, cochlearSignals, modulationSignals);


    // Cochlear statistics

    int cochlearSize = cochlearSignals.size();

    std::vector<int> correlationDifs { 1, 2, 3, 5, 8, 11, 16, 21 };

    std::vector<double> cochlearMeans(cochlearSize);
    std::vector<double> cochlearVariances(cochlearSize);

    double mean, variance;
    std::vector<double> realPart;
    for(int i = 0; i < cochlearSize; i++)
    {
        realPart = cochlearSignals[i].realPart();

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
                            cochlearSignals[i - j].realPart(), mean, cochlearMeans[i - j],
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

void TextureSynthesis::auditoryFilter(const Signal& signal,
        std::vector<Signal>& cochlearSignals,
            std::vector<std::vector<Signal>>& modulationSignals)
{
    FilterBank cochlearBank, modulationBank;
    generateCochlearBank(cochlearBank);
    generateModulationBank(modulationBank);

    cochlearSignals = cochlearBank.apply(signal);
    for(Signal cochlearSignal : cochlearSignals)
    {
        cochlearSignal.makeEnvelope();
        std::vector<Signal> outSignals = modulationBank.apply(cochlearSignal);
        modulationSignals.push_back(outSignals);
    }
}

void TextureSynthesis::generateCochlearBank(FilterBank& cochlearBank)
{
    int numFilters = 32;

    Filter *filter;

    for(int i = 1; i <= numFilters; i++)
    {
        filter = (Filter*)(new CochlearFilter(erbsInverse(i)));
        cochlearBank.addFilter(std::shared_ptr<Filter>(filter));
    }
}

void TextureSynthesis::generateModulationBank(FilterBank& modulationBank)
{
    int numFilters = 20;
    double scalar = 4.112 * exp(-7);
    double offset = 0.5 - scalar;

    Filter *filter;

    for(int i = 0; i < numFilters; i++)
    {
        filter = (Filter*)(new CochlearFilter(scalar * exp(i) + offset));
        modulationBank.addFilter(std::shared_ptr<Filter>(filter));
    }
}

double TextureSynthesis::erbsInverse(int numFilters)
{
    return 676170.4 / (47.06538 - exp(0.08950404 * numFilters)) - 14678.49;
}
