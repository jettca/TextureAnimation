#include "StatsGenerator.h"

#include "Synthesis/TextureFilterer.h"

using namespace TextureSynthesis;

void StatsGenerator::computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
        const std::vector<std::vector<Signal>>& modulationSignals,
        std::vector<double>& statistics)
{
    computeStatistics(cochlearEnvelopes, modulationSignals, statistics, nullptr, nullptr);
}

void StatsGenerator::computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
        const std::vector<std::vector<Signal>>& modulationSignals,
        std::vector<double>& statistics, const FilterBank& modBank,
        std::vector<std::vector<std::vector<double>>>& jacobian)
{
    computeStatistics(cochlearEnvelopes, modulationSignals, statistics, &modBank, &jacobian);
}

void StatsGenerator::computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
        const std::vector<std::vector<Signal>>& modulationSignals,
        std::vector<double>& statistics, const FilterBank* modBank,
        std::vector<std::vector<std::vector<double>>>* jacobian)
{
    // TODO: optimize vector handling

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
        statistics.insert(statistics.end(), moments.begin(), moments.end());
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
        double envVariance = cochlearVariances[i];
        double envMean = cochlearMeans[i];
        std::vector<double> envRealPart = cochlearEnvelopes[i].realPart();

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
                        powerGrad.push_back(computePowerGrad(envRealPart, modulationSignals[i][j],
                                    *(modBank->getFilter(j)), envMean, envVariance));
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
                                    modulationSignals[i][n], modulationSignals[j][n],
                                    *(modBank->getFilter(n)), modulationVariances[i][n],
                                    modulationVariances[j][n], env == i));
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
                            modulationVariances[i][n + 1], *(modBank->getFilter(n)),
                            *(modBank->getFilter(n + 1)), env == n);
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

    double varyingStdDev = sqrt(varyingVar) + 1e-10;
    double fixedStdDev = sqrt(fixedVar) + 1e-10;

    // TODO: why is one of the variances always zero?

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

std::vector<double> StatsGenerator::computePowerGrad(const std::vector<double>& data, const Signal& modSignal,
        const Filter& filter, double mean, double variance)
{
    int size = data.size();
    std::vector<double> grad(size);
    double sizeInv = 1.0 / size;

    Signal filteredModSignal(modSignal);
    double sqrMean = 0;
    for(int i = 0; i < size; i++)
    {
        sqrMean += pow(std::real(modSignal[i]), 2);
        filteredModSignal[i] *= sizeInv;
    }
    sqrMean *= sizeInv;
    filter.filter(filteredModSignal);

    variance += 1e-10;
    for(int i = 0; i < size; i++)
    {
        grad[i] = (2 * variance * std::real(filteredModSignal[i]) - 2 * sqrMean
            * (data[i] - mean) * sizeInv) / pow(variance, 2);
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

std::vector<double> StatsGenerator::c1ModulationCorrelationGrad(const Signal& modSignal1,
        const Signal& modSignal2, const Filter& filter, double variance1,
        double variance2, bool varyingData1)
{
    int size = std::min(modSignal1.size(), modSignal2.size());
    double sizeInv = 1.0 / size;
    std::vector<double> grad(size);

    const Signal *varying, *fixed;
    double varyingVar, fixedVar;
    if(varyingData1)
    {
        varying = &modSignal1; fixed = &modSignal2;
        varyingVar = variance1; fixedVar = variance2;
    }
    else
    {
        varying = &modSignal2; fixed = &modSignal1;
        varyingVar = variance2; fixedVar = variance1;
    }

    double varyingStdDev = sqrt(varyingVar) + 1e-10;
    double fixedStdDev = sqrt(fixedVar) + 1e-10;

    Signal filteredVarying(*varying);
    Signal filteredFixed(*fixed);

    filteredVarying.scale(sizeInv);
    filteredFixed.scale(sizeInv);

    filter.filter(filteredVarying);
    filter.filter(filteredFixed);

    double cbc = 0;
    for(int i = 0; i < size; i++)
        cbc += std::real((*varying)[i]) * std::real((*fixed)[i]) * sizeInv;
    cbc /= varyingStdDev * fixedStdDev;

    for(int i = 0; i < size; i++)
        grad[i] = std::real(filteredFixed[i]) / (varyingStdDev * fixedStdDev)
            - cbc / pow(varyingStdDev, 2) * std::real(filteredVarying[i]);

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
        const Signal& signal2, double variance1, double variance2, const Filter& filter1,
        const Filter& filter2, bool varyingData1)
{
    int size = std::min(signal1.size(), signal2.size());
    double sizeInv = 1.0 / size;
    std::vector<std::complex<double>> grad(size);

    const Signal *varying, *fixed;
    double varyingVar, fixedVar;
    const Filter *varyingFilter, *fixedFilter;
    if(varyingData1)
    {
        varying = &signal1; fixed = &signal2;
        varyingVar = variance1; fixedVar = variance2;
        varyingFilter = &filter1; fixedFilter = &filter2;
    }
    else
    {
        varying = &signal2; fixed = &signal1;
        varyingVar = variance2; fixedVar = variance1;
        varyingFilter = &filter2; fixedFilter = &filter1;
    }

    Signal varyingAnalytic(*varying);
    Signal fixedAnalytic(*fixed);
    varyingAnalytic.makeAnalytic();
    fixedAnalytic.makeAnalytic();

    Signal u(size, varyingAnalytic.sampleRate);
    for(int i = 0; i < size; i++)
    {
        double norm = std::norm(varyingAnalytic[i]);
        u[i] = 2 * pow(std::real(varyingAnalytic[i]), 2) / norm - norm;
    }

    std::vector<double> fa_real = fixedAnalytic.realPart();
    std::vector<double> fa_imag = fixedAnalytic.imaginaryPart();

    double va_magnitude = 0;
    double fa_magnitude = 0;
    for(int i = 0; i < size; i++)
    {
        va_magnitude += pow(sizeInv * std::real(varyingAnalytic[i]), 2);
        fa_magnitude += pow(sizeInv * fa_real[i], 2);
    }
    va_magnitude = sqrt(va_magnitude);
    fa_magnitude = sqrt(fa_magnitude);

    std::vector<double> stuff_real = computeShit(varyingAnalytic, fa_real, *varyingFilter);
    std::vector<double> stuff_imag = computeShit(varyingAnalytic, fa_imag, *varyingFilter);

    Signal stuff2(u);
    stuff2.scale(sizeInv);
    fixedFilter->filter(stuff2);
    Signal stuff2_imag(stuff2);
    stuff2_imag.makeAnalytic();
    stuff2_imag.scale(-1);

    double uReal = 0;
    double uImag = 0;
    for(int i = 0; i < size; i++)
    {
        uReal += std::real(u[i]) * std::real(fixedAnalytic[i]) * sizeInv;
        uImag += std::real(u[i]) * std::imag(fixedAnalytic[i]) * sizeInv;
    }

    varyingAnalytic.scale(sizeInv);
    fixedAnalytic.scale(sizeInv);
    varyingAnalytic.makeReal();
    fixedAnalytic.makeReal();
    varyingFilter->filter(varyingAnalytic);
    fixedFilter->filter(fixedAnalytic);

    std::complex<double> stuff3;
    for(int i = 0; i < size; i++)
    {
        grad[i] = std::complex<double>(stuff_real[i] + std::real(stuff2[i]),
                stuff_imag[i] + std::imag(stuff2_imag[i]));
        grad[i] *= va_magnitude * fa_magnitude;
        grad[i] -= std::complex<double>(uReal * fa_magnitude / va_magnitude,
                uImag * va_magnitude / fa_magnitude) * varyingAnalytic[i] * fixedAnalytic[i];
        grad[i] /= pow(va_magnitude, 2) + pow(fa_magnitude, 2);
    }

    return grad;
}


std::vector<double> StatsGenerator::computeShit(const Signal& varyingAnalytic,
        const std::vector<double>& fa, const Filter& varyingFilter)
{
    int size = varyingAnalytic.size();
    double sizeInv = 1.0 / size;
    std::vector<double> stuff(size);

    Signal toFilter1(size, varyingAnalytic.sampleRate);
    Signal toFilter2(size, varyingAnalytic.sampleRate);
    Signal toFilter3(size, varyingAnalytic.sampleRate);
    Signal toFilter4(size, varyingAnalytic.sampleRate);
    Signal toFilter5(size, varyingAnalytic.sampleRate);

    double va_real, va_imag, va_norm;
    for(int i = 0; i < size; i++)
    {
        va_real = std::real(varyingAnalytic[i]);
        va_imag = std::imag(varyingAnalytic[i]);
        va_norm = std::norm(varyingAnalytic[i]);

        toFilter1[i] = 4 * va_real / va_norm * sizeInv * fa[i];
        toFilter2[i] = 2 * pow(va_real, 2) / pow(va_norm, 3) * va_real * sizeInv * fa[i];
        toFilter3[i] = toFilter2[i] / va_real * va_imag;
        toFilter4[i] = toFilter1[i] * 0.25;
        toFilter5[i] = va_imag / va_norm * sizeInv * fa[i];
    }

    varyingFilter.filter(toFilter1);
    varyingFilter.filter(toFilter2);
    varyingFilter.filter(toFilter3);
    toFilter3.makeAnalytic();
    varyingFilter.filter(toFilter4);
    varyingFilter.filter(toFilter5);
    toFilter5.makeAnalytic();

    for(int i = 0; i < size; i++)
        stuff[i] = std::real(toFilter1[i]) - std::real(toFilter2[i]) + std::imag(toFilter3[i])
            - std::real(toFilter4[i]) + std::imag(toFilter5[i]);

    return stuff;
}
