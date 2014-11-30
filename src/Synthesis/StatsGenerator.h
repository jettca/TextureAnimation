#pragma once

#include <vector>
#include "Filtering/Signal.h"
#include "Synthesis/TextureFilterer.h"

namespace TextureSynthesis
{
    class StatsGenerator
    {
    public:
        void computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
                const std::vector<std::vector<Signal>>& modulationSignals,
                std::vector<double>& statistics);

        void computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
                const std::vector<std::vector<Signal>>& modulationSignals,
                std::vector<double>& statistics,
                std::vector<std::vector<std::vector<std::complex<double>>>> jacobian);

    private:

        void computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
                const std::vector<std::vector<Signal>>& modulationSignals,
                std::vector<std::complex<double>>& statistics,
                std::vector<std::vector<std::vector<std::complex<double>>>>* jacobian);
        
        double computeMean(const std::vector<double>& data);
        std::vector<double> computeMeanGrads(const std::vector<double>& data);

        double computeVariance(const std::vector<double>& data, double mean);
        std::vector<double> computeVarianceGrads(const std::vector<double>& data, double mean);

        double nthMoment(const std::vector<double>& data,
                int n, double mean, double variance);
        std::vector<double> nthMomentGrads(const std::vector<double>& data,
                int n, double mean, double variance);

        double crossCorrelation(const std::vector<double>& data1,
                const std::vector<double>& data2, double mean1, double mean2,
                double variance1, double variance2);
        std::vector<double> crossCorrelationGrads(const std::vector<double>& data1,
                const std::vector<double>& data2, double mean1, double mean2,
                double variance1, double variance2);

        double c1ModulationCorrelation(const std::vector<double>& data1,
                const std::vector<double>& data2, double variance1, double variance2);
        std::vector<double> c1ModulationCorrelationGrads(const std::vector<double>& data1,
                const std::vector<double>& data2, double variance1, double variance2);

        std::complex<double> c2ModulationCorrelation(const Signal& signal1,
                const Signal& signal2, double variance1, double variance2);
        std::vector<std::complex<double>> c2ModulationCorrelationGrads(const Signal& signal1,
                const Signal& signal2, double variance1, double variance2);

        double computePower(const std::vector<double>& data, double variance);
        std::vector<double> computePowerGrads(const std::vector<double>& data, double variance);
    };
}
