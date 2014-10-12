#pragma once

#include <vector>
#include "Signal.h"

namespace TextureSynthesis
{
    class StatsGenerator
    {
    public:
        static void computeStatistics(const std::vector<Signal>& cochlearSignals,
                const std::vector<std::vector<Signal>>& modulationSignals,
                std::vector<std::complex<double>>& statistics);

        static void computeStatistics(const Signal& signal,
                std::vector<std::complex<double>>& statistics);
    private:
        static double computeMean(const std::vector<double>& data);

        static double computeVariance(const std::vector<double>& data, double mean);

        static double nthMoment(const std::vector<double>& data,
                int n, double mean, double variance);

        static double crossCorrelation(const std::vector<double>& data1,
                const std::vector<double>& data2, double mean1, double mean2,
                double variance1, double variance2);

        static double c1ModulationCorrelation(const std::vector<double>& data1,
                const std::vector<double>& data2, double variance1, double variance2);

        static std::complex<double> c2ModulationCorrelation(Signal signal1, Signal signal2,
                double variance1, double variance2);

        static double computePower(const std::vector<double>& data, double variance);
    };
}
