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
                std::vector<std::complex<double>>& statistics);

        void computePartials(const std::vector<Signal>& cochlearEnvelopes,
                const std::vector<std::vector<Signal>>& modulationSignals,
                int curEnvelope, int curSample,
                std::vector<std::complex<double>>& partials);

    private:
        double computeMean(const std::vector<double>& data);

        double computeVariance(const std::vector<double>& data, double mean);

        double nthMoment(const std::vector<double>& data,
                int n, double mean, double variance);

        double crossCorrelation(const std::vector<double>& data1,
                const std::vector<double>& data2, double mean1, double mean2,
                double variance1, double variance2);

        double c1ModulationCorrelation(const std::vector<double>& data1,
                const std::vector<double>& data2, double variance1, double variance2);

        std::complex<double> c2ModulationCorrelation(const Signal& signal1,
                const Signal& signal2, double variance1, double variance2);

        double computePower(const std::vector<double>& data, double variance);
    };
}
