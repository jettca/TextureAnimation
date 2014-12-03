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
                std::vector<double>& statistics, const FilterBank& modBank,
                std::vector<std::vector<std::vector<double>>>& jacobian);

        void computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
                const std::vector<std::vector<Signal>>& modulationSignals,
                std::vector<double>& statistics, const FilterBank* modBank,
                std::vector<std::vector<std::vector<double>>>* jacobian);

    private:

        void computeStatistics(const std::vector<Signal>& cochlearEnvelopes,
                const std::vector<std::vector<Signal>>& modulationSignals,
                std::vector<std::complex<double>>& statistics,
                std::vector<std::vector<std::vector<std::complex<double>>>>* jacobian);
        
        double computeMean(const std::vector<double>& data);
        std::vector<double> computeMeanGrad(const std::vector<double>& data);

        double computeVariance(const std::vector<double>& data, double mean);

        double nthMoment(const std::vector<double>& data,
                int n, double mean, double variance);
        std::vector<double> nthMomentGrad(const std::vector<double>& data,
                int n, double variance, const std::vector<double>& moments);

        double crossCorrelation(const std::vector<double>& data1,
                const std::vector<double>& data2, double mean1, double mean2,
                double variance1, double variance2);
        std::vector<double> crossCorrelationGrad(const std::vector<double>& data1,
                const std::vector<double>& data2, double mean1, double mean2,
                double variance1, double variance2, bool varyingData1);

        double computePower(const std::vector<double>& data, double variance);
        std::vector<double> computePowerGrad(const std::vector<double>& data,
                const Signal& modSignal, const Filter& filter, double mean, double variance);

        double c1ModulationCorrelation(const std::vector<double>& data1,
                const std::vector<double>& data2, double variance1, double variance2);
        std::vector<double> c1ModulationCorrelationGrad(const Signal& modSignal1,
                const Signal& modSignal2, const Filter& filter, double variance1,
                double variance2, bool varyingData1);

        std::complex<double> c2ModulationCorrelation(const Signal& signal1,
                const Signal& signal2, double variance1, double variance2);
        std::vector<std::complex<double>> c2ModulationCorrelationGrad(const Signal& signal1,
                const Signal& signal2, double variance1, double variance2, const Filter& filter1,
                const Filter& filter2, bool varyingData1);

        std::vector<double> computeShit(const Signal& varyingAnalytic,
                const std::vector<double>& fa, const Filter& varyingFilter);
    };
}
