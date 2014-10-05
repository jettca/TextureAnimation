#pragma once

#include "FilterBank.h"

namespace TextureSynthesis
{
    void computeStatistics(const Signal& signal,
            std::vector<std::complex<double>>& statistics);
    void auditoryFilter(const Signal& signal, std::vector<Signal>& cochlearSignals,
                std::vector<std::vector<Signal>>& modulationSignals);
    void generateCochlearBank(FilterBank& cochlearBank);
    void generateModulationBank(FilterBank& modulationBank);
    double erbsInverse(int numFilters);
}
