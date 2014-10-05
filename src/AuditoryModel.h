#pragma once

#include "FilterBank.h"

namespace TextureSynthesis
{
    std::vector<std::vector<Signal>> auditoryFilter(Signal signal);
    void generateCochlearBank(FilterBank& cochlearBank);
    void generateModulationBank(FilterBank& modulationBank);
    double erbsInverse(int numFilters);
}
