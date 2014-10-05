#include "AuditoryModel.h"
#include "CochlearFilter.h"
#include <cmath>

using namespace TextureSynthesis;

std::vector<std::vector<Signal>> TextureSynthesis::auditoryFilter(Signal signal)
{
    FilterBank cochlearBank, modulationBank;
    generateCochlearBank(cochlearBank);
    generateModulationBank(modulationBank);

    std::vector<std::vector<Signal>> outSignals;

    std::vector<Signal> cochlearSignals = cochlearBank.apply(signal);
    for(Signal cochlearSignal : cochlearSignals)
    {
        cochlearSignal.makeEnvelope();
        std::vector<Signal> modulationSignals = modulationBank.apply(cochlearSignal);
        outSignals.push_back(modulationSignals);
    }

    return outSignals;
}

void TextureSynthesis::generateCochlearBank(FilterBank& cochlearBank)
{
    int numFilters = 30;

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
    }
}

double TextureSynthesis::erbsInverse(int numFilters)
{
    return 676170.4 / (47.06538 - exp(0.08950404 * numFilters)) - 14678.49;
}
