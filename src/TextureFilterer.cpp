#include "TextureFilterer.h"

#include "LowpassFilter.h"
#include "CochlearFilter.h"

using namespace TextureSynthesis;

void TextureFilterer::auditoryFilter(const Signal& signal,
        std::vector<Signal>& cochlearSignals,
        std::vector<std::vector<Signal>>& modulationSignals)
{
    FilterBank cochlearBank, modulationBank;
    generateCochlearBank(cochlearBank);
    generateModulationBank(modulationBank);

    cochlearBank.apply(signal, cochlearSignals);

    int numCochlearSignals = cochlearSignals.size();
    if(modulationSignals.size() < numCochlearSignals)
        modulationSignals.resize(numCochlearSignals);

    Signal subbandPhase(signal._signal.size(), signal._sampleRate);
    _downsamplers.clear();
    for(int i = 0; i < numCochlearSignals; i++)
    {
        cochlearSignals[i].makeEnvelope(subbandPhase);
        _subbandPhases.push_back(subbandPhase);

        _downsamplers.push_back(Downsampler(cochlearSignals[i], 400));
        cochlearSignals[i].pow(_cochlearExponent);
        modulationBank.apply(cochlearSignals[i], _someModulationSignals);
        modulationSignals[i] = _someModulationSignals;
    }
}

void TextureFilterer::recombine(std::vector<Signal>& cochlearEnvelopes,
        Signal& combinedSignal)
{
    int signalSize = combinedSignal._signal.size();
    for(int j = 0; j < signalSize; j++)
        combinedSignal._signal[j] = 0;

    int numEnvelopes = cochlearEnvelopes.size();
    for(int i = 0; i < numEnvelopes; i++)
    {
        cochlearEnvelopes[i].pow(-_cochlearExponent);
        _downsamplers[i].revert(cochlearEnvelopes[i]);
        for(int j = 0; j < signalSize; j++)
            combinedSignal._signal[j] += std::real(cochlearEnvelopes[i]._signal[j] *
                _subbandPhases[i]._signal[j]);
    }
}

void TextureFilterer::generateCochlearBank(FilterBank& cochlearBank)
{
    int numFilters = 32;

    Filter *filter;

    for(int i = 1; i <= numFilters; i++)
    {
        filter = (Filter*)(new CochlearFilter(erbsInverse(i)));
        cochlearBank.addFilter(std::shared_ptr<Filter>(filter));
    }
}

void TextureFilterer::generateModulationBank(FilterBank& modulationBank)
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

double TextureFilterer::erbsInverse(int numFilters)
{
    return 676170.4 / (47.06538 - exp(0.08950404 * numFilters)) - 14678.49;
}

std::vector<Signal> TextureFilterer::_someModulationSignals;
std::vector<Downsampler> TextureFilterer::_downsamplers;
std::vector<Signal> TextureFilterer::_subbandPhases;
double TextureFilterer::_cochlearExponent = 0.3;
