#include "Synthesis/TextureFilterer.h"

#include "Filtering/LowpassFilter.h"
#include "Filtering/CochlearFilter.h"

using namespace TextureSynthesis;

double TextureFilterer::targetDownsampleRate = 400;
int TextureFilterer::numCochlearEnvelopes = 32;
int TextureFilterer::numModulationSignals = 20;

TextureFilterer::TextureFilterer() :
    _downsamplers(),
    _subbandPhases()
{ }

void TextureFilterer::auditoryFilter(const Signal& signal,
        std::vector<Signal>& cochlearEnvelopes,
        std::vector<std::vector<Signal>>& modulationSignals,
        bool makeRecombinable)
{
    FilterBank cochlearBank;
    generateCochlearBank(cochlearBank);

    cochlearBank.apply(signal, cochlearEnvelopes);

    int numCochlearEnvelopes = cochlearEnvelopes.size();

    Signal subbandPhase(signal._signal.size(), signal._sampleRate);
    if(makeRecombinable)
    {
        _downsamplers.clear();
        _subbandPhases.clear();
    }
    for(int i = 0; i < numCochlearEnvelopes; i++)
    {
        cochlearEnvelopes[i].makeEnvelope(subbandPhase);
        if(makeRecombinable)
        {
            _subbandPhases.push_back(subbandPhase);
            _downsamplers.push_back(Downsampler(cochlearEnvelopes[i], targetDownsampleRate));
        }
        else
            Downsampler(cochlearEnvelopes[i], targetDownsampleRate);
        cochlearEnvelopes[i].pow(_cochlearExponent);
    }

    modulationFilter(cochlearEnvelopes, modulationSignals);
}

void TextureFilterer::modulationFilter(const std::vector<Signal>& cochlearEnvelopes,
        std::vector<std::vector<Signal>>& modulationSignals)
{
    FilterBank modulationBank;
    generateModulationBank(modulationBank);
    int numCochlearEnvelopes = cochlearEnvelopes.size();
    if(modulationSignals.size() < numCochlearEnvelopes)
        modulationSignals.resize(numCochlearEnvelopes);
    for(int i = 0; i < numCochlearEnvelopes; i++)
    {
        modulationBank.apply(cochlearEnvelopes[i], modulationSignals[i]);
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
    Filter *filter;

    for(int i = 1; i <= numCochlearEnvelopes; i++)
    {
        filter = (Filter*)(new CochlearFilter(erbsInverse(i)));
        cochlearBank.addFilter(std::shared_ptr<Filter>(filter));
    }
}

void TextureFilterer::generateModulationBank(FilterBank& modulationBank)
{
    double scalar = 4.112 * exp(-7);
    double offset = 0.5 - scalar;

    Filter *filter;

    for(int i = 0; i < numModulationSignals; i++)
    {
        filter = (Filter*)(new CochlearFilter(scalar * exp(i) + offset));
        modulationBank.addFilter(std::shared_ptr<Filter>(filter));
    }
}

double TextureFilterer::erbsInverse(int numFilters)
{
    return 676170.4 / (47.06538 - exp(0.08950404 * numFilters)) - 14678.49;
}
