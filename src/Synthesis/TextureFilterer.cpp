#include "Synthesis/TextureFilterer.h"

#include "Filtering/LowpassFilter.h"
#include "Filtering/CochlearFilter.h"

using namespace TextureSynthesis;

TextureFilterer::TextureFilterer() :
    _downsamplers(),
    _subbandPhases()
{ }

void TextureFilterer::auditoryFilter(const Signal& signal,
        std::vector<Signal>& cochlearEnvelopes,
        std::vector<std::vector<Signal>>& modulationSignals,
        bool makeRecombinable)
{
    // Generate and apply cochlear filter bank
    _cochlearBank.apply(signal, cochlearEnvelopes);

    // Initialize data for recombining envelopes
    Signal subbandPhase(signal.size(), signal.sampleRate);
    if(makeRecombinable)
    {
        _downsamplers.clear();
        _subbandPhases.clear();
    }
    
    // Convert filtered signals into downsampled envelopes
    for(int i = 0; i < _numCochlearEnvelopes; i++)
    {
        if(makeRecombinable)
        {
            cochlearEnvelopes[i].makeEnvelope(subbandPhase);
            _subbandPhases.push_back(subbandPhase);
            _downsamplers.push_back(Downsampler(cochlearEnvelopes[i], _targetDownsampleRate));
        }
        else
        {
            cochlearEnvelopes[i].makeEnvelope();
            Downsampler(cochlearEnvelopes[i], _targetDownsampleRate);
        }

        cochlearEnvelopes[i].pow(_cochlearExponent);
    }

    // Generate modulation signals
    modulationFilter(cochlearEnvelopes, modulationSignals);
}

void TextureFilterer::modulationFilter(const std::vector<Signal>& cochlearEnvelopes,
        std::vector<std::vector<Signal>>& modulationSignals) const
{
    // Resize output list if necessary
    if(modulationSignals.size() < _numCochlearEnvelopes)
        modulationSignals.resize(_numCochlearEnvelopes);

    // Apply filter bank to all envelopes
    for(int i = 0; i < _numCochlearEnvelopes; i++)
        _modulationBank.apply(cochlearEnvelopes[i], modulationSignals[i]);
}

void TextureFilterer::recombine(std::vector<Signal> cochlearEnvelopes,
        Signal& combinedSignal) const
{
    // Initialize combined signal to zeros
    int signalSize = combinedSignal.size();
    for(int j = 0; j < signalSize; j++)
        combinedSignal[j] = 0;

    // Invert power, upsample, and add envelopes and their phases to the combined signal
    for(int i = 0; i < _numCochlearEnvelopes; i++)
    {
        cochlearEnvelopes[i].pow(_oneOverCochlearExponent);
        _downsamplers[i].revert(cochlearEnvelopes[i]);
        for(int j = 0; j < signalSize; j++)
            cochlearEnvelopes[i][j] *= _subbandPhases[i][j];

        // Reapply filter and sum to the combined signal
        _cochlearBank.apply(cochlearEnvelopes[i], i);
        for(int j = 0; j < signalSize; j++)
            combinedSignal[j] += cochlearEnvelopes[i][j];
    }
}

FilterBank TextureFilterer::generateCochlearBank()
{
    FilterBank cochlearBank;

    std::shared_ptr<Filter> filter;
    for(int i = 1; i <= _numCochlearEnvelopes; i++)
    {
        filter = std::shared_ptr<Filter>((Filter*)(new CochlearFilter(erbsInverse(i))));
        cochlearBank.addFilter(std::shared_ptr<Filter>(filter));
    }

    return cochlearBank;
}

FilterBank TextureFilterer::generateModulationBank()
{
    FilterBank modulationBank;

    Filter *filter;

    for(int i = 0; i < _numModulationSignals; i++)
    {
        // TODO: make this actually correct
        filter = (Filter*)(new CochlearFilter(pow(2, i + 1)));
        modulationBank.addFilter(std::shared_ptr<Filter>(filter));
    }

    return modulationBank;
}

double TextureFilterer::erbsInverse(int numFilters)
{
    return 676170.4 / (47.06538 - exp(0.08950404 * numFilters)) - 14678.49;
}

int TextureFilterer::numCochlearEnvelopes()
{
    return _numCochlearEnvelopes;
}

int TextureFilterer::numModulationSignals()
{
    return _numModulationSignals;
}

double TextureFilterer::targetDownsampleRate()
{
    return _targetDownsampleRate;
}

const FilterBank& TextureFilterer::modulationBank()
{
    return _modulationBank;
}

FilterBank TextureFilterer::_cochlearBank = generateCochlearBank();
FilterBank TextureFilterer::_modulationBank = generateModulationBank();

int     TextureFilterer::_numCochlearEnvelopes      = 32;
int     TextureFilterer::_numModulationSignals      = 7;
double  TextureFilterer::_targetDownsampleRate      = 400;
double  TextureFilterer::_cochlearExponent          = -0.3;
double  TextureFilterer::_oneOverCochlearExponent   = 1 / _cochlearExponent;
