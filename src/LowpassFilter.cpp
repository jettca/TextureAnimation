#include "LowpassFilter.h"
#include <iostream>

using namespace TextureSynthesis;

LowpassFilter::LowpassFilter(double lowpassFrequency)
    : _lowpassFrequency(lowpassFrequency)
{ }

void LowpassFilter::filter(Aquila::SpectrumType& spectrum, double sampleRate)
{
    int size = spectrum.size();
    Aquila::SpectrumType filterSpectrum(size);

    int numlow = size * _lowpassFrequency / sampleRate;
    for(int i = 0; i < numlow; i++)
        filterSpectrum[i] = 1;
    for(int i = numlow; i < size; i++)
        filterSpectrum[i] = 0;

    std::transform(
        std::begin(spectrum),
        std::end(spectrum),
        std::begin(filterSpectrum),
        std::begin(spectrum),
        [] (Aquila::ComplexType x, Aquila::ComplexType y) { return x * y; }
    );
}
