#include "CochlearFilter.h"
#include <math.h>
#include <iostream>

using namespace TextureSynthesis;

CochlearFilter::CochlearFilter(double centerFrequency)
    : _centerFrequency(centerFrequency)
{ }

void CochlearFilter::filter(Aquila::SpectrumType& spectrum, double sampleRate)
{
    int size = spectrum.size();
    Aquila::SpectrumType filterSpectrum(size);
    double frequency;

    double f = _centerFrequency / 1000;
    double bandwidth3db = 6.23 * f * f + 93.39 * f + 28.52;
    double bandwidth = bandwidth3db * 3.0 / 2.0;
    double angularFrequency = M_PI / bandwidth;

    int cosStart = int(size * (_centerFrequency - bandwidth / 2) / sampleRate);
    int cosEnd = int(size * (_centerFrequency + bandwidth / 2) / sampleRate);
    
    // Clip start and end
    cosStart = std::max(std::min(cosStart, size), 0);
    cosEnd = std::max(std::min(cosEnd, size), 0);

    for(int i = 0; i < cosStart; i++)
        filterSpectrum[i] = 0;
    for(int i = cosStart; i < cosEnd; i++)
    {
        frequency = sampleRate / size * i;
        filterSpectrum[i] = cos(angularFrequency * (frequency - _centerFrequency));
    }
    for(int i = cosEnd; i < size; i++)
        filterSpectrum[i] = 0;

    std::transform(
        std::begin(spectrum),
        std::end(spectrum),
        std::begin(filterSpectrum),
        std::begin(spectrum),
        [] (Aquila::ComplexType x, Aquila::ComplexType y) { return x * y; }
    );
}
