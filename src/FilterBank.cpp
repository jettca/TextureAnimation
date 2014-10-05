#include "FilterBank.h"

#include <complex>

using namespace TextureSynthesis;

FilterBank::FilterBank()
    : _filters()
{ }

void FilterBank::addFilter(std::shared_ptr<Filter> filter)
{
    _filters.push_back(filter);
}

std::vector<Signal> FilterBank::apply(Signal signal)
{
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(signal._signal.size());
    Aquila::SpectrumType spectrum = fft->fft(signal._signal);
    return apply(spectrum, signal._sampleRate);
}

std::vector<Signal> FilterBank::apply(Aquila::SpectrumType spectrum, double sampleRate)
{
    std::vector<Signal> signals;

    Signal signal(spectrum.size(), sampleRate);
    Aquila::SpectrumType filteredSpectrum;
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(spectrum.size());
    
    for(int i = 0; i < _filters.size(); i++)
    {
        filteredSpectrum = Aquila::SpectrumType(spectrum);
        _filters.at(i)->filter(filteredSpectrum, sampleRate);

        fft->ifft(filteredSpectrum, signal._signal);
        signals.push_back(signal);
    }

    return signals;
}
