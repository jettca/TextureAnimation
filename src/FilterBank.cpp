#include "FilterBank.h"

using namespace TextureSynthesis;

FilterBank::FilterBank()
    : _filters()
{ }

void FilterBank::addFilter(std::shared_ptr<Filter> filter)
{
    _filters.push_back(filter);
}

std::vector<Signal> FilterBank::apply(Aquila::SignalSource signal)
{
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(signal.length());
    Aquila::SpectrumType spectrum = fft->fft(signal.toArray());
    return apply(spectrum, signal.getSampleFrequency());
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

        fft->ifft(filteredSpectrum, &signal._samples[0]);
        signals.push_back(signal);
    }

    return signals;
}
