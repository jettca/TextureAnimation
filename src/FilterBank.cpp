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

void FilterBank::apply(const Signal& signal, std::vector<Signal>& outSignals)
{
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(signal._signal.size());

    Aquila::SpectrumType spectrum = fft->fft(signal._signal);

    apply(spectrum, signal._sampleRate, outSignals);
}

void FilterBank::apply(const Aquila::SpectrumType& spectrum, double sampleRate,
        std::vector<Signal>& outSignals)
{
    Signal signal(spectrum.size(), sampleRate);
    Aquila::SpectrumType filteredSpectrum;
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(spectrum.size());

    if(outSignals.size() < _filters.size())
        for(int i = 0; i < _filters.size(); i++)
            outSignals.push_back(Signal(spectrum.size(), sampleRate));
    
    for(int i = 0; i < _filters.size(); i++)
    {
        filteredSpectrum = spectrum;

        _filters.at(i)->filter(filteredSpectrum, sampleRate);

        fft->ifft(filteredSpectrum, signal._signal);
        outSignals[i].set(signal);
    }
}
