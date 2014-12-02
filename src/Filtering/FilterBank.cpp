#include "Filtering/FilterBank.h"

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
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(signal.size());
    Aquila::SpectrumType spectrum = fft->fft(signal.samples());
    apply(spectrum, signal.sampleRate, outSignals);
}

void FilterBank::apply(const Aquila::SpectrumType& spectrum, double sampleRate,
        std::vector<Signal>& outSignals)
{
    Signal signal(spectrum.size(), sampleRate);
    Aquila::SpectrumType filteredSpectrum;
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(spectrum.size());

    if(outSignals.size() < _filters.size())
    {
        for(int i = 0; i < _filters.size(); i++)
            outSignals.push_back(Signal(spectrum.size(), sampleRate));
    }
    
    for(int i = 0; i < _filters.size(); i++)
    {
        filteredSpectrum = spectrum;

        _filters.at(i)->filter(filteredSpectrum, sampleRate);

        fft->ifft(filteredSpectrum, signal.samples());
        outSignals[i].set(signal);
    }
}

void FilterBank::apply(Signal& signal, int index)
{
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(signal.size());
    Aquila::SpectrumType spectrum = fft->fft(signal.samples());
    _filters.at(index)->filter(spectrum, signal.sampleRate);
    fft->ifft(spectrum, signal.samples());
}

const std::shared_ptr<Filter> FilterBank::getFilter(int index) const
{
    return _filters.at(index);
}
