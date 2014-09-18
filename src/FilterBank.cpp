#include "FilterBank.h"

using namespace TextureSynthesis;

FilterBank::FilterBank()
    : filters()
{ }

void FilterBank::addFilter(Aquila::SpectrumType filter)
{
    filters.push_back(filter);
}

std::vector<Wave> FilterBank::apply(Wave wave)
{
    std::vector<Wave> waves;
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(wave.size());    
    Aquila::SpectrumType spectrum = fft->fft(&wave[0]);
    
    Aquila::SpectrumType filteredSpectrum;
    Wave filteredWave(wave.size());

    for(Aquila::SpectrumType& filter : filters)
    {
        filteredSpectrum = spectrum;
        int size = filteredSpectrum.size();
        for(int i = 0; i < size; i++)
        {
            filteredSpectrum[i] *= filter[i];
        }
        
        fft->ifft(filteredSpectrum, &filteredWave[0]);
        waves.push_back(filteredWave);
    }

    return waves;
}
