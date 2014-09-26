#include "Filter.h"

using namespace TextureSynthesis;

Aquila::SpectrumType Filter::filter(Aquila::SignalSource signal)
{
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(signal.length());

    Aquila::SpectrumType spectrum = fft->fft(signal.toArray());
    filter(spectrum, signal.getSampleFrequency());

    return spectrum;
}
