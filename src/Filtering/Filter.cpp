#include "Filtering/Filter.h"

using namespace TextureSynthesis;

void Filter::filter(const Signal& signal, Aquila::SpectrumType& spectrum)
{
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(signal.size());

    spectrum = fft->fft(signal.samples());
    filter(spectrum, signal.sampleRate);
}

void Filter::filter(Signal& signal)
{
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(signal.size());
    Aquila::SpectrumType spectrum = fft->fft(signal.samples());
    filter(spectrum, signal.sampleRate);
    fft->ifft(spectrum, signal.samples());
}
