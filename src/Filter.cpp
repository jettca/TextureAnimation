#include "Filter.h"

using namespace TextureSynthesis;

void Filter::filter(const Signal& signal, Aquila::SpectrumType& spectrum)
{
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(signal._signal.size());

    spectrum = fft->fft(signal._signal);
    filter(spectrum, signal._sampleRate);
}

void Filter::filter(Signal& signal)
{
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(signal._signal.size());
    Aquila::SpectrumType spectrum = fft->fft(signal._signal);
    filter(spectrum, signal._sampleRate);
    fft->ifft(spectrum, signal._signal);
}
