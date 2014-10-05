#include "Signal.h"

#include <iostream>

using namespace TextureSynthesis;

Signal::Signal(int length, double sampleRate)
    : _signal(length), _sampleRate(_sampleRate)
{ }

Signal::Signal(Aquila::SignalSource source)
    : _signal(source.length()), _sampleRate(source.getSampleFrequency())
{
    int signalSize = _signal.size();
    for(int i = 0; i < signalSize; i++)
        _signal[i] = source.sample(i);
}

void Signal::makeAnalytic()
{
    int signalSize = _signal.size();
    std::shared_ptr<Aquila::Fft> fft = Aquila::FftFactory::getFft(signalSize);
    Aquila::SpectrumType spectrum = fft->fft(_signal);

    for(int i = (signalSize + 1) / 2; i < signalSize; i++)
        spectrum[i] = 0;

    fft->ifft(spectrum, _signal);
}

void Signal::makeEnvelope()
{
    makeAnalytic();
    int signalSize = _signal.size();
    for(int i = 0; i < signalSize; i++)
    {
        _signal[i] = std::norm(_signal[i]);
    }
}

void Signal::pow(double a)
{
    int signalSize = _signal.size();
    for(int i = 0; i < signalSize; i++)
    {
        _signal[i] = std::pow(_signal[i], a);
    }
}

std::vector<double> Signal::realPart()
{
    int signalSize = _signal.size();
    std::vector<double> realSignal(signalSize);
    for(int i = 0; i < signalSize; i++)
        realSignal[i] = std::real(_signal[i]);

    return realSignal;
}

std::vector<double> Signal::imaginaryPart()
{
    int signalSize = _signal.size();
    std::vector<double> imaginarySignal(signalSize);
    for(int i = 0; i < signalSize; i++)
        imaginarySignal[i] = std::imag(_signal[i]);

    return imaginarySignal;
}
