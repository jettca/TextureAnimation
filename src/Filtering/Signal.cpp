#include "Signal.h"

#include <iostream>

using namespace TextureSynthesis;

Signal::Signal(int length, double sampleRate)
    : _signal(length), _sampleRate(sampleRate)
{ }

Signal::Signal(const Aquila::SignalSource& source)
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
        _signal[i] = std::norm(_signal[i]);
}

void Signal::makeEnvelope(Signal& phase)
{
    makeAnalytic();

    int signalSize = _signal.size();

    if(phase._signal.size() != signalSize)
        phase._signal.resize(signalSize);

    phase._sampleRate = _sampleRate;

    std::complex<double> envVal;
    std::complex<double> zero(0, 0);
    for(int i = 0; i < signalSize; i++)
    {
        envVal = std::norm(_signal[i]);
        if(envVal == zero)
            phase._signal[i] = 0;
        else
            phase._signal[i] = _signal[i] / envVal;
        _signal[i] = envVal;
    }
}

void Signal::pow(double a)
{
    int signalSize = _signal.size();
    for(int i = 0; i < signalSize; i++)
    {
        if(std::real(_signal[i]) != 0 && std::imag(_signal[i]) != 0)
            _signal[i] = std::pow(_signal[i], a);
    }
}

std::vector<double> Signal::realPart() const
{
    int signalSize = _signal.size();
    std::vector<double> realSignal(signalSize);
    for(int i = 0; i < signalSize; i++)
        realSignal[i] = std::real(_signal[i]);

    return realSignal;
}

std::vector<double> Signal::imaginaryPart() const
{
    int signalSize = _signal.size();
    std::vector<double> imaginarySignal(signalSize);
    for(int i = 0; i < signalSize; i++)
        imaginarySignal[i] = std::imag(_signal[i]);

    return imaginarySignal;
}

void Signal::set(const Signal& signal)
{
    _sampleRate = signal._sampleRate;

    int oldSize = _signal.size();
    int newSize = signal._signal.size();

    if(newSize != oldSize)
        _signal.resize(newSize);

    for(int i = 0; i < newSize; i++)
        _signal[i] = signal._signal[i];
}
