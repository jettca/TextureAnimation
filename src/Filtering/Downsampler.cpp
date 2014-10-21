#include "Filtering/Downsampler.h"
#include "Filtering/LowpassFilter.h"

using namespace TextureSynthesis;

Downsampler::Downsampler(Signal& signal, double targetRate)
    : _fineStructure(signal),
    _newSampleRate(signal._sampleRate /
        pow(2, (int)(log(signal._sampleRate / targetRate) / log(2))))
{
    LowpassFilter filter(_newSampleRate / 2);
    filter.filter(signal);

    int oldSize = signal._signal.size();
    int newSize = (int)(_newSampleRate / signal._sampleRate * oldSize);

    for(int i = 0; i < signal._signal.size(); i++)
        _fineStructure._signal[i] -= signal._signal[i];

    double sizeRatio = signal._signal.size() / newSize;

    for(int i = 0; i < newSize; i++)
        signal._signal[i] = signal._signal.at((int)(sizeRatio * i));

    signal._signal.resize(std::min(newSize, (int)signal._signal.size()));
    signal._sampleRate = _newSampleRate;
}

void Downsampler::revert(Signal& signal)
{
    int revertedSize = _fineStructure._signal.size();
    signal._signal.resize(revertedSize);
    signal._sampleRate = _fineStructure._sampleRate;

    double sizeRatio = signal._signal.size() / revertedSize;

    for(int i = revertedSize - 1; i >= 0; i--)
        signal._signal[i] = signal._signal[sizeRatio * i] + _fineStructure._signal[i];
}

double Downsampler::actualSampleRate()
{
    return _newSampleRate;
}

std::pair<int, double> Downsampler::newSizeAndRate(int oldSize, double oldRate,
        double targetRate)
{
    double newRate = oldRate / pow(2, (int)(log(oldRate / targetRate) / log(2)));
    int newSize = (int)(newRate / oldRate * oldSize);

    return std::pair<int, double>(newSize, newRate);
}
