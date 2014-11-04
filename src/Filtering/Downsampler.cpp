#include "Filtering/Downsampler.h"
#include "Filtering/LowpassFilter.h"

using namespace TextureSynthesis;

Downsampler::Downsampler(Signal& signal, double targetRate):
    _fineStructure(signal),
    _newSampleRate(signal.sampleRate /
        pow(2, (int)(log(signal.sampleRate / targetRate) / log(2))))
{
    int oldSize = signal.size();

    // Apply a lowpass filter to obtain the fine structure
    LowpassFilter filter(_newSampleRate / 2);
    filter.filter(signal);
    for(int i = 0; i < oldSize; i++)
        _fineStructure[i] -= signal[i];

    // Fill values into new positions
    int newSize = (int)(_newSampleRate / signal.sampleRate * oldSize);
    double sizeRatio = oldSize / newSize;
    for(int i = 0; i < newSize; i++)
        signal[i] = signal[(int)(sizeRatio * i)];

    // Set the new size and sample rate
    signal.resize(std::min(newSize, oldSize));
    signal.sampleRate = _newSampleRate;
}

void Downsampler::revert(Signal& signal) const
{
    // Resize and set the sample rate on the signal
    int revertedSize = _fineStructure.size();
    signal.resize(revertedSize);
    signal.sampleRate = _fineStructure.sampleRate;

    // Move the values to their new positions and add the fine structure
    double sizeRatio = signal.size() / revertedSize;
    for(int i = revertedSize - 1; i >= 0; i--)
        signal[i] = signal[sizeRatio * i] + _fineStructure[i];
}

double Downsampler::actualSampleRate() const
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
