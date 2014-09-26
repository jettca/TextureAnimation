#include "aquila/aquila.h"
#include "FilterBank.h"
#include "LowpassFilter.h"
#include "CochlearFilter.h"
#include "AudioDevice.h"

using namespace TextureSynthesis;

static bool _donePlaying = false;

void donePlayingCallback()
{
    _donePlaying = true;
}

int main()
{
    const std::size_t SIZE = pow(2, 15);
    const Aquila::FrequencyType sampleFreq = 22050, f1 = 96, f2 = 813;

    Aquila::SineGenerator sine1(sampleFreq);
    sine1.setAmplitude(32).setFrequency(f1).generate(SIZE);
    Aquila::SineGenerator sine2(sampleFreq);
    sine2.setAmplitude(8).setFrequency(f2).setPhase(0.75).generate(SIZE);
    Aquila::SignalSource sum = sine1 + sine2;

    CochlearFilter *filter1P = new CochlearFilter(96);
    LowpassFilter *filter2P = new LowpassFilter(500);

    FilterBank filterBank;
    filterBank.addFilter(std::shared_ptr<Filter>((Filter*)filter1P));
    filterBank.addFilter(std::shared_ptr<Filter>((Filter*)filter2P));

    std::vector<Signal> signals = filterBank.apply(sum);

    AudioDevice audioDevice;
    audioDevice.play(signals[0], donePlayingCallback);

    while(!_donePlaying)
        SDL_Delay(100);

    SDL_CloseAudio();

    return 0;
} 
