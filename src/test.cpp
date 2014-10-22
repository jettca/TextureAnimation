#include "aquila/aquila.h"
#include "Synthesis/AudioDevice.h"
#include "Synthesis/TextureSynthesizer.h"
#include "Filtering/Downsampler.h"

#include <iostream>
#include <cmath>
#include <thread>
#include <chrono>
#include <fenv.h>

using namespace TextureSynthesis;

static bool stillPlaying = true;

void donePlayingCallback()
{
    stillPlaying = false;
}

int main(int argc, char **argv)
{
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    bool makeWaveFile = false;
    std::string outfile;
    if(argc < 2 || argc > 3)
    {
        std::cerr << "Usage: TextureAnimation <input file> [<output file>]\n";
        return 1;
    }
    else if(argc == 3)
    {
        makeWaveFile = true;
        outfile = argv[2];
    }

    Aquila::WaveFile input(argv[1]);
    int sourceLen = std::min(pow(2, (int)(log(input.getSamplesCount()) / log(2))), pow(2, 12));
    double sampleRate = input.getSampleFrequency();
    
    int numChannels = input.getChannelsNum();
    int maxValue = pow(2, input.getBitsPerSample() - 1) - 1;
    double channelAvg;
    Signal sourceSignal(sourceLen, sampleRate);
    for(int i = 0; i < sourceLen; i++)
    {
        channelAvg = 0;
        for(int j = 0; j < numChannels; j++)
            channelAvg += input.sample(i) / maxValue;
        sourceSignal._signal[i] = channelAvg / numChannels;
    }

    TextureSynthesizer synthesizer(sourceSignal);
    Signal outSignal(sourceLen, sampleRate);
    synthesizer.synthesize(outSignal);

    for(int i = 0; i < outSignal._signal.size(); i++)
    {
        if(std::real(outSignal._signal[i]) >= 1)
            outSignal._signal[i] = maxValue;
        else if(std::real(outSignal._signal[i]) <= -1)
            outSignal._signal[i] = -maxValue;
        else
            outSignal._signal[i] *= maxValue;
    }

    if(makeWaveFile)
    {
        Aquila::SignalSource output(outSignal.realPart(), sampleRate);
        Aquila::WaveFile::save(output, outfile);
    }
    else
    {
        AudioDevice audioDevice;
        audioDevice.play(outSignal, donePlayingCallback);

        while(stillPlaying)
            std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }

    return 0;
}
