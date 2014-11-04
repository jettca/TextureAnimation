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

// Handles playback loop
static bool stillPlaying = true;
void donePlayingCallback()
{
    stillPlaying = false;
}

int main(int argc, char **argv)
{
    // TODO: fix bad sound ??
    feenableexcept(FE_INVALID | FE_OVERFLOW);   // enable floating point exceptions

    // Handle input args
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

    // Load input file into Signal object with power of 2 length
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
        sourceSignal[i] = channelAvg / numChannels;
    }

    // Synthesize new signal
    TextureSynthesizer synthesizer(sourceSignal);
    Signal outSignal(sourceLen, sampleRate);
    synthesizer.synthesize(outSignal);

    // Clip new signal for playback
    for(int i = 0; i < outSignal.size(); i++)
    {
        if(std::real(outSignal[i]) >= 1)
            outSignal[i] = maxValue;
        else if(std::real(outSignal[i]) <= -1)
            outSignal[i] = -maxValue;
        else
            outSignal[i] *= maxValue;
    }

    // Play back or write file
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
