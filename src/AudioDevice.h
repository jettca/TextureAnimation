#include <SDL2/SDL.h>
#include <SDL2/SDL_mixer.h>

#include "Signal.h"

namespace TextureSynthesis
{
    class UserData
    {
    public:
        UserData(const Signal *signal, std::function<void()> callback);
        const Signal *_signalP;
        std::function<void()> _callback;
        int _signalPos;
    };

    class AudioDevice
    {
    public:
        AudioDevice();
        void play(const Signal& signal, std::function<void()> callback);

    private:
        std::shared_ptr<UserData> _curUserDataP;

        static void writeAudio(void *_, Uint8 *stream, int len);
    };
}
