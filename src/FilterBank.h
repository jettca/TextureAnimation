#include "aquila/aquila.h"
#include <vector>

namespace TextureSynthesis
{
    typedef std::vector<Aquila::SampleType> Wave;

    class FilterBank
    {
    public:
        FilterBank();

        void addFilter(Aquila::SpectrumType filter);
        std::vector<Wave> apply(Wave wave);

    private:
        std::vector<Aquila::SpectrumType> filters;
    };
}
