#include "../util/sampling_utils.h"

namespace lightwave
{
    Distribution2D::Distribution2D(const float *func, int nu, int nv)
    {
        pConditionalV.reserve(nv);
        for (int v = 0; v < nv; ++v)
        {
            // Compute conditional sampling distribution for v tilda
            pConditionalV.emplace_back(new Distribution1D(&func[v * nu], nu));
        }
        // Compute marginal sampling distribution p tilda
        std::vector<float> marginalFunc;
        marginalFunc.reserve(nv);
        for (int v = 0; v < nv; ++v)
            marginalFunc.push_back(pConditionalV[v]->funcInt);
        pMarginal.reset(new Distribution1D(&marginalFunc[0], nv));
    }
} // namespace lightwave
