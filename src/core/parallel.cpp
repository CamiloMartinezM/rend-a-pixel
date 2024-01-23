#include <lightwave/parallel.hpp>

namespace lightwave
{
    thread_local int ThreadIndex = -1;

    thread_local int MaxThreads = std::thread::hardware_concurrency();
} // namespace lightwave