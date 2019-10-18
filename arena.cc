#include "local_arena_allocator.h"

// Reserve 2GB on each nodelet for satisfying allocations
replicated emu::local_arena emu::g_arena((1UL<<31));
