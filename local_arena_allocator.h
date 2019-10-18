#ifndef LOCAL_ARENA_ALLOCATOR_H
#define LOCAL_ARENA_ALLOCATOR_H

#include <emu_c_utils/emu_c_utils.h>
#include <cassert>
#include <cstddef>
#include <memory>

namespace emu {

    class local_arena
    {
    protected:
	typedef unsigned char uchar;
	// Value originally returned from mw_mallocrepl
	void * buffer;
	// Pointer to next available byte on this nodelet
	uchar * next_chunk;
	// Size of each stripe, in elements
	size_t size;
    public:

	// @param stripe_size: Number of bytes to reserve on each nodelet
	local_arena(size_t stripe_size)
	{
            assert(stripe_size > 0);
            // Allocate a stripe on each nodelet
            void * ptr = mw_mallocrepl(stripe_size);
            assert(ptr); // FIXME
            // Save the pointer in all copies of the allocator
            // since we don't know which one will get freed
            mw_replicated_init((long*)&buffer, (long)ptr);
            mw_replicated_init((long*)&size, (long)stripe_size);

            // Let the nth copy of next_chunk point to the stripe on the nth nodelet
            for (long nlet = 0; nlet < NODELETS(); ++nlet) {
                uchar * nth_chunk = (uchar*)mw_get_nth(buffer, nlet);
                uchar ** nth_ptr = (uchar**)mw_get_nth(&next_chunk, nlet);
                *nth_ptr = nth_chunk;
            }
	}

	~local_arena()
	{
            mw_free(buffer);
	}

	void *
            allocate(size_t n, const void * hint)
	{
            if (hint) { MIGRATE((void*)hint); }
            // Serial version
            void * ptr = next_chunk;
            next_chunk += n;
            // Parallel-safe version
            // T * ptr = (T*)ATOMIC_ADDMS((long*)next_chunk, (long)(n * sizeof(T)));
            return ptr;
	}
    };

// Reserve 2GB on each nodelet for satisfying allocations
    extern replicated local_arena g_arena;

    template<typename T>
        class local_arena_allocator
    {
    protected:
	local_arena & arena;
    public:
	typedef T value_type;

        // Declare friendship with all templated versions of this class
        template<typename U>
            friend class local_arena_allocator;

        // Default constructor, uses the global replicated arena
    local_arena_allocator() : arena(g_arena) {}

        // Alternate constructor for specifying a different arena to use
	explicit
            local_arena_allocator(local_arena & arena) : arena(arena) {}

        // Templated copy-constructor, for use with rebind
        template<typename U>
            local_arena_allocator(const local_arena_allocator<U>& other) : arena(other.arena) {}

        template<typename U>
            bool operator== (const local_arena_allocator<U>& other)
        {
	    return arena == other.arena;
        }

	T *
            allocate(size_t n)
	{
            return static_cast<T*>(
                arena.allocate(n * sizeof(T), nullptr)
		);
	}

	T *
            allocate(size_t n, const void * hint)
	{
            return static_cast<T*>(
                arena.allocate(n * sizeof(T), hint)
		);
	}

	void
            deallocate(T *, size_t)
	{
            // Not implemented
	}
};

} // end namespace emu

#endif // LOCAL_ARENA_ALLOCATOR_H
