#ifndef CARDINAL_CORE_PARALLEL
#define CARDINAL_CORE_PARALLEL

#include <thread>

//// Safety
//----------
// This is a thin unsafe wrapper around a thread array
// The flag is set to false if any threads fail to join
// The caller takes responsiblity for ensuring:
// * the flag has a lifetime that outlives the struct
// * the flag is checked and handled appropriately
// Threads that fail to join will be detached

struct threads
{
	std::thread * tasks;
	int count;
	bool * ok;

	explicit threads(int n, bool * flag) : 
		tasks(new std::thread[n]), count(n), ok(flag) 
	{
		*ok = true;
	}

	~threads()
	{
		join_all();
		delete[] tasks;
	}

	threads(const threads&) = delete;

	threads(threads&&) noexcept = delete;

	threads& operator=(const threads&) = delete;

	threads& operator=(threads&&) noexcept = delete;

	void join_all() noexcept
	{
		for ( int i = 0; i < count; ++i )
		{
			if ( !tasks[i].joinable() )
				continue;
			try {
				tasks[i].join();
			} catch (...) {
				*ok = false;
				tasks[i].detach();
			}
		}
	}
};

#endif // CARDINAL_CORE_PARALLEL
