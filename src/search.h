#ifndef CARDINAL_SEARCH
#define CARDINAL_SEARCH

//// Comparison
//--------------

// signed absolute or relative difference
template<typename T>
double sdiff(T x, T ref, bool relative = false)
{
	if ( relative )
		return static_cast<double>(x - ref) / ref;
	else
		return static_cast<double>(x - ref);
}

#endif // CARDINAL_SEARCH
