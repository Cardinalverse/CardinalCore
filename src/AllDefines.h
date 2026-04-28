#ifndef CARDINAL_CORE_ALLDEFINES
#define CARDINAL_CORE_ALLDEFINES

//// Infinity, NaN, and NA
//-------------------------

// return missing value<T>
template<typename T>
T my_NA();

template<> inline
int my_NA<int>()
{
	return NA_INTEGER;
}

template<> inline
double my_NA<double>()
{
	return NA_REAL;
}

// check if missing value<T>
template<typename T>
bool my_isNA(T x);

template<> inline
bool my_isNA<int>()
{
	return x == NA_INTEGER || x == NA_LOGICAL;
}

template<> inline
bool my_isNA<double>()
{
	return ISNA(x) || ISNAN(x);
}

// return finite minimum value<T>
template<typename T>
T my_MinVal();

// return finite maximum value<T>
template<typename T>
T my_MaxVal();


#endif // CARDINAL_CORE_ALLDEFINES
