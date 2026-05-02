#ifndef CARDINAL_CORE_UTILS
#define CARDINAL_CORE_UTILS

//// Utility
//-----------

template<typename T>
void fill_buffer(
	T * buffer,      // buffer to fill
	size_t size,     // size of buffer
	T init = 0,      // value to use to initialize
	T increment = 0, // increment init for each item
	int stride = 1)
{
	for ( size_t i = 0; i < size; i++ )
	{
		buffer[stride * i] = init;
		init += increment;
	}
}

//// Structs
//-----------

template<typename T>
struct matrix {
	const T * data;
	size_t nrows;
	size_t ncols;
	ptrdiff_t row_stride;
	ptrdiff_t col_stride;
};

#endif // CARDINAL_CORE_UTILS
