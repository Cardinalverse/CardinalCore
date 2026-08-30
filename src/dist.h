#ifndef CARDINAL_CORE_DIST
#define CARDINAL_CORE_DIST

#include <memory>
#include "core.h"
#include "kernels.h"

//// Distance
//-----------
// Search utilities

// P-norms
enum Norm {
	L1,   // L1 norm:  sum(|x_i|)   -> Manhattan distance
	L2,   // L2 norm:  sum(|x_i|^2) -> Euclidean distance
	LInf, // Max norm: max(|x_i|)   -> Maximum distance
};

// Get the p-norm of a vector
template<Norm P, Vec V>
auto norm(V x) noexcept
{
	if constexpr ( P == L1 )
		return sum(abs(x));
	else if constexpr ( P == L2 )
		return pow(sum(x * x), 0.5);
	else if constexpr ( P == LInf )
		return max(abs(x));
	else
		static_assert(dependent_false<V>, "unsupported norm");
}

// Get a Minkowski distance between two vectors
template<Num T = double, Vec L, Vec R>
T dist(L lhs, R rhs, Norm p = L2) noexcept
{
	auto d = coerce<T>(lhs) - coerce<T>(rhs);
	switch(p) {
		case L1: return norm<L1>(d);
		case L2: return norm<L2>(d);
		case LInf: return norm<LInf>(d);
	}
}

#endif // CARDINAL_CORE_DIST
