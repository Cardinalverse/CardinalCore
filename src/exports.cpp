
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include "order.h"
#include "search.h"
#include "kernels.h"

extern "C" {

//// Quicksort and Quickselect
//----------------------------

SEXP do_qdiff(SEXP x, SEXP ref, SEXP relative)
{
	if ( TYPEOF(x) != TYPEOF(ref) )
		Rf_error("'x' and 'ref' must have the same data type");
	if ( LENGTH(x) != LENGTH(ref) )
		Rf_error("'x' and 'ref' must have the same length");
	SEXP dx = PROTECT(Rf_allocVector(REALSXP, LENGTH(x)));
	for ( R_len_t i = 0; i < LENGTH(x); ++i )
	{
		switch(TYPEOF(x))
		{
			case INTSXP:
				REAL(dx)[i] = diff(
					INTEGER_ELT(x, i),
					INTEGER_ELT(ref, i),
					Rf_asLogical(relative));
				break;
			case REALSXP:
				REAL(dx)[i] = diff(
					REAL_ELT(x, i),
					REAL_ELT(ref, i),
					Rf_asLogical(relative));
				break;
			default:
				Rf_error("'x' and 'ref' must be integer or double");
		}
	}
	UNPROTECT(1);
	return dx;
}

SEXP do_qorder(SEXP x)
{
	SEXP index = PROTECT(Rf_allocVector(INTSXP, LENGTH(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			qsort_index(
				vec<int>::from(index).seqfill(0),
				vec<int>::from(x));
			break;
		case REALSXP:
			qsort_index(
				vec<int>::from(index).seqfill(0),
				vec<double>::from(x));
			break;
		default:
			Rf_error("'x' must be integer or double");
	}
	UNPROTECT(1);
	return index;
}

SEXP do_qselect(SEXP x, SEXP k)
{
	SEXP order = PROTECT(Rf_allocVector(TYPEOF(x), LENGTH(k)));
	SEXP index = PROTECT(Rf_allocVector(INTSXP, LENGTH(x)));
	for ( R_len_t i = 0; i < LENGTH(k); ++i )
	{
		switch(TYPEOF(x))
		{
			case INTSXP:
				INTEGER(order)[i] = qselect_index(
					vec<int>::from(index).seqfill(0),
					vec<int>::from(x),
					INTEGER_ELT(k, i));
				break;
			case REALSXP:
				REAL(order)[i] = qselect_index(
					vec<int>::from(index).seqfill(0),
					vec<double>::from(x),
					INTEGER_ELT(k, i));
				break;
			default:
				Rf_error("'x' must be integer or double");
		}
	}
	UNPROTECT(2);
	return order;
}

SEXP do_qmedian(SEXP x)
{
	double median = NA_REAL;
	switch(TYPEOF(x))
	{
		case INTSXP:
			median = qmedian(vec<int>::from(x));
			break;
		case REALSXP:
			median = qmedian(vec<double>::from(x));
			break;
		default:
			Rf_error("'x' must be integer or double");
	}
	return Rf_ScalarReal(median);
}

SEXP do_qmad(SEXP x, SEXP center, SEXP constant)
{
	double mad = NA_REAL;
	switch(TYPEOF(x))
	{
		case INTSXP:
			mad = qmad(
				vec<int>::from(x),
				Rf_asReal(center),
				Rf_asReal(constant));
			break;
		case REALSXP:
			mad = qmad(
				vec<double>::from(x),
				Rf_asReal(center),
				Rf_asReal(constant));
			break;
		default:
			Rf_error("'x' must be integer or double");
	}
	return Rf_ScalarReal(mad);
}

//// Search and nearest neighbors
//-------------------------------

SEXP do_bsearch(
	SEXP query, 
	SEXP ref, 
	SEXP tolerance, 
	SEXP relative, 
	SEXP nearest, 
	SEXP nomatch)
{
	if ( TYPEOF(query) != TYPEOF(ref) )
		Rf_error("'query' and 'ref' must have the same data type");
	SEXP index = PROTECT(Rf_allocVector(INTSXP, LENGTH(query)));
	switch(TYPEOF(ref))
	{
		case INTSXP:
			binary_search(
				vec<int>::from(index),
				vec<int>::from(query),
				vec<int>::from(ref),
				Rf_asReal(tolerance),
				Rf_asLogical(relative),
				Rf_asLogical(nearest),
				Rf_asInteger(nomatch));
			break;
		case REALSXP:
			binary_search(
				vec<int>::from(index),
				vec<double>::from(query),
				vec<double>::from(ref),
				Rf_asReal(tolerance),
				Rf_asLogical(relative),
				Rf_asLogical(nearest),
				Rf_asInteger(nomatch));
			break;
		default:
			Rf_error("'ref' must be integer or double");
	}
	UNPROTECT(1);
	return index;
}

SEXP do_kdtree_build(SEXP data)
{
	SEXP tree = PROTECT(Rf_allocVector(VECSXP, 4));
	SEXP left = PROTECT(Rf_allocVector(INTSXP, Rf_nrows(data)));
	SEXP right = PROTECT(Rf_allocVector(INTSXP, Rf_nrows(data)));
	ptrdiff_t root;
	SET_VECTOR_ELT(tree, 0, data);
	SET_VECTOR_ELT(tree, 1, left);
	SET_VECTOR_ELT(tree, 2, right);
	SET_VECTOR_ELT(tree, 3, Rf_ScalarInteger(NA_INTEGER));
	SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
	SET_STRING_ELT(names, 0, Rf_mkChar("data"));
	SET_STRING_ELT(names, 1, Rf_mkChar("left"));
	SET_STRING_ELT(names, 2, Rf_mkChar("right"));
	SET_STRING_ELT(names, 3, Rf_mkChar("root"));
	Rf_setAttrib(tree, R_NamesSymbol, names);
	Rf_setAttrib(tree, R_ClassSymbol, Rf_mkString("kdtree"));
	switch(TYPEOF(data))
	{
		case INTSXP:
			root = kdtree<int,int>::from(tree).build();
			break;
		case REALSXP:
			root = kdtree<int,double>::from(tree).build();
			break;
		default:
			Rf_error("'data' must be integer or double");
	}
	SET_VECTOR_ELT(tree, 3, Rf_ScalarInteger(root));
	UNPROTECT(4);
	return tree;
}

SEXP do_kdtree_range_search(
	SEXP query,
	SEXP tree,
	SEXP tolerance,
	SEXP relative,
	SEXP num_threads)
{
	SEXP data = VECTOR_ELT(tree, 0);
	if ( TYPEOF(query) != TYPEOF(data) )
		Rf_error("'query' and 'data' must have the same data type");
	if ( Rf_ncols(query) != Rf_ncols(data) )
		Rf_error("'query' and 'data' must have the same number of cols");
	if ( LENGTH(tolerance) != Rf_ncols(data) )
		Rf_error("length of 'tolerance' must match ncol(data)");
	if ( LENGTH(relative) != Rf_ncols(data) )
		Rf_error("length of 'relative' must match ncol(data)");
	SEXP counts = PROTECT(Rf_allocVector(INTSXP, Rf_nrows(query)));
	switch(TYPEOF(query))
	{
		case INTSXP:
			compute(
				range_counts{
					kdtree<int,int>::from(tree),
					vec<int>::from(counts),
					mat<int>::from(query),
					vec<double>::from(tolerance),
					vec<int>::from(relative)},
				Rf_asInteger(num_threads));
			break;
		case REALSXP:
			compute(
				range_counts{
					kdtree<int,double>::from(tree),
					vec<int>::from(counts),
					mat<double>::from(query),
					vec<double>::from(tolerance),
					vec<int>::from(relative)},
				Rf_asInteger(num_threads));
			break;
	}
	SEXP offset = PROTECT(Rf_allocVector(INTSXP, 1 + XLENGTH(counts)));
	int * poffset = INTEGER(offset);
	int * pcounts = INTEGER(counts);
	for ( ptrdiff_t i = 0; i < XLENGTH(offset); ++i )
	{
		if ( i == 0 )
			poffset[i] = 0;
		else
			poffset[i] = poffset[i - 1] + pcounts[i - 1];
	}
	SEXP index = PROTECT(Rf_allocVector(INTSXP, poffset[XLENGTH(counts)]));
	switch(TYPEOF(query))
	{
		case INTSXP:
			compute(
				range_searches{
					kdtree<int,int>::from(tree),
					rag<int,int>::from(index, offset),
					mat<int>::from(query),
					vec<double>::from(tolerance),
					vec<int>::from(relative)},
				Rf_asInteger(num_threads));
			break;
		case REALSXP:
			compute(
				range_searches{
					kdtree<int,double>::from(tree),
					rag<int,int>::from(index, offset),
					mat<double>::from(query),
					vec<double>::from(tolerance),
					vec<int>::from(relative)},
				Rf_asInteger(num_threads));
			break;
	}
	vec<int>::from(index) += rep{1, XLENGTH(index)};
	SEXP hits = PROTECT(Rf_allocVector(VECSXP, 3));
	SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
	SET_VECTOR_ELT(hits, 0, index);
	SET_VECTOR_ELT(hits, 1, offset);
	SET_VECTOR_ELT(hits, 2, counts);
	SET_STRING_ELT(names, 0, Rf_mkChar("index"));
	SET_STRING_ELT(names, 1, Rf_mkChar("offset"));
	SET_STRING_ELT(names, 2, Rf_mkChar("counts"));
	Rf_setAttrib(hits, R_NamesSymbol, names);
	UNPROTECT(5);
	return hits;
}

//// Matrix statistics
//---------------------

SEXP do_col_sums(SEXP x, SEXP num_threads)
{
	SEXP sums = PROTECT(Rf_allocVector(REALSXP, Rf_ncols(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
		{
			compute(
				col_sums<int>{
					vec<double>::from(sums).fill(0),
					mat<int>::from(x)},
				Rf_asInteger(num_threads));
			break;
		}
		case REALSXP:
		{
			compute(
				col_sums<double>{
					vec<double>::from(sums).fill(0),
					mat<double>::from(x)},
				Rf_asInteger(num_threads));
			break;
		}
	}
	UNPROTECT(1);
	return sums;
}

//// Test expressions
//--------------------

SEXP do_test_expression(SEXP x, SEXP index)
{
	SEXP result = PROTECT(Rf_allocVector(TYPEOF(x), XLENGTH(index)));
	switch(TYPEOF(x))
	{
		case INTSXP:
		{
			test_expression(
				vec<int>::from(result),
				vec<int>::from(x),
				vec<int>::from(index));
			break;
		}
		case REALSXP:
		{
			test_expression(
				vec<double>::from(result),
				vec<double>::from(x),
				vec<int>::from(index));
			break;
		}
	}
	UNPROTECT(1);
	return result;
}

} // extern "C"
