
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include "order.h"
#include "search.h"
#include "kernels.h"

extern "C" {

//// Quicksort and Quickselect
//----------------------------

SEXP do_qorder(SEXP x)
{
	SEXP index = PROTECT(Rf_allocVector(INTSXP, LENGTH(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			qsort_index(
				r_vec<int>(index).fill_seq(),
				r_vec<int>(x));
			break;
		case REALSXP:
			qsort_index(
				r_vec<int>(index).fill_seq(),
				r_vec<double>(x));
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
					r_vec<int>(index).fill_seq(),
					r_vec<int>(x),
					INTEGER_ELT(k, i));
				break;
			case REALSXP:
				REAL(order)[i] = qselect_index(
					r_vec<int>(index).fill_seq(),
					r_vec<double>(x),
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
			median = qmedian(r_vec<int>(x));
			break;
		case REALSXP:
			median = qmedian(r_vec<double>(x));
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
				r_vec<int>(x),
				Rf_asReal(center),
				Rf_asReal(constant));
			break;
		case REALSXP:
			mad = qmad(
				r_vec<double>(x),
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
	SEXP table,
	SEXP tolerance,
	SEXP relative,
	SEXP referent,
	SEXP nomatch)
{
	if ( TYPEOF(query) != TYPEOF(table) )
		Rf_error("'query' and 'table' must have the same data type");
	SEXP index = PROTECT(Rf_allocVector(INTSXP, LENGTH(query)));
	switch(TYPEOF(table))
	{
		case INTSXP:
			bsearch(
				r_vec<int>(index),
				r_vec<int>(query),
				r_vec<int>(table),
				Rf_asReal(tolerance),
				Rf_asLogical(relative),
				static_cast<Ref>(Rf_asInteger(referent)),
				Rf_asInteger(nomatch));
			break;
		case REALSXP:
			bsearch(
				r_vec<int>(index),
				r_vec<double>(query),
				r_vec<double>(table),
				Rf_asReal(tolerance),
				Rf_asLogical(relative),
				static_cast<Ref>(Rf_asInteger(referent)),
				Rf_asInteger(nomatch));
			break;
		default:
			Rf_error("'query' and 'table' must be integer or double");
	}
	r_vec<int>(index) += 1;
	UNPROTECT(1);
	return index;
}

SEXP do_rsearch(
	SEXP query,
	SEXP table,
	SEXP tolerance,
	SEXP relative,
	SEXP referent,
	SEXP nomatch)
{
	if ( TYPEOF(query) != TYPEOF(table) )
		Rf_error("'query' and 'table' must have the same data type");
	SEXP ends = PROTECT(Rf_allocMatrix(INTSXP, LENGTH(query), 2));
	switch(TYPEOF(table))
	{
		case INTSXP:
			rsearch(
				r_mat<int>(ends).col(0),
				r_mat<int>(ends).col(1),
				r_vec<int>(query),
				r_vec<int>(table),
				Rf_asReal(tolerance),
				Rf_asLogical(relative),
				static_cast<Ref>(Rf_asInteger(referent)),
				Rf_asInteger(nomatch));
			break;
		case REALSXP:
			rsearch(
				r_mat<int>(ends).col(0),
				r_mat<int>(ends).col(1),
				r_vec<double>(query),
				r_vec<double>(table),
				Rf_asReal(tolerance),
				Rf_asLogical(relative),
				static_cast<Ref>(Rf_asInteger(referent)),
				Rf_asInteger(nomatch));
			break;
		default:
			Rf_error("'query' and 'table' must be integer or double");
	}
	r_mat<int>(ends).col(0) += 1;
	SEXP colnames = PROTECT(Rf_allocVector(STRSXP, 2));
	SET_STRING_ELT(colnames, 0, Rf_mkChar("start"));
	SET_STRING_ELT(colnames, 1, Rf_mkChar("end"));
	SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dimnames, 0, R_NilValue);
	SET_VECTOR_ELT(dimnames, 1, colnames);
	Rf_setAttrib(ends, R_DimNamesSymbol, dimnames);
	UNPROTECT(3);
	return ends;
}

SEXP do_kdtree_build(SEXP table)
{
	SEXP tree = PROTECT(Rf_allocVector(VECSXP, 4));
	SEXP left = PROTECT(Rf_allocVector(INTSXP, Rf_nrows(table)));
	SEXP right = PROTECT(Rf_allocVector(INTSXP, Rf_nrows(table)));
	ptrdiff_t root;
	SET_VECTOR_ELT(tree, 0, table);
	SET_VECTOR_ELT(tree, 1, left);
	SET_VECTOR_ELT(tree, 2, right);
	SET_VECTOR_ELT(tree, 3, Rf_ScalarInteger(NA_INTEGER));
	SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
	SET_STRING_ELT(names, 0, Rf_mkChar("table"));
	SET_STRING_ELT(names, 1, Rf_mkChar("left"));
	SET_STRING_ELT(names, 2, Rf_mkChar("right"));
	SET_STRING_ELT(names, 3, Rf_mkChar("root"));
	Rf_setAttrib(tree, R_NamesSymbol, names);
	Rf_setAttrib(tree, R_ClassSymbol, Rf_mkString("kdtree"));
	switch(TYPEOF(table))
	{
		case INTSXP:
			root = kdtree<int,int>::from(tree).build();
			break;
		case REALSXP:
			root = kdtree<int,double>::from(tree).build();
			break;
		default:
			Rf_error("'table' must be integer or double");
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
	SEXP referent,
	SEXP num_threads)
{
	SEXP table = VECTOR_ELT(tree, 0);
	if ( TYPEOF(query) != TYPEOF(table) )
		Rf_error("'query' and 'table' must have the same data type");
	if ( Rf_ncols(query) != Rf_ncols(table) )
		Rf_error("'query' and 'table' must have the same number of cols");
	if ( LENGTH(tolerance) != Rf_ncols(table) )
		Rf_error("length of 'tolerance' must match ncol(table)");
	if ( LENGTH(relative) != Rf_ncols(table) )
		Rf_error("length of 'relative' must match ncol(table)");
	SEXP counts = PROTECT(Rf_allocVector(INTSXP, Rf_nrows(query)));
	switch(TYPEOF(query))
	{
		case INTSXP:
			compute(
				range_counts{
					r_vec<int>(counts),
					kdtree<int,int>::from(tree),
					r_mat<int>(query),
					r_vec<double>(tolerance),
					r_vec<int>(relative),
					static_cast<Ref>(Rf_asInteger(referent)),
				},
				Rf_asInteger(num_threads));
			break;
		case REALSXP:
			compute(
				range_counts{
					r_vec<int>(counts),
					kdtree<int,double>::from(tree),
					r_mat<double>(query),
					r_vec<double>(tolerance),
					r_vec<int>(relative),
					static_cast<Ref>(Rf_asInteger(referent)),
				},
				Rf_asInteger(num_threads));
			break;
	}
	SEXP offset = PROTECT(Rf_allocVector(INTSXP, 1 + XLENGTH(counts)));
	int * poffset = INTEGER(offset);
	int * pcounts = INTEGER(counts);
	for ( ptrdiff_t i = 0; i < XLENGTH(offset); ++i ) {
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
					r_vecs_pack<int,int>(index, offset),
					kdtree<int,int>::from(tree),
					r_mat<int>(query),
					r_vec<double>(tolerance),
					r_vec<int>(relative),
					static_cast<Ref>(Rf_asInteger(referent)),
				},
				Rf_asInteger(num_threads));
			break;
		case REALSXP:
			compute(
				range_searches{
					r_vecs_pack<int,int>(index, offset),
					kdtree<int,double>::from(tree),
					r_mat<double>(query),
					r_vec<double>(tolerance),
					r_vec<int>(relative),
					static_cast<Ref>(Rf_asInteger(referent)),
				},
				Rf_asInteger(num_threads));
			break;
	}
	r_vec<int>(index) += 1;
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

SEXP do_kdtree_knn_search(
	SEXP query,
	SEXP tree,
	SEXP k,
	SEXP p,
	SEXP num_threads)
{
	SEXP table = VECTOR_ELT(tree, 0);
	if ( TYPEOF(query) != TYPEOF(table) )
		Rf_error("'query' and 'table' must have the same data type");
	if ( Rf_ncols(query) != Rf_ncols(table) )
		Rf_error("'query' and 'table' must have the same number of cols");
	SEXP offset = PROTECT(Rf_allocVector(INTSXP, 1 + XLENGTH(k)));
	int * poffset = INTEGER(offset);
	int * pcounts = INTEGER(k);
	for ( ptrdiff_t i = 0; i < XLENGTH(offset); ++i ) {
		if ( i == 0 )
			poffset[i] = 0;
		else
			poffset[i] = poffset[i - 1] + pcounts[i - 1];
	}
	SEXP index = PROTECT(Rf_allocVector(INTSXP, poffset[XLENGTH(k)]));
	SEXP dists = PROTECT(Rf_allocVector(REALSXP, poffset[XLENGTH(k)]));
	switch(TYPEOF(query))
	{
		case INTSXP:
			compute(
				knn_searches{
					r_vecs_pack<int,int>(index, offset),
					r_vecs_pack<double,int>(dists, offset),
					kdtree<int,int>::from(tree),
					r_mat<int>(query),
					static_cast<Norm>(Rf_asInteger(p)),
				},
				Rf_asInteger(num_threads));
			break;
		case REALSXP:
			compute(
				knn_searches{
					r_vecs_pack<int,int>(index, offset),
					r_vecs_pack<double,int>(dists, offset),
					kdtree<int,double>::from(tree),
					r_mat<double>(query),
					static_cast<Norm>(Rf_asInteger(p)),
				},
				Rf_asInteger(num_threads));
			break;
	}
	r_vec<int>(index) += 1;
	SEXP hits = PROTECT(Rf_allocVector(VECSXP, 4));
	SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
	SET_VECTOR_ELT(hits, 0, index);
	SET_VECTOR_ELT(hits, 1, dists);
	SET_VECTOR_ELT(hits, 2, offset);
	SET_VECTOR_ELT(hits, 3, k);
	SET_STRING_ELT(names, 0, Rf_mkChar("index"));
	SET_STRING_ELT(names, 1, Rf_mkChar("distance"));
	SET_STRING_ELT(names, 2, Rf_mkChar("offset"));
	SET_STRING_ELT(names, 3, Rf_mkChar("counts"));
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
				col_sums{
					r_vec<double>(sums).fill(0),
					r_mat<int>(x)},
				Rf_asInteger(num_threads));
			break;
		}
		case REALSXP:
		{
			compute(
				col_sums{
					r_vec<double>(sums).fill(0),
					r_mat<double>(x)},
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
				r_vec<int>(result),
				r_vec<int>(x),
				r_vec<int>(index));
			break;
		}
		case REALSXP:
		{
			test_expression(
				r_vec<double>(result),
				r_vec<double>(x),
				r_vec<int>(index));
			break;
		}
	}
	UNPROTECT(1);
	return result;
}

} // extern "C"
