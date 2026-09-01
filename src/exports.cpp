
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include "core.h"
#include "kernels.h"
#include "order.h"
#include "search.h"
#include "signal.h"
#include "peaks.h"

//// Exports
//-----------

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
	add1(r_vec<int>(index));
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
	add1(r_mat<int>(ends).col(0));
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
		default:
			Rf_error("'query' and 'table' must be integer or double");
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
		default:
			Rf_error("'query' and 'table' must be integer or double");
	}
	add1(r_vec<int>(index));
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
		default:
			Rf_error("'query' and 'table' must be integer or double");
	}
	add1(r_vec<int>(index));
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

//// Signal processing
//--------------------

SEXP do_filt1_mean(SEXP y, SEXP k)
{
	SEXP yout = PROTECT(Rf_allocVector(REALSXP, XLENGTH(y)));
	switch(TYPEOF(y))
	{
		case INTSXP:
			filt1_mean(
				r_vec<double>(yout),
				r_vec<int>(y),
				Rf_asInteger(k));
			break;
		case REALSXP:
			filt1_mean(
				r_vec<double>(yout),
				r_vec<double>(y),
				Rf_asInteger(k));
			break;
		default:
			Rf_error("'y' must be integer or double");
	}
	UNPROTECT(1);
	return yout;
}

SEXP do_filt1_conv(SEXP y, SEXP w)
{
	SEXP yout = PROTECT(Rf_allocVector(REALSXP, XLENGTH(y)));
	switch(TYPEOF(y))
	{
		case INTSXP:
			filt1_conv(
				r_vec<double>(yout),
				r_vec<int>(y),
				r_vec<double>(w));
			break;
		case REALSXP:
			filt1_conv(
				r_vec<double>(yout),
				r_vec<double>(y),
				r_vec<double>(w));
			break;
		default:
			Rf_error("'y' must be integer or double");
	}
	UNPROTECT(1);
	return yout;
}

//// Peak processing
//------------------

SEXP do_peaks_find(SEXP y, SEXP k)
{
	int count = 0;
	switch(TYPEOF(y))
	{
		case INTSXP:
			count = peaks{r_vec<int>(y), Rf_asInteger(k)}.count();
			break;
		case REALSXP:
			count = peaks{r_vec<double>(y), Rf_asInteger(k)}.count();
			break;
		default:
			Rf_error("'y' must be integer or double");
	}
	SEXP index = PROTECT(Rf_allocVector(INTSXP, count));
	switch(TYPEOF(y))
	{
		case INTSXP:
			peaks{r_vec<int>(y), Rf_asInteger(k)}
				.index_into(r_vec<int>(index));
			break;
		case REALSXP:
			peaks{r_vec<double>(y), Rf_asInteger(k)}
				.index_into(r_vec<int>(index));
			break;
	}
	add1(r_vec<int>(index));
	UNPROTECT(1);
	return index;
}

SEXP do_peaks_sums(SEXP y, SEXP k)
{
	int count = 0;
	switch(TYPEOF(y))
	{
		case INTSXP:
			count = peaks{r_vec<int>(y), Rf_asInteger(k)}.count();
			break;
		case REALSXP:
			count = peaks{r_vec<double>(y), Rf_asInteger(k)}.count();
			break;
		default:
			Rf_error("'y' must be integer or double");
	}
	SEXP index = PROTECT(Rf_allocVector(INTSXP, count));
	SEXP left_end = PROTECT(Rf_allocVector(INTSXP, count));
	SEXP right_end = PROTECT(Rf_allocVector(INTSXP, count));
	SEXP sums = PROTECT(Rf_allocVector(REALSXP, count));
	switch(TYPEOF(y))
	{
		case INTSXP:
			peaks{r_vec<int>(y), Rf_asInteger(k)}
				.sums_into(
					r_vec<int>(index),
					r_vec<int>(left_end),
					r_vec<int>(right_end),
					r_vec<double>(sums));
			break;
		case REALSXP:
			peaks{r_vec<double>(y), Rf_asInteger(k)}
				.sums_into(
					r_vec<int>(index),
					r_vec<int>(left_end),
					r_vec<int>(right_end),
					r_vec<double>(sums));
			break;
	}
	add1(r_vec<int>(index));
	add1(r_vec<int>(left_end));
	add1(r_vec<int>(right_end));
	SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
	SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
	SET_VECTOR_ELT(out, 0, index);
	SET_VECTOR_ELT(out, 1, left_end);
	SET_VECTOR_ELT(out, 2, right_end);
	SET_VECTOR_ELT(out, 3, sums);
	SET_STRING_ELT(names, 0, Rf_mkChar("index"));
	SET_STRING_ELT(names, 1, Rf_mkChar("left_end"));
	SET_STRING_ELT(names, 2, Rf_mkChar("right_end"));
	SET_STRING_ELT(names, 3, Rf_mkChar("sum"));
	Rf_setAttrib(out, R_NamesSymbol, names);
	UNPROTECT(6);
	return out;
}

SEXP do_peaks_prominences(SEXP y, SEXP k, SEXP wlen)
{
	int count = 0;
	switch(TYPEOF(y))
	{
		case INTSXP:
			count = peaks{r_vec<int>(y), Rf_asInteger(k)}.count();
			break;
		case REALSXP:
			count = peaks{r_vec<double>(y), Rf_asInteger(k)}.count();
			break;
		default:
			Rf_error("'y' must be integer or double");
	}
	SEXP index = PROTECT(Rf_allocVector(INTSXP, count));
	SEXP left_base = PROTECT(Rf_allocVector(INTSXP, count));
	SEXP right_base = PROTECT(Rf_allocVector(INTSXP, count));
	SEXP prominences = PROTECT(Rf_allocVector(REALSXP, count));
	switch(TYPEOF(y))
	{
		case INTSXP:
			peaks{r_vec<int>(y), Rf_asInteger(k)}
				.prominences_into(
					r_vec<int>(index),
					r_vec<int>(left_base),
					r_vec<int>(right_base),
					r_vec<double>(prominences),
					Rf_asInteger(wlen));
			break;
		case REALSXP:
			peaks{r_vec<double>(y), Rf_asInteger(k)}
				.prominences_into(
					r_vec<int>(index),
					r_vec<int>(left_base),
					r_vec<int>(right_base),
					r_vec<double>(prominences),
					Rf_asInteger(wlen));
			break;
	}
	add1(r_vec<int>(index));
	add1(r_vec<int>(left_base));
	add1(r_vec<int>(right_base));
	SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
	SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
	SET_VECTOR_ELT(out, 0, index);
	SET_VECTOR_ELT(out, 1, left_base);
	SET_VECTOR_ELT(out, 2, right_base);
	SET_VECTOR_ELT(out, 3, prominences);
	SET_STRING_ELT(names, 0, Rf_mkChar("index"));
	SET_STRING_ELT(names, 1, Rf_mkChar("left_base"));
	SET_STRING_ELT(names, 2, Rf_mkChar("right_base"));
	SET_STRING_ELT(names, 3, Rf_mkChar("prominence"));
	Rf_setAttrib(out, R_NamesSymbol, names);
	UNPROTECT(6);
	return out;
}

SEXP do_peaks_areas(SEXP y, SEXP x, SEXP k)
{
	int count = 0;
	switch(TYPEOF(y))
	{
		case INTSXP:
			count = peaks{r_vec<int>(y), Rf_asInteger(k)}.count();
			break;
		case REALSXP:
			count = peaks{r_vec<double>(y), Rf_asInteger(k)}.count();
			break;
		default:
			Rf_error("'y' must be integer or double");
	}
	SEXP index = PROTECT(Rf_allocVector(INTSXP, count));
	SEXP left_end = PROTECT(Rf_allocVector(INTSXP, count));
	SEXP right_end = PROTECT(Rf_allocVector(INTSXP, count));
	SEXP areas = PROTECT(Rf_allocVector(REALSXP, count));
	switch(TYPEOF(y))
	{
		case INTSXP:
			peaks{r_vec<int>(y), Rf_asInteger(k)}
				.areas_into(
					r_vec<int>(index),
					r_vec<int>(left_end),
					r_vec<int>(right_end),
					r_vec<double>(areas),
					r_vec<int>(x));
			break;
		case REALSXP:
			peaks{r_vec<double>(y), Rf_asInteger(k)}
				.areas_into(
					r_vec<int>(index),
					r_vec<int>(left_end),
					r_vec<int>(right_end),
					r_vec<double>(areas),
					r_vec<double>(x));
			break;
	}
	add1(r_vec<int>(index));
	add1(r_vec<int>(left_end));
	add1(r_vec<int>(right_end));
	SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
	SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
	SET_VECTOR_ELT(out, 0, index);
	SET_VECTOR_ELT(out, 1, left_end);
	SET_VECTOR_ELT(out, 2, right_end);
	SET_VECTOR_ELT(out, 3, areas);
	SET_STRING_ELT(names, 0, Rf_mkChar("index"));
	SET_STRING_ELT(names, 1, Rf_mkChar("left_end"));
	SET_STRING_ELT(names, 2, Rf_mkChar("right_end"));
	SET_STRING_ELT(names, 3, Rf_mkChar("area"));
	Rf_setAttrib(out, R_NamesSymbol, names);
	UNPROTECT(6);
	return out;
}

SEXP do_peaks_widths(SEXP y, SEXP x, SEXP k, SEXP fmax)
{
	int count = 0;
	switch(TYPEOF(y))
	{
		case INTSXP:
			count = peaks{r_vec<int>(y), Rf_asInteger(k)}.count();
			break;
		case REALSXP:
			count = peaks{r_vec<double>(y), Rf_asInteger(k)}.count();
			break;
		default:
			Rf_error("'y' must be integer or double");
	}
	SEXP index = PROTECT(Rf_allocVector(INTSXP, count));
	SEXP left_ips = PROTECT(Rf_allocVector(REALSXP, count));
	SEXP right_ips = PROTECT(Rf_allocVector(REALSXP, count));
	SEXP widths = PROTECT(Rf_allocVector(REALSXP, count));
	switch(TYPEOF(y))
	{
		case INTSXP:
			peaks{r_vec<int>(y), Rf_asInteger(k)}
				.widths_into(
					r_vec<int>(index),
					r_vec<double>(left_ips),
					r_vec<double>(right_ips),
					r_vec<double>(widths),
					r_vec<int>(x),
					Rf_asReal(fmax));
			break;
		case REALSXP:
			peaks{r_vec<double>(y), Rf_asInteger(k)}
				.widths_into(
					r_vec<int>(index),
					r_vec<double>(left_ips),
					r_vec<double>(right_ips),
					r_vec<double>(widths),
					r_vec<double>(x),
					Rf_asReal(fmax));
			break;
	}
	add1(r_vec<int>(index));
	add1(r_vec<double>(left_ips));
	add1(r_vec<double>(right_ips));
	SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
	SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
	SET_VECTOR_ELT(out, 0, index);
	SET_VECTOR_ELT(out, 1, left_ips);
	SET_VECTOR_ELT(out, 2, right_ips);
	SET_VECTOR_ELT(out, 3, widths);
	SET_STRING_ELT(names, 0, Rf_mkChar("index"));
	SET_STRING_ELT(names, 1, Rf_mkChar("left_ips"));
	SET_STRING_ELT(names, 2, Rf_mkChar("right_ips"));
	SET_STRING_ELT(names, 3, Rf_mkChar("width"));
	Rf_setAttrib(out, R_NamesSymbol, names);
	UNPROTECT(6);
	return out;
}

SEXP do_peaks_centroids(SEXP y, SEXP x, SEXP k)
{
	int count = 0;
	switch(TYPEOF(y))
	{
		case INTSXP:
			count = peaks{r_vec<int>(y), Rf_asInteger(k)}.count();
			break;
		case REALSXP:
			count = peaks{r_vec<double>(y), Rf_asInteger(k)}.count();
			break;
		default:
			Rf_error("'y' must be integer or double");
	}
	SEXP index = PROTECT(Rf_allocVector(INTSXP, count));
	SEXP left_end = PROTECT(Rf_allocVector(INTSXP, count));
	SEXP right_end = PROTECT(Rf_allocVector(INTSXP, count));
	SEXP centroids = PROTECT(Rf_allocVector(REALSXP, count));
	switch(TYPEOF(y))
	{
		case INTSXP:
			peaks{r_vec<int>(y), Rf_asInteger(k)}
				.centroids_into(
					r_vec<int>(index),
					r_vec<int>(left_end),
					r_vec<int>(right_end),
					r_vec<double>(centroids),
					r_vec<int>(x));
			break;
		case REALSXP:
			peaks{r_vec<double>(y), Rf_asInteger(k)}
				.centroids_into(
					r_vec<int>(index),
					r_vec<int>(left_end),
					r_vec<int>(right_end),
					r_vec<double>(centroids),
					r_vec<double>(x));
			break;
	}
	add1(r_vec<int>(index));
	add1(r_vec<int>(left_end));
	add1(r_vec<int>(right_end));
	SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
	SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
	SET_VECTOR_ELT(out, 0, index);
	SET_VECTOR_ELT(out, 1, left_end);
	SET_VECTOR_ELT(out, 2, right_end);
	SET_VECTOR_ELT(out, 3, centroids);
	SET_STRING_ELT(names, 0, Rf_mkChar("index"));
	SET_STRING_ELT(names, 1, Rf_mkChar("left_end"));
	SET_STRING_ELT(names, 2, Rf_mkChar("right_end"));
	SET_STRING_ELT(names, 3, Rf_mkChar("centroid"));
	Rf_setAttrib(out, R_NamesSymbol, names);
	UNPROTECT(6);
	return out;
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
		default:
			Rf_error("'x' must be integer or double");
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
		default:
			Rf_error("'x' must be integer or double");
	}
	UNPROTECT(1);
	return result;
}

} // extern "C"
