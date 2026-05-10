
#include "CardinalCore.h"

extern "C" {

//// Quicksort and Quickselect
//----------------------------

SEXP do_qdiff(SEXP x, SEXP ref, SEXP relative)
{
	if ( TYPEOF(x) != TYPEOF(ref) )
		Rf_error("'x' and 'ref' must have the same data type");
	if ( LENGTH(x) != LENGTH(ref) )
		Rf_error("'x' and 'ref' must have the same length");
	SEXP result = PROTECT(Rf_allocVector(REALSXP, LENGTH(x)));
	for ( R_len_t i = 0; i < LENGTH(x); ++i )
	{
		switch(TYPEOF(x))
		{
			case INTSXP:
				REAL(result)[i] = diff<int>(
					INTEGER_ELT(x, i),
					INTEGER_ELT(ref, i),
					Rf_asLogical(relative));
				break;
			case REALSXP:
				REAL(result)[i] = diff<double>(
					REAL_ELT(x, i),
					REAL_ELT(ref, i),
					Rf_asLogical(relative));
				break;
			default:
				Rf_error("'x' and 'ref' must be integer or double");
		}
	}
	UNPROTECT(1);
	return result;
}

SEXP do_qselect(SEXP x, SEXP k)
{
	SEXP result = PROTECT(Rf_allocVector(TYPEOF(x), LENGTH(k)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			quick_select<int,int,int>(
				as_vctr<int>(x),
				as_vctr<int>(k),
				INTEGER(result));
			break;
		case REALSXP:
			quick_select<double,int,int>(
				as_vctr<double>(x),
				as_vctr<int>(k),
				REAL(result));
			break;
		default:
			Rf_error("'x' must be integer or double");
	}
	UNPROTECT(1);
	return result;
}

SEXP do_qorder(SEXP x)
{
	SEXP result = PROTECT(Rf_allocVector(INTSXP, LENGTH(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			quick_order<int,int>(
				as_vctr<int>(x),
				INTEGER(result));
			break;
		case REALSXP:
			quick_order<double,int>(
				as_vctr<double>(x),
				INTEGER(result));
			break;
		default:
			Rf_error("'x' must be integer or double");
	}
	UNPROTECT(1);
	return result;
}

SEXP do_qmedian(SEXP x)
{
	switch(TYPEOF(x))
	{
		case INTSXP:
			return Rf_ScalarReal(quick_median<int,int>(
				as_vctr<int>(x)));
			break;
		case REALSXP:
			return Rf_ScalarReal(quick_median<double,int>(
				as_vctr<double>(x)));
		default:
			Rf_error("'x' must be integer or double");
	}
}

SEXP do_qmad(SEXP x, SEXP center, SEXP scale)
{
	switch(TYPEOF(x))
	{
		case INTSXP:
			return Rf_ScalarReal(quick_mad<int,int>(
				as_vctr<int>(x),
				Rf_asReal(center),
				Rf_asReal(scale)));
		case REALSXP:
			return Rf_ScalarReal(quick_mad<double,int>(
				as_vctr<double>(x),
				Rf_asReal(center),
				Rf_asReal(scale)));
		default:
			Rf_error("'x' must be integer or double");
	}
}

SEXP do_bsearch(
	SEXP x, 
	SEXP data, 
	SEXP tolerance, 
	SEXP relative, 
	SEXP nearest, 
	SEXP nomatch)
{
	if ( TYPEOF(x) != TYPEOF(data) )
		Rf_error("'x' and 'data' must have the same data type");
	SEXP result = PROTECT(Rf_allocVector(INTSXP, LENGTH(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
			binary_search<int,int>(
				as_vctr<int>(x),
				as_vctr<int>(data),
				INTEGER(result),
				Rf_asReal(tolerance),
				Rf_asLogical(relative),
				Rf_asLogical(nearest),
				Rf_asInteger(nomatch));
			break;
		case REALSXP:
			binary_search<double,int>(
				as_vctr<double>(x),
				as_vctr<double>(data),
				INTEGER(result),
				Rf_asReal(tolerance),
				Rf_asLogical(relative),
				Rf_asLogical(nearest),
				Rf_asInteger(nomatch));
			break;
		default:
			Rf_error("'x' must be integer or double");
	}
	UNPROTECT(1);
	return result;
}

SEXP do_col_sums(SEXP x, SEXP num_threads)
{
	SEXP result = PROTECT(Rf_allocVector(REALSXP, Rf_ncols(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
		{
			col_sums<int>(
				as_matrix<int>(x), 
				REAL(result), 
				Rf_asInteger(num_threads));
		}
		case REALSXP:
		{
			col_sums<double>(
				as_matrix<double>(x), 
				REAL(result), 
				Rf_asInteger(num_threads));
		}
	}
	UNPROTECT(1);
	return result;
}

SEXP do_col_scatter_sums(
	SEXP x, 
	SEXP group,
	SEXP ngroups,
	SEXP num_threads)
{
	SEXP result = PROTECT(Rf_allocMatrix(REALSXP, 
		Rf_asInteger(ngroups), Rf_ncols(x)));
	switch(TYPEOF(x))
	{
		case INTSXP:
		{
			col_scatter_sums<int>(
				as_matrix<int>(x), 
				INTEGER(group),
				Rf_asInteger(ngroups),
				REAL(result), 
				Rf_asInteger(num_threads));
		}
		case REALSXP:
		{
			col_scatter_sums<double>(
				as_matrix<double>(x), 
				INTEGER(group),
				Rf_asInteger(ngroups),
				REAL(result), 
				Rf_asInteger(num_threads));
		}
	}
	UNPROTECT(1);
	return result;
}

} // extern "C"
