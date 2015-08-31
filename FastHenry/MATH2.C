/* Enrico, asinh() and atanh() functions, not provided */
/* by the standard Microsoft libraries */

#include <math.h>
#include <stdio.h>

#define LN2 6.93147180559945286227e-01
#define _2POW28 268435456.0

// not needed any more, now math.h defines these functions

/*
double asinh(double arg)
{
	double farg, fres;

	// the base formula is asinh(x) = ln(x+sqrt(x^2+1))
	//
	// however, since asinh(-x) = -asinh(x), to be sure
	// to have symmetric results, it's better to use
	// asinh(x) = sign(x)*ln(|x|+sqrt(x^2+1))
	//
	// moreover, for large x the results tend to be not
	// accurate; therefore, for large x (that is, |x| > 2^28),
	// it's better to use the approximate formula
	// asinh(x) = sign(x)*(ln(|x|) + ln(2))
	//
	// for intermediate values, it's better to use the formula
	// asinh(x) = sign(x)*ln(2|x|+1/(|x|+sqrt(x^2+1))),
	// since if the second term tend to be zero, the first
	// term still stays; this formula is obtained from the
	// base one by adding and subtracting terms
	// (i.e. add and subtract x+sqrt(x^2+1) to this 
	// formula and you'll obtain the base one)
	//
	// on the other hand, for small x (|x|<2) it would be better
	// to use a specialized formula for ln(1+y) (where
	// y is sqrt(x^2+1)) since the standard log is not
	// very accurate in this case. The specialized formula
	// is called log1p on Unix systems but on Windows there's
	// no equivalent and is not easy to implement, so
	// it is not used here

	farg = fabs(arg);
	// if |arg| > 2^28, use approximation formula
	if(farg > _2POW28) {
		fres = log(farg) + LN2;
	}
	else {
		fres = log(2*farg + 1/(farg + sqrt(farg * farg + 1)));
	}

	if(arg > 0)
		return(fres);
	else
		return(-fres);
}

double atanh(double arg)
{
	double farg, fres;

	// the base formula is atanh(x) = 0.5*ln((1+x)/(1-x))
	//
	// however, since atanh(-x) = -atanh(x), to be sure
	// to have symmetric results, it's better to use
	// atanh(x) = sign(x)*0.5*ln((1+|x|)/(1-|x|))
	//
	// note that for small x (as is the case here, since |x| must
	// be < 1), it would be better to use a specialized formula 
	// for ln(1+y), elaborating the base formula, adding and 
	// subtracting terms, to get 0.5*ln(1+2*x/(1-x)), so that
	// y is 2*x/(1-x), since the standard log is not
	// very accurate in this case. The specialized formula
	// is called log1p on Unix systems but on Windows there's
	// no equivalent and is not easy to implement, so
	// it is not used here

	farg = fabs(arg);
	// check argument
	if( farg >= 1.0 ) {
		viewprintf(stderr, "Warning: bad atanh argument %e!\n", arg);
		return(0);
	}

	fres = (1 + farg) / (1 - farg);

	// should never happen!
	if( fres < 0 ) {
		viewprintf(stderr, "Warning: bad atanh argument %e!\n", arg);
		return(0);
	}

	if(arg > 0)
		return(log(fres)/2);
	else
		return(-log(fres)/2);
}

*/