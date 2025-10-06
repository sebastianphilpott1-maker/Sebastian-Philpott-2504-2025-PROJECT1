#############################################################################
#############################################################################
#
# This file implements polynomial division for dense polynomials.
#                                                                               
#############################################################################
#############################################################################

# We won't re-implement any of these functions for dense polynomials, the abstract versions will 
# produce the correct result.

# div_rem_mod_p(num::PolynomialDense{C, D}, den::PolynomialDense{C, D}, prime::Integer)
# div_mod_p(num::Polynomial, den::Polynomial, prime::Integer)
# rem_mod_p(num::Polynomial, den::Polynomial, prime::Integer)