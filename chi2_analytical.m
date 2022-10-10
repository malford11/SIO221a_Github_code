function chi2_analytical(vals,n)
%function chi2_analytical(vals,n)
%Return the PDF of a chi2 variable with n degrees of freedom, for the
%values specified in vals.
chi2_analytical = 1./2.^(n/2)./gamma(n/2).*exp(-vals/2).*vals.^(n/2-1);

