# q-special-functions

programs for q-special functions and q-series

Thank you for checking my programs.

Programs in this repository are made to calculate q-special functions and q-series with guaranteed accuracy. All q-special functions are computed by calculating infinite q-Pochhammer symbols and q-hypergeometric series with guaranteed accuaracy. 

Note: The function "complex_nbd" is defined in Heine.hpp. This works like "cintval" in INTLAB.

If you use my programs in publications, please include the following reference.

Daisuke Kanaizumi, Programs for q-special functions and q-series,
https://github.com/Daisuke-Kanaizumi/q-special-functions

You do NOT have to inform your usage to me, but I am very glad to hear how my programs helped you.

(Requirements)

Before using my C++ programs, please install these libraries.

kv library -a C++ library for verified numerical computation 

http://verifiedby.me/kv/index-e.html

https://github.com/mskashi/kv

Boost C++ Libraries http://www.boost.org/

If you want to use my MATLAB programs, please install INTLAB.

INTLAB - INTerval LABoratory
The Matlab/Octave toolbox for Reliable Computing
http://www.ti3.tu-harburg.de/rump/intlab/

DISCLAIMER: Extensive tests have been performed to ensure reliability of the algorithms (under Windows 10 with Visual Studio 2015, 2017 and Ubuntu 14.04). However, neither an error-free processor nor an error-free program can be guaranteed.

(q-special functions in this repository)

elliptic gamma function

little a-Jacobi function 

Heine hypergeometric function

q-Airy functions (HKW, KMNOY, IKMMS, Ramanujan)

q-Appell function (\Phi_1)

q-Barnes G function

q-Bessel functions (Jackson, Hahn-Exton, little, modified)
Jackson`s 2nd q-Bessel function can be computed by using the integration formulas derived by Rahman (1987) and Ismail-Zhang (2016,arXiv)

q-beta function

q-Bi function 

q-digamma

q-logarithm

q-gamma functions (Jackson, Moak)

q-Lauricella function (type D)

quantum dilogarithm Li_2(x;q) (defined by Kirillov)

quantum polylogarithm

Ramanujan psi sum

Ramanujan theta function

symmetric q-gamma, symmetric q-beta

(q-series in this repository)
Bailey Mod9, Dyson Mod27, Gessel-Stanton, Gollnitz-Gordon, Jackson-Slater ,Lebesgue, Rogers Mod14, Rogers-Ramanujan, Rogers-Selberg

(References)

S.M. Rump: INTLAB - INTerval LABoratory. 
In Tibor Csendes, editor, Developments in Reliable Computing, pages 77-104. 
Kluwer Academic Publishers, Dordrecht, 1999.

Mourad E. H. Ismail,
CLASSICAL AND QUANTUM ORTHOGONAL POLYNOMIALS IN ONE VARIABLE,
Encyclopedia of Mathematics and its Applications 98,
CAMBRIDGE UNIVERSITY PRESS, 2005.

George E. Andrews, Richard Askey, Ranjan Roy,
SPECIAL FUNCTIONS,
Encyclopedia of Mathematics and its Applications 71,
CAMBRIDGE UNIVERSITY PRESS, 1999.

Wolfram Mathworld http://mathworld.wolfram.com/

George Gasper, Mizan Rahman,
Basic hypergeometric series,
Encyclopedia of Mathematics and its Applications 96,
CAMBRIDGE UNIVERSITY PRESS, 2004.

Victor Kac, Pokman Cheung,
Quantum calculus,
Springer Science & Business Media, 2001.
