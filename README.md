# q-special-functions

programs for q-special functions and q-series

Thank you for checking my programs.

Programs in this repository are made to calculate q-special functions and q-series with guaranteed accuracy. All q-special functions are computed by calculating infinite q-Pochhammer symbols and q-hypergeometric series with guaranteed accuaracy. 

If you use my programs in publications, please include the following reference.

Daisuke Kanaizumi, Verified numerical computation methods for q-special functions and q-series,
Bachelor thesis, Department of Applied Mathematics, Waseda University (Japanese document)

You do NOT have to inform your usage to me, but I am very glad to hear how my programs helped you.

I have a plan to translate my thesis in English, but this may take more time. I apologize for the inconvenience.

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

If you want to use my Python programs, please install pint.

pint-interval type Library for Python 3.x
https://github.com/o108minmin/pint

DISCLAIMER: Extensive tests have been performed to ensure reliability of the algorithms (under Windows 10 with Visual Studio 2015, 2017 and Ubuntu 14.04). However, neither an error-free processor nor an error-free program can be guaranteed.

(q-special functions in this repository)

elliptic gamma function

little a-Jacobi function 

Heine hypergeometric function

q-Airy functions (HKW, KMNOY, IKMMS, Ramanujan)

q-Appell function (\Phi_1)

q-Barnes G function

q-Bessel functions (Jackson, Hahn-Exton, little, modified)

q-beta function

q-beta integrals (Askey-Wilson, Nassrallah-Rahman, Ramanujan)

q-Bi function 

q-digamma

q-logarithm

q-gamma functions (Jackson, Moak)

quantum dilogarithm Li_2(x;q) (defined by Kirillov)

quantum polylogarithm

Ramanujan psi sum

Ramanujan theta function

symmetric q-gamma, symmetric q-beta

(q-series in this repository)

Bailey Mod9

Dyson Mod27

Gessel-Stanton

Gollnitz-Gordon

Jackson-Slater

Lebesgue

Rogers Mod14

Rogers-Ramanujan

Rogers-Selberg

(References)

S.M. Rump: INTLAB - INTerval LABoratory. 
In Tibor Csendes, editor, Developments in Reliable Computing, pages 77-104. 
Kluwer Academic Publishers, Dordrecht, 1999.

George. E. Andrews, Mircea Merca, 
The truncated pentagonal number theorem,
Journal of Combinatorial Theory, Series A,

M. A. Olshanetsky, V. B. K. Rogov,
The Modified q-Bessel Functions and the q-Bessel Macdonald Functions

Mourad E. H. Ismail,
CLASSICAL AND QUANTUM ORTHOGONAL POLYNOMIALS IN ONE VARIABLE,
Encyclopedia of Mathematics and its Applications 98,
CAMBRIDGE UNIVERSITY PRESS

George E. Andrews, Richard Askey, Ranjan Roy,
SPECIAL FUNCTIONS,
Encyclopedia of Mathematics and its Applications 71,
CAMBRIDGE UNIVERSITY PRESS

Don Zagier, The Dilogarithm Function

Shin Isojima, Ultradiscrete limit of Bessel function type solutions of the Painleve III equation

Yousuke Ohyama, Particular solutions of q-Painleve equations and q-hypergeometric equations

Bouzeffour, New Addition Formula for the Little q-Bessel Functions, arXiv, 2013

Ruiming Zhang, Plancherel-Rotach asymptotics for certain basic hypergeometric series, Advances in Mathematics, 2008

Fredrik Johansson, Computing hypergeometric functions rigorously, arXiv, 2016

Brahim and Sidomou, On Some Symmetric q-Special Functions, 2013

Wolfram Mathworld http://mathworld.wolfram.com/

Savage, Sills, On an identity of Gessel and Stanton and new little Gollnitz identities, 2009 

M. A. Bershtein, A. I. Shechechkin, q-deformed Painlev\`e \tau function and q-deformed conformal blocks, arXiv, 2016

Isojima, Konno, Mimura, Murata, Satsuma, Ultradiscrete Painlev\`e II equation and a special function solution, JOURNAL OF PHYSICS A: MATHEMATICAL AND THEORETICAL, 2011

Ahmed Fitouhi, Fethi Bouzeffour, Wafa Binous, Expansion and asymptotic in terms of basic Bessel functions, Applied Mathematics and Computation 188 (2007) 2034–2044

H. T Koelink, Hansen-Lommel Orthogonality Relations for Jackson`s q-Bessel functions, Journal of Mathematical Analysis and Applications 175, 425-437 (1993)

A. B. Olde Daalhuis, Asymptotic Expansions for q-Gamma, q-Exponential and q-Bessel Functions, Journal of Mathematical Analysis and Applications 186, 896-913 (1994)

Mourad E. H. Ismail, Changgui Zhang, Zeros of entire functions and a problem of Ramanujan, Advances in Mathematics 209 (2007) 363–380
