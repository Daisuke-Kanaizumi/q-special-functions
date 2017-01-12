# q-special-functions

programs for q-special functions

Thank you for checking my programs.

Programs in this repository are made to calculate q-special functions with guaranteed accuracy. All q-special functions are computed by calculating infinite q-Pochhammer symbols with guaranteed accuaracy. Japanese Documentations (which include details for the methods that I used) are uploaded to Speakerdeck.

https://speakerdeck.com/daisuke15

I will translate them to English (as soon as possible). Please wait for a while.

If you use my programs in publications, please include the following reference.

Daisuke Kanaizumi, Applications of verified numerical computation methods for q-special functions, q-series, hypergeometric functions and error functions,
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

DISCLAIMER: Extensive tests have been performed to ensure reliability of the algorithms (under Windows 10 with Visual Studio 2015, and Linux). However, neither an error-free processor nor an error-free program can be guaranteed.

(Functions in this repository)

Bailey-Daum sum 

q-Bessel functions (Jackson, Hahn-Exton)

q-beta function

q-beta integrals (Askey-Wilson, Nassrallah-Rahman, Ramanujan)

q-gamma functions (Jackson, Moak)

q-Ramanujan function

quantum dilogarithm Li_2(x;q) (defined by Kirillov)

Ramanujan psi sum

Ramanujan theta function

(References)

S.M. Rump: INTLAB - INTerval LABoratory. 
In Tibor Csendes, editor, Developments in Reliable Computing, pages 77-104. 
Kluwer Academic Publishers, Dordrecht, 1999.

George. E. Andrews, Mircea Merca, 
The truncated pentagonal number theorem,
Journal of Combinatorial Theory, Series A

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

Wolfram Mathworld http://mathworld.wolfram.com/
