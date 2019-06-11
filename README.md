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

ACKNOWLEDGEMENTS: The author thanks to 
* Masahide Kashiwagi (http://verifiedby.me/) & Kouta Sekine (https://271.jp/index.php) for their  technical supports and advices.
* Everyone at the Waseda Soliton Lab for their kindness, friendship and advices.
* Kenichi Maruno for giving guidances about basic hypergeometric series and providing documents/computers.

(q-special functions in this repository)

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

(References for Validated Numerics and Interval Arithmetic)

S.M. Rump: INTLAB - INTerval LABoratory. 
In Tibor Csendes, editor, Developments in Reliable Computing, pages 77-104. 
Kluwer Academic Publishers, Dordrecht, 1999.

Kulisch, U., Lohner, R., & Facius, A. (Eds.). (2012). Perspectives on Enclosure Methods. Springer Science & Business Media.

Tucker, W. (2011). Validated numerics: a short introduction to rigorous
computations. Princeton University Press.

Mayer, G. (2017). Interval analysis: and automatic result verification (Vol. 65).
Walter de Gruyter GmbH & Co KG.

Alefeld, G., & Herzberger, J. (2012). Introduction to Interval Computation.
Academic Press.

Moore, R. E., Kearfott, R. B., & Cloud, M. J. (2009). Introduction to Interval
Analysis (Vol. 110). SIAM.

Jaulin, L., Kieffer, M., Didrit, O., & Walter, E. (2001). Applied Interval Analysis. Springer, London.

(References for q-Special Functions)

George Gasper, Mizan Rahman,
Basic hypergeometric series,
Encyclopedia of Mathematics and its Applications 96,
CAMBRIDGE UNIVERSITY PRESS, 2004.

Mourad E. H. Ismail,
CLASSICAL AND QUANTUM ORTHOGONAL POLYNOMIALS IN ONE VARIABLE,
Encyclopedia of Mathematics and its Applications 98,
CAMBRIDGE UNIVERSITY PRESS, 2005.

George E. Andrews, Richard Askey, Ranjan Roy,
SPECIAL FUNCTIONS,
Encyclopedia of Mathematics and its Applications 71,
CAMBRIDGE UNIVERSITY PRESS, 1999.

Wolfram Mathworld http://mathworld.wolfram.com/

DLMF: Digital Library of Mathematical Functions https://dlmf.nist.gov/

Ahmad, B., Ntouyas, S., Tariboon, J.
Quantum Calculus (2016): New Concepts, Impulsive IVPs and BVPs, Inequalities. World Scientific, Singapore. 

Victor Kac, Pokman Cheung,
Quantum calculus,
Springer Science & Business Media, 2001.

Ernst, T. (2012). A Comprehensive Treatment of q-calculus. Springer
Science & Business Media.

Ernst, T. (2000). The History of q-Calculus and a New Method. Department
of Mathematics, Uppsala University.

Koekoek, R., Lesky, P. A., & Swarttouw, R. F. (2010). Hypergeometric
Orthogonal Polynomials and their q-Analogues. Springer Science & Business
Media.

Andrews, G. E. (1986). q-Series: Their Development and Application in
Analysis, Number Theory, Combinatorics, Physics, and Computer Algebra
(No. 66). American Mathematical Society.

Exton, H. (1983). q-Hypergeometric Functions and Applications. Horwood.

Slater, L. J. (1966). Generalized Hypergeometric Functions, CAMBRIDGE UNIVERSITY PRESS.

Bailey, W. N. (1964). Generalized Hypergeometric Series, CAMBRIDGE UNIVERSITY PRESS.

Swarttouw, R. F. (1992), The Hahn-Exton q-Bessel Function, PhD Thesis,
Delft Technical University.
https://repository.tudelft.nl/islandora/object/uuid%3A0bfd7d96-bb29-4846-8023-6242ce15e18f?collection=research

Koelink, E. (2018). q-Special Functions, Basic Hypergeometric Series and
Operators. arXiv preprint arXiv:1808.03441.

Koornwinder, T. H. (2005). q-Special Functions, an Overview. arXiv preprint
math/0511148.

(References about numerical methods for mathematical functions)

Muller, (2016). Elementary Functions -Algorithms and Implementation (3rd ed.). Birkhauser.

Gil, A., Segura, J., & Temme, N. M. (2007). Numerical methods for special functions (Vol. 99). SIAM.

Marcellán, F. (2006). Orthogonal polynomials and special functions: computation and applications (No. 1883). Springer Science & Business Media.
