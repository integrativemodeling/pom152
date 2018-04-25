echo "
1 |  1.59 | x1 1.59 (1.01, -0.50)
 4958   | 1.000 (1.000, 1.000) | ../dat/modeling5_963r.dat (0.000)
"
tar xzf ../../pdb/modeling5.tar.gz modeling5_963r.pdb
mkdir 1state
mv modeling*.pdb 1state
#exit -1

echo "
1 |  1.31 | x1 1.31 (1.01, -0.48)
  232   | 0.660 (0.746, 0.063) | ../dat/modeling1_309r.dat (0.119)
 3886   | 0.340 (0.681, 0.123) | ../dat/modeling4_899r.dat (0.026)
"
tar xzf ../../pdb/modeling1.tar.gz modeling1_309r.pdb
tar xzf ../../pdb/modeling4.tar.gz modeling4_899r.pdb
mkdir 2states
mv modeling*.pdb 2states
#exit -1


echo "
1 |  1.28 | x1 1.28 (1.01, -0.46)
  876   | 0.148 (0.175, 0.023) | ../dat/modeling1_88r.dat (0.296)
 3886   | 0.636 (0.655, 0.097) | ../dat/modeling4_899r.dat (0.301)
 3892   | 0.216 (0.178, 0.063) | ../dat/modeling4_903r.dat (0.006)
"
tar xzf ../../pdb/modeling1.tar.gz modeling1_88r.pdb
tar xzf ../../pdb/modeling4.tar.gz modeling4_899r.pdb
tar xzf ../../pdb/modeling4.tar.gz modeling4_903r.pdb
mkdir 3states
mv modeling*.pdb 3states
#exit -1


echo "
1 |  1.25 | x1 1.25 (1.01, -0.46)
  876   | 0.088 (0.166, 0.028) | ../dat/modeling1_88r.dat (0.844)
 1020   | 0.570 (0.431, 0.233) | ../dat/modeling2_118r.dat (0.317)
 1822   | 0.100 (0.098, 0.010) | ../dat/modeling2_840r.dat (0.146)
 3771   | 0.243 (0.308, 0.083) | ../dat/modeling4_795r.dat (0.137)
"
tar xzf ../../pdb/modeling1.tar.gz modeling1_88r.pdb
tar xzf ../../pdb/modeling2.tar.gz modeling2_118r.pdb
tar xzf ../../pdb/modeling2.tar.gz modeling2_840r.pdb
tar xzf ../../pdb/modeling4.tar.gz modeling4_795r.pdb
mkdir 4states
mv modeling*.pdb 4states
exit -1


echo "
7 |  1.12 | x1 1.12 (1.01, -0.50)
 4817   | 0.438 (0.351, 0.052) | ../dat2012_open/modeling21_836r.dat (0.114)
 9921   | 0.111 (0.142, 0.034) | ../dat2012_open/modeling22_92r.dat (0.114)
12507   | 0.231 (0.240, 0.020) | ../dat2012_open/modeling23_3257r.dat (0.114)
17528   | 0.146 (0.307, 0.056) | ../dat2012_open/modeling24_3276r.dat (0.967)
48476   | 0.074 (0.074, 1.000) | ../dat2012_closed/modeling35_4129r.dat (0.001)
"
tar xzf ../../pdb/modeling21.tar.gz modeling21_836r.pdb
tar xzf ../../pdb/modeling22.tar.gz modeling22_92r.pdb
tar xzf ../../pdb/modeling23.tar.gz modeling23_3257r.pdb
tar xzf ../../pdb/modeling24.tar.gz modeling24_3276r.pdb
tar xzf ../../pdb/modeling35.tar.gz modeling35_4129r.pdb
mkdir 5states
mv modeling*.pdb 5states
#exit -1

