echo "
1 |  3.87 | x1 3.87 (1.01, -0.50)
 2234   | 1.000 (1.000, 1.000) | ../dat/modeling3_310r.dat (0.000)
"
tar xzf ../../pdb/modeling3.tar.gz modeling3_310r.pdb
mkdir 1state
mv modeling*.pdb 1state
#exit -1

echo "
1 |  3.67 | x1 3.67 (1.01, -0.50)
 2711   | 0.449 (0.501, 0.123) | ../dat/modeling3_740r.dat (0.013)
 3886   | 0.551 (0.661, 0.125) | ../dat/modeling4_899r.dat (0.035)
"
tar xzf ../../pdb/modeling3.tar.gz modeling3_740r.pdb
tar xzf ../../pdb/modeling4.tar.gz modeling4_899r.pdb
mkdir 2states
mv modeling*.pdb 2states
#exit -1


echo "
1 |  3.64 | x1 3.64 (1.01, -0.50)
 1417   | 0.362 (0.338, 0.051) | ../dat/modeling2_476r.dat (0.532)
 3886   | 0.388 (0.545, 0.075) | ../dat/modeling4_899r.dat (0.362)
 4963   | 0.250 (0.363, 0.138) | ../dat/modeling5_968r.dat (0.234)
"
tar xzf ../../pdb/modeling2.tar.gz modeling2_476r.pdb
tar xzf ../../pdb/modeling4.tar.gz modeling4_899r.pdb
tar xzf ../../pdb/modeling5.tar.gz modeling5_968r.pdb
mkdir 3states
mv modeling*.pdb 3states
exit -1


echo "
1 |  1.13 | x1 1.13 (1.01, -0.50)
  100   | 0.325 (0.334, 0.021) | ../dat2012_open/modeling21_1090r.dat (0.229)
 3310   | 0.324 (0.376, 0.015) | ../dat2012_open/modeling21_3980r.dat (0.229)
12460   | 0.258 (0.228, 0.022) | ../dat2012_open/modeling23_3214r.dat (0.229)
39241   | 0.093 (0.095, 0.014) | ../dat2012_closed/modeling33_4818r.dat (0.001)
"
tar xzf ../../pdb/modeling21.tar.gz modeling21_1090r.pdb
tar xzf ../../pdb/modeling21.tar.gz modeling21_3980r.pdb
tar xzf ../../pdb/modeling23.tar.gz modeling23_3214r.pdb
tar xzf ../../pdb/modeling33.tar.gz modeling33_4818r.pdb
mkdir 4states
mv modeling*.pdb 4states
#exit -1


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

