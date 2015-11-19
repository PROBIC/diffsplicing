#!/bin/bash

paste t1_r1 base > t1_r1.exp
paste t1_r2 base > t1_r2.exp
paste t1_r3 base > t1_r3.exp

paste t2_r1 base > t2_r1.exp
paste t2_r2 base > t2_r2.exp
paste t2_r3 base > t2_r3.exp

paste t3_r1 base > t3_r1.exp
paste t3_r2 base > t3_r2.exp
paste t3_r3 base > t3_r3.exp

paste t4_r1 base > t4_r1.exp
paste t4_r2 base > t4_r2.exp
paste t4_r3 base > t4_r3.exp

paste t5_r1 base > t5_r1.exp
paste t5_r2 base > t5_r2.exp
paste t5_r3 base > t5_r3.exp

paste t6_r1 base > t6_r1.exp
paste t6_r2 base > t6_r2.exp
paste t6_r3 base > t6_r3.exp

paste t7_r1 base > t7_r1.exp
paste t7_r2 base > t7_r2.exp
paste t7_r3 base > t7_r3.exp

paste t8_r1 base > t8_r1.exp
paste t8_r2 base > t8_r2.exp
paste t8_r3 base > t8_r3.exp

paste t9_r1 base > t9_r1.exp
paste t9_r2 base > t9_r2.exp
paste t9_r3 base > t9_r3.exp

paste t10_r1 base > t10_r1.exp
paste t10_r2 base > t10_r2.exp
paste t10_r3 base > t10_r3.exp

sed 's/\t/ /g' t1_r1.exp > t1_r1.exps
sed 's/\t/ /g' t1_r2.exp > t1_r2.exps
sed 's/\t/ /g' t1_r3.exp > t1_r3.exps

sed 's/\t/ /g' t2_r1.exp > t2_r1.exps
sed 's/\t/ /g' t2_r2.exp > t2_r2.exps
sed 's/\t/ /g' t2_r3.exp > t2_r3.exps

sed 's/\t/ /g' t3_r1.exp > t3_r1.exps
sed 's/\t/ /g' t3_r2.exp > t3_r2.exps
sed 's/\t/ /g' t3_r3.exp > t3_r3.exps

sed 's/\t/ /g' t4_r1.exp > t4_r1.exps
sed 's/\t/ /g' t4_r2.exp > t4_r2.exps
sed 's/\t/ /g' t4_r3.exp > t4_r3.exps

sed 's/\t/ /g' t5_r1.exp > t5_r1.exps
sed 's/\t/ /g' t5_r2.exp > t5_r2.exps
sed 's/\t/ /g' t5_r3.exp > t5_r3.exps

sed 's/\t/ /g' t6_r1.exp > t6_r1.exps
sed 's/\t/ /g' t6_r2.exp > t6_r2.exps
sed 's/\t/ /g' t6_r3.exp > t6_r3.exps

sed 's/\t/ /g' t7_r1.exp > t7_r1.exps
sed 's/\t/ /g' t7_r2.exp > t7_r2.exps
sed 's/\t/ /g' t7_r3.exp > t7_r3.exps

sed 's/\t/ /g' t8_r1.exp > t8_r1.exps
sed 's/\t/ /g' t8_r2.exp > t8_r2.exps
sed 's/\t/ /g' t8_r3.exp > t8_r3.exps

sed 's/\t/ /g' t9_r1.exp > t9_r1.exps
sed 's/\t/ /g' t9_r2.exp > t9_r2.exps
sed 's/\t/ /g' t9_r3.exp > t9_r3.exps

sed 's/\t/ /g' t10_r1.exp > t10_r1.exps
sed 's/\t/ /g' t10_r2.exp > t10_r2.exps
sed 's/\t/ /g' t10_r3.exp > t10_r3.exps

