function D = Hex8_TangentD(E_tan, nu)
coef = E_tan / ((1 + nu) * (1 - 2*nu));

D = coef * [1-nu   nu     nu     0              0              0;
            nu     1-nu   nu     0              0              0;
            nu     nu     1-nu   0              0              0;
            0      0      0      (1-2*nu)/2     0              0;
            0      0      0      0              (1-2*nu)/2     0;
            0      0      0      0              0              (1-2*nu)/2];
end