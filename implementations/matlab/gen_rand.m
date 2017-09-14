function value = gen_rand(LB, HB)
    range = HB - LB + 1;
    value = rem(commonRandom(), range) + LB;
end
