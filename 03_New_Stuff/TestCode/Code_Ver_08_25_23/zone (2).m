function return_value = zone(dec)     % find the zone where the BD is pointing %
 
    if (dec >= 80.0)                 return_value=1 ; return ; end
    if (dec <  80.0 && dec >=  75.0) return_value=2 ; return ; end
    if (dec <  75.0 && dec >=  70.0) return_value=3 ; return ; end
    if (dec <  70.0 && dec >=  65.0) return_value=4 ; return ; end
    if (dec <  65.0 && dec >=  60.0) return_value=5 ; return ; end
    if (dec <  60.0 && dec >=  55.0) return_value=6 ; return ; end
    if (dec <  55.0 && dec >=  50.0) return_value=7 ; return ; end
    if (dec <  50.0 && dec >=  45.0) return_value=8 ; return ; end
    if (dec <  45.0 && dec >=  40.0) return_value=9 ; return ; end
    if (dec <  40.0 && dec >=  35.0) return_value=10 ; return ; end
    if (dec <  35.0 && dec >=  30.0) return_value=11 ; return ; end
    if (dec <  30.0 && dec >=  25.0) return_value=12 ; return ; end
    if (dec <  25.0 && dec >=  20.0) return_value=13 ; return ; end
    if (dec <  20.0 && dec >=  15.0) return_value=14 ; return ; end
    if (dec <  15.0 && dec >=  10.0) return_value=15 ; return ; end
    if (dec <  10.0 && dec >=   5.0) return_value=16 ; return ; end
    if (dec <   5.0 && dec >=   0.0) return_value=17 ; return ; end
    if (dec <   0.0 && dec >=  -5.0) return_value=18 ; return ; end
    if (dec <  -5.0 && dec >= -10.0) return_value=19 ; return ; end
    if (dec < -10.0 && dec >= -15.0) return_value=20 ; return ; end
    if (dec < -15.0 && dec >= -20.0) return_value=21 ; return ; end
    if (dec < -20.0 && dec >= -25.0) return_value=22 ; return ; end
    if (dec < -25.0 && dec >= -30.0) return_value=23 ; return ; end
    if (dec < -30.0 && dec >= -35.0) return_value=24 ; return ; end
    if (dec < -35.0 && dec >= -40.0) return_value=25 ; return ; end
    if (dec < -40.0 && dec >= -45.0) return_value=26 ; return ; end
    if (dec < -45.0 && dec >= -50.0) return_value=27 ; return ; end
    if (dec < -50.0 && dec >= -55.0) return_value=28 ; return ; end
    if (dec < -55.0 && dec >= -60.0) return_value=29; return ; end
    if (dec < -60.0 && dec >= -65.0) return_value=30 ; return ; end
    if (dec < -65.0 && dec >= -70.0) return_value=31 ; return ; end
    if (dec < -70.0 && dec >= -75.0) return_value=32 ; return ; end
    if (dec < -75.0 && dec >= -80.0) return_value=33 ; return ; end
    if (dec < -80.0 )                return_value=34 ; return ; end

return_value = -1;
end