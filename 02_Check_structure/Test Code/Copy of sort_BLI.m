function BLI =sort_BLI(ij, BLI) % sorting BLI from smallest ang. dist.
for i=1:1:ij
    for j=i+1:1:ij % smaller BLI first */
        if(BLI(j).IBL > BLI(i).IBL) 
            temporary.starnum = BLI(i).starnum ;
            temporary.IBL     = BLI(i).IBL     ;
            temporary.L       = BLI(i).L    ;
            temporary.mag     = BLI(i).mag     ;

            BLI(i).starnum = BLI(j).starnum ;
            BLI(i).IBL     = BLI(j).IBL     ;
            BLI(i).L       = BLI(j).L    ;
            BLI(i).mag     = BLI(j).mag     ;

            BLI(j).starnum = temporary.starnum ;
            BLI(j).IBL     = temporary.IBL     ;
            BLI(j).L       = temporary.L       ;
            BLI(j).mag     = temporary.mag     ;
        end
    end
end
       
