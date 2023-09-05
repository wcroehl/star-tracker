function starcell2(fn_dir)
    %arrays
    cellid = zeros(1,3);
    numstars = zeros(1,8);
    starnum = zeros(1,7);
    
    for i = 0:1:33
        fscanf(fn_dir, "%s %s", cellid, numstars);
        scell2(i).cell_id = atoi(cellid);
        scell2(i).num_stars = atoi(numstars);
        scell2(i).star_num = 4*scell(i).numstars;
        
        if scell(i).star_num == NULL
            fprintf(stdout, " FAIL\n");
        end
        for j = 0:1:scell2(i).num_stars
            fscanf(fn_dir, "%s", starnum);
            scell2(i).star_num(j) = atoi(starnum);
        end
    end
    
end