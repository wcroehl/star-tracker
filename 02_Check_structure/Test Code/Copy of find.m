function [rtn] = find(t, crf_count, obs_count, obs_ra, obs_dec, obs_mag,...
    id_count, id_star, id_body, vici, Mp, outdist, t_start, t_end)
    rtn = 0;
    sig = 0;
    count = sig;
    if t > t_start && t < t_end
        fprintf("FIND: obs_count=%3d  id_count=%3d\n", obs_count, id_count);
    end
    
    for i=0:1:obs_count-1        
        find_num = 40; %sizeof(int) = 4
        
        if t > t_start && t < t_end 
            fprintf(stdout, "Here1: *id_count=%3d\n", id_count) ; 
        end
        
        iptr = find_num;
        
        if t > t_start && t < t_end
            fprintf(stdout, "FIND: t=%15.6f c_star[%2d]->mag = %7.2f\n", t, i, c_star(i)); 
        end
        
        for j=0:1:crf_count-1
            for k=0:1:id_count-1 
                if crf(j) == id_star(k)
                    sig = 8 ;
                end
            end
            if sig == 8
                sig = 0;
                continue;
            end
            
            vec_dot = c_star(i) * crf(j) + c_star(i) * crf(j)...
            + c_star(i) * crf(j);
            if vec_dot > 1.0
                vec_dot =  1.0 ;
            end
            if vec_dot < -1.0
                vec_dot = -1.0 ;
            end
            arc_dist = acos(vec_dot)*180.0/PI*3600.0 ;
            mag_diff = c_star(i)-crf(j);
            
            if (t > t_start && t < t_end)     
                fprintf(stdout, "FIND: t=%15.6f  crf[%2d]->num=%5d  mag=%5.2f ", t, j, crf(j), crf(i)) ;
                fprintf(stdout, "arc_dist (arcsec) =%10.4f\n", arc_dist) ;
            end
            
            if(abs(arc_dist) < vici && abs(mag_diff) < mTOL0) %found one candidate: was VICINITY
                iptr = stars(crf(j)).cat_num ;
                fprintf(outdist, "1  %15.6f %5d %15.8f \n", t, iptr, arc_dist);
                count=count+1;
                iptr=iptr+1;
            end
        end
        find_num2 = 50; %sizeof(int)=4
        switch (count)
            case 0 %False observation or New observation
                count0 = 0;
                chosen_num = crf(floor(crf_count/2.0));
                
                rtn = new_obs(t, t_start, t_end, count0, find_num2, chosen_num, ...
                    c_star+i, obs_ra+i, obs_dec+i, obs_mag+i,...  
                    id_star, id_count, vici, Mp, outdist);
                
                if rtn ~= -1
                    fprintf(stdout,"starID_dm()::new_obs() finished abnormally.\n");
                    return;
                end
                switch count0
                    case 0
                        sig = 0 ; % No identified star is reported */
                        break ;
                    case 1 
                        sig = 10; %New star : *find_num2  */
                        chosen_num = find_num2 ;
                        break ;
                    otherwise
                        sig = 20; %Multiple candidates : to choose one */
                        chosen_num = multi_obs(t, (c_star+i), find_num2, count0) ; 
                end
                break;
                
            case 1
                sig = 1;
                chosen_num = find_num;
                break;
                
            otherwise
                sig = 2;
                chosen_num = multi_obs(t, (c_star+i), find_num, count);
        end
        if(sig ~= 0)
            id_star(id_count) = chosen_num ;
            id_count = id_count+1;
            id_body(id_count) = i ;  %for c_star   
        end
        count = 0 ;    
    end
    return;
end