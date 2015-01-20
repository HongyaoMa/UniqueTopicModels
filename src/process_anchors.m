function anchor_inds_out = process_anchors(anchors_in)

        k = length(anchors_in);
        
        anchor_modk = mod(anchors_in, k);
        anchor_modk(anchor_modk == 0) = k;
        [~, ind_sort] = sort(anchor_modk);
        anchor_inds_out = anchors_in(ind_sort);
end