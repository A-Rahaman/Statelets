function slices = Slicer(cand_len,TotalWindows)
slices = {};
sl_count = 1;
lm = 1;
   while (sl_count<=TotalWindows)
        %fprintf("Slicing\n\n")
        if((sl_count+cand_len)<=TotalWindows)
        slices{lm} = sl_count:(sl_count+cand_len)-1;
        sl_count = sl_count+cand_len;
        %fprintf("Slice = %d Sl_count = %d\n\n",lm,sl_count);
        else
            slices{lm} = sl_count:TotalWindows;
            sl_count = sl_count+length(slices{lm});
            %fprintf("Last Slice %d\n",lm);
        end
        lm = lm + 1;
   end
end