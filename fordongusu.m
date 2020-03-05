for i=1:420
deneme(t_peak_start_pos(i):t_peak_final_pos(i))=ecgfiltered_last(t_peak_start_pos(i):t_peak_final_pos(i));
end

en1=sum(sum(sc(:,1:60e3)))
en2=sum(sum(sc(:,60e3+1:120e3)))
en3=sum(sum(sc(:,120e3+1:180e3)))