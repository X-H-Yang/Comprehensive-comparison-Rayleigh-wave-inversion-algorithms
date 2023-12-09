function drawDispersionsCompareField_5(dispersions_R_inverted,...
    f,modes_num_vec,index_vec,myFontSize)


dispersions_R_inverted_cell = cell(1,length(modes_num_vec));
f_cell = cell(1,length(modes_num_vec));

index_start = 1;
for i = 1:1:length(modes_num_vec)
    temp = index_vec{i};
    f_cell{i} = f(temp(1):temp(2));
    dispersions_R_inverted_cell{i} = dispersions_R_inverted(...
        index_start:index_start+length(f_cell{i})-1);
    index_start = index_start + length(f_cell{i});
end

my_linewidth = 1.8;
plot(f_cell{1},dispersions_R_inverted_cell{1},'Color',[255 0 0]/255,...
    'LineWidth',my_linewidth);
if length(modes_num_vec) > 1
    for i = 1:1:length(modes_num_vec)
        plot(f_cell{i},dispersions_R_inverted_cell{i},'Color',...
            [255 0 0]/255,'LineWidth',my_linewidth);
    end
end


xlabel('Frequency [Hz]','FontSize',myFontSize);
ylabel('Phase velocity [m/s]','FontSize',myFontSize);
set(gca,'FontName','Times New Roman','FontSize',myFontSize);
end
