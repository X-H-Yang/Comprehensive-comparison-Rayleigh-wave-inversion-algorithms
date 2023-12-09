function drawProfile_fun(Y_hat,h_true,Vs_true)

max_depth = 100;

layers_num = length(h_true) + 1;
h_inverted = Y_hat(1:layers_num-1);
Vs_inverted = Y_hat(layers_num:end);

% draw inverted profile
h_inverted_plot = [0 repelem(cumsum(h_inverted),2) max_depth];
Vs_inverted_plot = repelem(Vs_inverted,2);


% draw true profile
h_true_plot = [0 repelem(cumsum(h_true),2) max_depth];
Vs_true_plot = repelem(Vs_true,2);


plot(Vs_true_plot,h_true_plot,'k','Linewidth',1.8);
hold on
plot(Vs_inverted_plot,h_inverted_plot,'r','Linewidth',1.8);


set(gca,'xaxislocation','top');

% box on
set(gca,'ydir','reverse')


