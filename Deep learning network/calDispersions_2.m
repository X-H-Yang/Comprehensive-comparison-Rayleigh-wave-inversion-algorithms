function dispersion_R_inverted = calDispersions_2(Vs_profile,Vp,den,f,modes_num_vec,index_vec)

% Acknowledgement: The forward modeling program used to generate 
%                  theoretical Rayleigh wave dispersion curves in this 
%                  study was obtained from the  website 
%                  (https://github.com/eespr/MuLTI) provided by 
%                  Killingbeck et al. (2018)
%
% Killingbeck et al. (2018): Killingbeck, S. F., Livermore, P. W., 
%                            Booth, A. D., & West, L. J. (2018). Multimodal 
%                            layered transdimensional inversion of seismic 
%                            dispersion curves with depth constraints. 
%                            Geochemistry, Geophysics, Geosystems, 19(12), 
%                            4957-4971.

df = f(2)-f(1);

pointNum = 0;
for jj = 1:1:length(modes_num_vec)
    temp = modes_num_vec(jj);
    point_temp = index_vec{temp};
    star_point = point_temp(1);
    end_point = point_temp(2);
    
    pointNum = pointNum + (end_point-star_point) + 1;

end

Vs_profile = Vs_profile(:)';
Vp = Vp(:)';
den = den(:)';

layers_num = floor(length(Vs_profile)/2) + 1;

h = Vs_profile(1:layers_num-1);
h2 = [h 0];
Vs = Vs_profile(layers_num:end);

dispersion_R_inverted = [];

if sum(Vs_profile<=0)
    dispersion_R_inverted = zeros(1,pointNum);
    
else
    try
        out = gpdc(h2,Vp,Vs,den,'fV',f);
        out2 = rdivide(1, out(:, 2:end));
        for jj = 1:1:length(modes_num_vec)
            temp = modes_num_vec(jj);
            point_temp = index_vec{temp};
            star_point = point_temp(1);
            end_point = point_temp(2);
            dispersion_R_inverted = [dispersion_R_inverted out2(star_point:end_point,temp)'];

        end
        dispersion_R_inverted(isnan(dispersion_R_inverted)) = 0;
    catch
        dispersion_R_inverted = zeros(1,pointNum);
    end
end

end

