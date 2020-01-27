function [coherence] = compute_coherence(residual, sens, Ncoils, option)
coherence = zeros(size(residual));

%option == 1 : l1 sum of inner products
%option == 2 : l2 sum of inner products

if option == 1
    for i = 1:Ncoils
        sens_i = squeeze(sens(:,:,i));
        coherence = coherence + sqrt(size(sens,1)*size(sens,2))*abs(Nfft2(conj(sens_i).*residual));
    end
end
   
if option == 2
    for i = 1:Ncoils
        sens_i = squeeze(sens(:,:,i));
        coherence = coherence + (sqrt(size(sens,1)*size(sens,2))*abs(Nfft2(conj(sens_i).*residual))).^2;
    end
    coherence = sqrt(coherence);
end

if option == 3
    coherence = zeros(size(sens,1),size(sens,2));
    for i = 1:size(sens,1)
        for j = 1:size(sens,2)
            coherence(i,j) = compute_projection(residual,sens,i,j);            
        end
    end
end

if option == 4
    for i = 1:Ncoils
    end
end
return;