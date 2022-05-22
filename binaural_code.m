%% code/instructions to reproduce Figure 9

% First, load the binaural MEG dataset/matrix 
% Here, the data matrix is stored in "data_meg".
%
% The dataset we used was obtained from Brainstorm and can be downloaded
% from http://neuroimage.usc.edu/brainstorm

% obtain SVD basis U
[U,S,V] = svd(data_meg,'econ');
 
% perform the QR algorithm so that 1:all sensors are selected
k=1;
for p = 1:size(data_meg,1) 
    r = p; % here we set number of modes = number of sensors selected
    if r < size(U,2)
        UU = U(:,1:r);
    else
        UU = U;
    end
    
    % QR to obtain pivots
    if (p <= r)
        [Q,R,pivots] = qr(UU','vector');
    elseif (p > r)
        [Q,R,pivots] = qr(UU*UU','vector');
    end
    pivots1{k} = pivots(1:p);

    % permutation matrix
    Pdata{k} = zeros(p,size(UU,1)); 
    for j = 1:p
       Pdata{k}(j,pivots1{k}(j)) = 1; 
    end
    
    % QR reconstructions using 1:all sensors.
    % We used the reconstruction with 30 sensors, i.e. Xrecon{30}
    % Then, it was fed back to Brainstorm to plot Figure 9(b)
    Xrecon{k} = real(UU*pinv(Pdata{k}*UU)*Pdata{k}*data_meg);
    relerr(k) = norm(data_meg - Xrecon{k})/norm(data_meg);

    k=k+1;
end

% singular values in %
S_perc = 100*diag(S)./sum(diag(S)); % singular values in %

% plot Figure 9(a)
figure(1)
colororder({'k','r'})
yyaxis left
semilogy(relerr(1:100),'k','LineWidth',1.5)
xlabel('Total # sensors/modes used')
ylabel('Relative error')
yyaxis right
semilogy(S_perc(1:100),'r--','LineWidth',1.5)
ylabel('Singular value (%)')
legend('QR selected sensors','Singular value')
title('Relative error (binaural signal)')
