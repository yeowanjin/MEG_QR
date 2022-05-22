%% code/instructions to reproduce Figures 3-8 and part of Figure 2

% First, load the 32 dipole datasets/matrices for the 200 nAm or 2000 nAm 
% 306 sensor Elekta phantom measurements
% Here, the data matrices are stored in ``data1{dipole#}''.
%
% The dataset we used was obtained from Brainstorm and can be downloaded
% from http://neuroimage.usc.edu/brainstorm

% QR for each dataset of the 32 dipolar sources
for loopy = 1:32
    
    % obtain SVD basis U
    [U,S,V] = svd(data1{loopy},'econ');
    
    % perform the QR algorithm so that 1:all sensors are selected
    k = 1;
    for p = 1:size(data1{loopy},1) % p = # sensors, r = # modes
        r = p; % here we set number of modes = number of sensors selected
        if r < size(U,2)
            UU = U(:,1:r);
        else
            UU = U;
        end
        
        % SVD reconstructions using 1:all sensors.
        % We used the reconstruction with 2 sensors, i.e.
        % Xrecon{dipole#,2}
        Xrecon_svd{loopy,k} = UU*S(1:size(UU,2),1:size(UU,2))*V(:,1:size(UU,2))';
        
        % QR to obtain pivots
        if (p <= r)
            [Q,R,pivots] = qr(UU','vector');
        elseif (p > r)
            [Q,R,pivots] = qr(UU*UU','vector');
        end
        pivots1{loopy,k} = pivots(1:p); 

        % permutation matrix
        Pdata{k} = zeros(p,size(UU,1));
        for j = 1:p
           Pdata{k}(j,pivots1{loopy,k}(j)) = 1; 
        end
        
        % QR reconstructions using 1:all sensors.
        % We used the reconstruction with 2 and 35 sensors, i.e. 
        % Xrecon{dipole#,2} and Xrecon{dipole#,35}
        % Then, they were fed back to Brainstorm to obtain field plots for
        % Figures 4, 5 and data for Figures 6, 7, 8
        Xrecon{loopy,k} = real(UU*pinv(Pdata{k}*UU)*Pdata{k}*data1);
        relerr_306{loopy}(k) = norm(data1 - Xrecon{loopy,k})/norm(data1);
                        
        k = k + 1;
    end

    % singular values in %
    S_perc{loopy} = 100*diag(S)/sum(diag(S));

end

% average over all 32 datasets
ave_306 = 0;
ave_S_perc = 0;
for l = 1:loopy
   ave_306 = ave_306 + relerr_306{l};
   ave_S_perc = ave_S_perc + S_perc{l};
end
ave_306 = ave_306/32;
ave_S_perc = ave_S_perc/32;

% plot black, red curves for Figure 2
figure(1)
colororder({'k','r'})
yyaxis left
semilogy(ave_306(1:100),'k','LineWidth',1.5)
xlabel('Total # sensors/modes used')
ylabel('Relative error')
yyaxis right
semilogy(ave_S_perc(1:100),'r--','LineWidth',1.5)
ylabel('Singular value (%)')
legend('QR selected sensors','Singular value')


% To recreate Figure 3, plot imagesc of the following:
% First column: data1{dipole#}
% Second column: Xrecon{dipole#,2}, Xrecon_svd{loopy,2}
% Third column: U(:,1)*S(1,1)*V(:,1)', U(:,2)*S(2,2)*V(:,2)' for dipole#