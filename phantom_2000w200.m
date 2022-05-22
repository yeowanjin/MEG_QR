%% code/instructions to reproduce blue curve of Figure 2 (200 nAm signal reconstruction with 2000 nAm basis)

% First, load the 32 dipole datasets for both 200 nAm and 2000 nAm 
% 306 sensor Elekta phantom measurements
% Here, the data matrices are stored respectively in
% ``data_a200_dip_306{dipole#}'' and % ``data_a2000_dip_306{dipole#}''
% 
% Also, save and load QR pivots/selected sensor positions from
% the m file, phantom_2002000.m.
% Here, this is stored as pivots_dip{dipole#,#sensors}
%
% The dataset we used was obtained from Brainstorm and can be downloaded
% from http://neuroimage.usc.edu/brainstorm

for loopy = 1:32
    
    % obtain SVD basis for 2000 nAm signal
    [U,S,V]=svd(data_a2000_dip_306{loopy},'econ');

    % use 2000 nAm SVD basis and p number of 2000 nAm QR pivots to 
    % reconstruct full-state 200 nAm sigal with p measurements
    for p = 1:306
        pivots = pivots_dip{loopy,p};
        UU = U(:,1:p);
        
        % permutation matrix
        Pdata = zeros(p,size(UU,1));
        for j = 1:p
           Pdata(j,pivots(j)) = 1; 
        end
        
        Xrecon{loopy,p} = real(UU*pinv(Pdata*UU)*Pdata*data_a200_dip_306);
        relerr_306{loopy}(p) = norm(data_a200_dip_306 - Xrecon{loopy,p})/norm(data_a200_dip_306);
        
    end    
end

% average over all 32 datasets
ave_306 = 0;
for l = 1:loopy
   ave_306 = ave_306 + relerr_306{l};
end
ave_306 = ave_306(1:100)/32;

% plot blue for Figure 2
figure
semilogy(ave_306,'b','LineWidth',1.5)
legend('using 2000 nAm basis')