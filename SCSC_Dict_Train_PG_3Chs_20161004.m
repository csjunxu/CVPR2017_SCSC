% semi-coupled weighted sparse coding
clear;
addpath('Data');
addpath('SPAMS'); 
addpath('Utilities');

load Data/params_PG.mat;
task = 'BID';

for j = 1:par.Patch_Channel 
    modelname = sprintf('../DSCDL_BID/Data/GMM_PG_%d_10_8x8_64_20161003T094301.mat',j);
    eval(['load ' modelname]);
    D = cell(par.Patch_Channel,par.cls_num);
    W = cell(par.Patch_Channel,par.cls_num);
    for i = 1 : par.cls_num
        XN_t = double(Xn{i});
        XC_t = double(Xc{i});
        fprintf('SCSC: Channel: %d, Cluster: %d\n', j, i);
        % clean
        [GMM_Dc,GMM_Sc,~] = svd(model.covs(:,:,i));
        GMM_Sc = diag(GMM_Sc);
        % noisy
        [GMM_Dn,GMM_Sn,~] = svd(XN_t*XN_t'/size(XN_t,2));
        GMM_Sn = diag(GMM_Sn);
        % paired dictionary and weights
        SCSC.DC{j,i} = GMM_Dc;
        SCSC.DN{j,i} = GMM_Dn;
        SCSC.WC{j,i} = sqrt(GMM_Sc);
        SCSC.WN{j,i} = sqrt(GMM_Sn);
        % paired coefficient and mapping
        Alphac = GMM_Dc'*XC_t;
        Alphan = GMM_Dn'*XN_t;
        [Alphac, Alphan, Uc, Un, Pn, f] = SCSC_DAP_Learning(Alphac, Alphan, XC_t, XN_t, GMM_Dc, GMM_Dn, par);
        SCSC.UC{j,i} = Uc;
        SCSC.UN{j,i} = Un;
        SCSC.PN{j,i} = Pn;
        SCSC.f{j,i} = f;
%         Dict_BID = sprintf('Data/SCSC_PG_3_10_8x8_64_%s_20161003.mat',task);
        save Data/SCSC_PG_3_10_8x8_64_BID_20161003.mat SCSC '-v7.3'
    end
end