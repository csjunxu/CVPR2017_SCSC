clear;
addpath('Data');
addpath('Utilities');
addpath('SPAMS');

GT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_MeanImage\';
GT_fpath = fullfile(GT_Original_image_dir, '*.png');
TT_Original_image_dir = 'C:\Users\csjunxu\Desktop\CVPR2017\cc_Results\Real_NoisyImage\';
TT_fpath = fullfile(TT_Original_image_dir, '*.png');
GT_im_dir  = dir(GT_fpath);
TT_im_dir  = dir(TT_fpath);
im_num = length(TT_im_dir);

params = 'Data/params_PG.mat';
load(params,'par','param');
Dict_SR_backup = 'Data/SCSC_PG_3_10_8x8_64_BID_20161003.mat';
load(Dict_SR_backup,'SCSC');

for lambda2 = [0.0001 0.001 0.01 0.1]
    param.lambda2 = lambda2;
    for lambda = [0.001 0.01 0.1]
        param.lambda = lambda;
        for sqrtmu = [0.01 0.1 0.5]
            par.sqrtmu = sqrtmu;
            PSNR = [];
            SSIM = [];
            CCPSNR = [];
            CCSSIM = [];
            for i = 1:im_num
                IMin = im2double(imread(fullfile(TT_Original_image_dir,TT_im_dir(i).name) ));
                IM_GT = im2double(imread(fullfile(GT_Original_image_dir, GT_im_dir(i).name)));
                S = regexp(TT_im_dir(i).name, '\.', 'split');
                IMname = S{1};
                [h,w,ch] = size(IMin);
                fprintf('%s: \n',TT_im_dir(i).name);
                CCPSNR = [CCPSNR csnr( IMin*255,IM_GT*255, 0, 0 )];
                CCSSIM = [CCSSIM cal_ssim( IMin*255, IM_GT*255, 0, 0 )];
                fprintf('The initial PSNR = %2.4f, SSIM = %2.4f. \n', CCPSNR(end), CCSSIM(end));
                %% 3 channels
                IMout = zeros(size(IMin));
                for cc = 1:ch
                    modelname = sprintf('../DSCDL_BID/Data/GMM_PG_%d_10_8x8_64_20161003T094301.mat',cc);
                    eval(['load ' modelname]);
                    par.cc= cc;
                    par.nInnerLoop = 1;
                    IMin_cc = IMin(:,:,cc);
                    IM_GT_cc = IM_GT(:,:,cc);
                    fprintf('Channel %d: The initial PSNR = %2.4f, SSIM = %2.4f. \n', cc, csnr( IMin_cc*255,IM_GT_cc*255, 0, 0 ), cal_ssim( IMin_cc*255, IM_GT_cc*255, 0, 0 ));
                    IMout_cc = SCSC_PG_3Chs_BID(IMin_cc,IM_GT_cc,model,SCSC,par,param);
                    fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', csnr( IMout_cc*255, IM_GT_cc*255, 0, 0 ),cal_ssim( IMout_cc*255, IM_GT_cc*255, 0, 0 ));
                    IMout(:,:,cc) = IMout_cc;
                end
                %% output
                PSNR = [PSNR csnr( IMout*255, IM_GT*255, 0, 0 )];
                SSIM = [SSIM cal_ssim( IMout*255, IM_GT*255, 0, 0 )];
                fprintf('The final PSNR = %2.4f, SSIM = %2.4f. \n', PSNR(end), SSIM(end));
                %% output
                %                 imwrite(IMout, ['../cc_Results/Real_SCSC/SCSC_PG_3Chs_BID_' IMname '.png']);
            end
            mPSNR = mean(PSNR);
            mSSIM = mean(SSIM);
            mCCPSNR = mean(CCPSNR);
            mCCSSIM = mean(CCSSIM);
            save(['Real_SCSC_PG_3Chs_BID_' num2str(lambda2) '_' num2str(lambda) '_'  num2str(sqrtmu) '_'  '.mat'],'PSNR','mPSNR','SSIM','mSSIM','CCPSNR','mCCPSNR','CCSSIM','mCCSSIM');
        end
    end
end