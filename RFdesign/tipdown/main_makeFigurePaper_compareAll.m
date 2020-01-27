% generate the figure that compared 8 methods for the joint design paper.
close all;
rFOV = 0;
for doOpt = 0:1
    for iniTraj = 1:4
        if rFOV == 1
            fname = ['forPaper_rFOV_' 'doOpt' num2str(doOpt) '_ktraj' num2str(iniTraj)];
            load(fname);
            pulseLength = 4e-3*length(gx);
            figure,im('notick',M_bloch,[0 max(abs(d_flip(:)))],'cbar',sprintf('nrmse: %.3f', nrmse));
            colormap(gca, 'default');
        else
            fname = ['forPaper_prephasing_' 'doOpt' num2str(doOpt) '_ktraj' num2str(iniTraj)];
            load(fname);
            pulseLength = 4e-3*length(gx);
            figure,im('notick',diff_bloch,[0 max(abs(d_flip(:)))],'cbar',sprintf('error image; nrmse: %.3f', nrmse));
            colormap(gca, 'default');
        end
        
        
        if doOpt == 0
            switch iniTraj
                case 1
                    % ylabel(sprintf('SoS \n(%.1f ms)',pulseLength));
                    ylabel('SoS');
                case 2
                    % ylabel(sprintf('SPINS \n(%.1f ms)',pulseLength));
                    ylabel('SPINS');
                case 3
                    % ylabel(sprintf('KT-points \n(%.1f ms)',pulseLength));
                    ylabel('KT-points');
                case 4
                    % ylabel(sprintf('extended KT-points \n(%.1f ms)',pulseLength));
                    ylabel('extended KT-points');
            end
        end
        
        if iniTraj == 4
            switch doOpt
                case 0
                    xlabel('initialization');
                case 1
                    xlabel('after interior-point optimization');
            end
        end
        fig=gcf;
        set(findall(fig,'-property','FontSize'),'FontSize',18)
        set(findall(fig,'-property','Font'),'Font','Helvetica')
        
        
        print('-depsc', fname);
        %copyfile([fname '.eps'], '~/Dropbox/MRI_research/mypapers/ktrajDesign/figs_paper/');
        copyfile([fname '.eps'], '~/Dropbox/MRI_research/mypapers/thesis/ktrajDesign/figs_paper/');

    end
end

%% plot desired pattern and fieldmap rFOV

if 1
load(['forPaper_rFOV_' 'doOpt' num2str(0) '_ktraj' num2str(1)]);
figure,im('notick',d_flip,[0 max(abs(d_flip(:)))],'cbar',sprintf('target excitation pattern'));
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','Font'),'Font','Helvetica')
print -depsc rfov_target;
copyfile rfov_target.eps ~/Dropbox/MRI_research/mypapers/ktrajDesign/figs_paper/

%colormap(gca, 'default');
figure,im('notick',b0map.*roi,'cbar',sprintf('B0 field map'));
colormap(gca, 'default');
set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(findall(gcf,'-property','Font'),'Font','Helvetica')
print -depsc rfov_fieldmap;
copyfile rfov_fieldmap.eps ~/Dropbox/MRI_research/mypapers/ktrajDesign/figs_paper/
end


%% plot ktraj
if 1
    for iniTraj = 1:4
        switch iniTraj
            case 1
                load k_sos;
            case 2
                load k_spins;
            case 3
                load k_ktpoints_rFOV;
            case 4
                load k_extendedKT_rFOV;
        end
        k3d = k3d(end-981:end,:);
        figure,h = plot3(k3d(:,1),k3d(:,2),k3d(:,3));
        xlabel('kx (cycle/cm)'); ylabel('ky (cycle/cm)'); zlabel('kz (cycle/cm)');  
        limits = 0.2; 
 %       xlim([-limits limits]); ylim([-limits limits]); zlim([-limits limits]);
        xtick([-limits, 0, limits]); ytick([-limits, 0, limits]); ztick([-limits, 0, limits]);
        fig=gcf;
        set(findall(fig,'-property','FontSize'),'FontSize',16)
        set(findall(fig,'-property','Font'),'Font','Helvetica')
        set(h, 'LineWidth', 2.0)
        kfname = ['ktraj2_' num2str(iniTraj)];
        
        print('-depsc', kfname);
        %copyfile([kfname '.eps'], '~/Dropbox/MRI_research/mypapers/ktrajDesign/figs_paper/');
  copyfile([kfname '.eps'], '~/Dropbox/MRI_research/mypapers/thesis/ktrajDesign/figs_paper/');
    end
end


