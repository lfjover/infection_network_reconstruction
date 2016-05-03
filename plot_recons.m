function plot_recons(Mtilde, Mrec, H, W,timeSeries)
    clf
    width = 19;
    height = 20;
    setfigure(width,height,10,5)

    if timeSeries
        ts = load(timeSeries);
        [t,x] = ts.allRuns{1,:};

        h = x(:,1:10);
        v = x(:,11:end);

        fs = 14;
        ax1 = axes('units','centimeters','pos',[ 2.5 15 6 4]);

        plot(t, v,'linewidth',1.5);
        title('Virus','fontsize',fs)
        xlim([t(1) t(end)])
        ylabel('Density (particles/ml)','fontsize',fs)
        xlabel('Time (h)','fontsize',fs)
        ax2 = axes('units','centimeters','pos',[ 11 15 6 4]);

        plot(t, h,'linewidth',1.5);
        title('Host','fontsize',fs)
        xlabel('Time (h)','fontsize',fs)
        xlim([t(1) t(end)])

    end

    [nr,nc] =size(W);

    ax3 = axes('units','centimeters','pos',[ 2.5 8 7 4]);
    imagesc(W)
    set(gca,'ytick',[],'fontsize',12)
    colormap(ax3,jet)
    colorbar

    ax4 = axes('units','centimeters','pos',[ 11 8 7 4]);
    imagesc(H);
    set(gca,'ytick',[],'fontsize',12)
    colormap(ax4,flipud(hot))
    colorbar

    axes('units','centimeters','position', [4 12.6 1 1],'visible','off')
    text(0,0,'$W$','fontsize',40,'interpreter','latex')   
    axes('units','centimeters','position', [12.5 12.6 1 1],'visible','off')
    text(0,0,'$H$','fontsize',40,'interpreter','latex')


    % Mtilde and Mrec good single experiment

    axes('units','centimeters','pos',[ 2.5 .5 5 5])
    imagesc(Mtilde)
    set(gca,'ytick',[],'xtick',[])
    colormap jet
    map = colormap;
    map(1,:) = 1;
    colormap(map)
    s = caxis;


    axes('units','centimeters','pos',[ 11 .5 5 5]) 
    imagesc(Mrec)
    caxis(s);
    set(gca,'ytick',[],'xtick',[])
    
    axes('units','centimeters','position', [4.2 6.2 1 1],'visible','off')
    text(0,0,'$\tilde{M}$','fontsize',35,'interpreter','latex')
    axes('units','centimeters','position', [12.5 6.2  1 1],'visible','off')
    text(0,0,'$\tilde{M}_{rec}$',...
        'fontsize',35,'interpreter','latex')
    
    %a) b) c)
    axes('units','centimeters','position', [0.4 19.5 1 1],'visible','off')
    text(0,0,'a','fontsize',35,'interpreter','latex')
    axes('units','centimeters','position', [0.4 12.5 1 1],'visible','off')
    text(0,0,'b','fontsize',35,'interpreter','latex')
    axes('units','centimeters','position', [0.4  6 1 1],'visible','off')
    text(0,0,'c','fontsize',35,'interpreter','latex')
    
    
    
