function [] = vis_clusters(H, L, A_tsne, folder, alg_name, H_status) 
    %gene expression or A
    get_images(H, A_tsne, strcat(folder,'/gc.png'),H_status);
           
    %cell location or L
    get_images(H, L, strcat(folder,'/cc.png'),H_status);
           
    %H_A
    imwrite(cat(2,imread(strcat(folder,'/cc.png')), ...
      imread(strcat(folder,'/gc.png'))), ...
      strcat(folder,'/',alg_name,'.png'));
  
    delete(strcat(folder,'/gc.png'));
    delete(strcat(folder,'/cc.png'));
  
    %silhouettes
%     silhouette(L.',mems,'Euclidean');
%     xlim([-1 1]);
%     saveas(gcf, strcat(folder,'/cc_sil.png'));
%     clf;
%     silhouette(A_tsne.',mems);
%     xlim([-1 1]);
%     saveas(gcf, strcat(folder,'/gc_sil.png'));
end

function [mems] = get_images(H, S_tsne, fname, H_status)
    % if H is a matrix
    if H_status == 0
        [~,mems] = max(H);
    % if H is a vector
    else
        mems = H;
    end
        
    clf;
    cmap = hsv(max(mems));
    gscatter(S_tsne(:,1), S_tsne(:,2), mems, cmap, [], 12);
    set(legend,'visible','off');
%     set(gca,'XTick',[], 'YTick', []);
    xlim([min(S_tsne(:,1))-1 max(S_tsne(:,1))+1]);
    ylim([min(S_tsne(:,2))-1 max(S_tsne(:,2))+1]);
    xlabel("X");
    ylabel("Y");
    axis square;
    saveas(gcf,fname);
end