% Given A, an m-by-n data matrix, and L, a 2-by-n location matrix,
% use block coordinate descent to compute 
%
% min_{W_A,W~_A,H_A,W_L,H_L} ||A-W_A*H_A||_F^2 + alpha||H_A - H_A circ H_L||_F^2
% subject to W_A, H_A>=0
function [H_A, H_L, W_A, W_L, norm_diffs, nmf_ids, j_ids] = scHybridNMF(A, L, k, H_A, H_L, W_A, W_L, alpha)
    % Get optional arguments.
    num_iters = 300;
    tol = 1e-3;
    
    %initializes norm differences
    norm_diffs = zeros(num_iters+1,5);
    proj_norms = zeros(num_iters,1);
    
    %get m and n
    [m,n] = size(A);
    
    %initialize W_A (m-by-k), H_A (k-by-n) by nmf
    if isempty(W_A) && isempty(H_A)
        [W_A,H_A] = nmf(A,k,'type','sparse');
    end
    [~,nmf_ids] = max(H_A);

    %initialize W_L (2-by-k), H_L (k-by-n) by kmeans
    if isempty(W_L) && isempty(H_L)
        [j_ids,W_L] = kmeans(L.',k);
        H_L = double(bsxfun(@eq, j_ids(:), 1:max(j_ids)));
        H_L = H_L.';
        W_L = W_L.';
    else
        [~,j_ids] = max(H_L);
    end
    C = ones(size(H_L)) - H_L;
    
    %set parameters
    if isempty(alpha)
        alpha = norm(A,'fro')^2 / norm(H_A,'fro')^2;
    end
    
    %calculate initial norms
    norm_diffs(1,1) = norm(A - W_A * H_A, 'fro')^2;
    norm_diffs(1,2) = alpha * norm(H_A - H_A .* H_L, 'fro')^2;
    norm_diffs(1,3) = norm_diffs(1,1) + norm_diffs(1,2);
    norm_diffs(1,4) = norm(H_A - H_A .* H_L, 'fro')^2;
    norm_diffs(1,5) = norm_diffs(1,1) + norm_diffs(1,4);
    
    for i = 1 : num_iters
        %bcd for H_A
        left_H_A = W_A.' * W_A;
        right_H_A = W_A.' * A;
        for j = 1:k
            j_group = find(j_ids == j);
            C_j_group = ones(k,1);
            C_j_group(j) = 0;
            left_H_A_j = left_H_A + alpha * diag(C_j_group);
            H_A(:,j_group) = nnlsm_blockpivot(left_H_A_j, right_H_A(:,j_group), 1, H_A(:,j_group));
        end

        %bcd for W_A
        left_W_A = H_A * H_A.';
        right_W_A = H_A * A.';
        W_A = nnlsm_blockpivot(left_W_A, right_W_A, 1, W_A.');
        W_A = W_A.';   
        
        %normalize H_A and W_A
        norm2=sum(W_A,1);
        toNormalize = norm2>0;
        W_A(:,toNormalize) = W_A(:,toNormalize)./repmat(norm2(toNormalize),m,1);
        H_A(toNormalize,:) = H_A(toNormalize,:).*repmat(norm2(toNormalize)',1,n);
        
        %update norms 
        norm_diffs(i+1,1) = norm(A - W_A * H_A, 'fro')^2;
        norm_diffs(i+1,2) = alpha * norm(H_A - H_A .* H_L, 'fro')^2;
        norm_diffs(i+1,3) = norm_diffs(i+1,1) + norm_diffs(i+1,2);
        norm_diffs(i+1,4) = norm(H_A - H_A .* H_L, 'fro')^2;
        norm_diffs(i+1,5) = norm_diffs(i+1,1) + norm_diffs(i+1,4);
        
        grad = [2 * (right_H_A - left_H_A * H_A + alpha * C .* eye(size(C))), 2 * (right_W_A - left_W_A * W_A.')]; 
        if i == 1
            proj_norms(i) = norm(grad(grad<=0|[H_A,W_A.']>0));
            continue;
        else
            proj_norms(i) = norm(grad(grad<=0|[H_A,W_A.']>0));
        end
        
        if proj_norms(i) < tol * proj_norms(i)
            num_iters = i;
            break;
        end
        if proj_norms(i-1) - proj_norms(i) < tol && i > 5
            num_iters = i;
            break;
        end
    end
    
    %display objective function plot
    clf;
    plot(1:num_iters,norm_diffs(1:num_iters,1));
    hold on;
    plot(1:num_iters,norm_diffs(1:num_iters,2));
    plot(1:num_iters,norm_diffs(1:num_iters,3),'k-.');
    xlim([1 num_iters]);
    legend('$||A-W_AH_A||_F^2$','$\alpha||H_A-H_A \circ H_L||_F^2$', ...
       'Objective Function Value','Interpreter','latex','fontsize',8);
    fname = join([folder,'/',int2str(k),'_objectives.png']);
    saveas(gcf,fname);
    clf;
    
    %display norm difference plot
    plot(1:num_iters,norm_diffs(1:num_iters,1));
    hold on;
    plot(1:num_iters,norm_diffs(1:num_iters,4));
    plot(1:num_iters,norm_diffs(1:num_iters,5),'k-.');
    xlim([1 num_iters]);
    legend('$||A-W_AH_A||_F^2$','$||H_A-H_A \circ H_L||_F^2$', ...
           'Sum of Norm Differences','Interpreter','latex','fontsize',8);
    fname = join([folder,'/',int2str(k),'_norm_differences.png']);
    saveas(gcf,fname);
    clf;
end