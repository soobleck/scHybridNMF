% Given A, an m-by-n data matrix, and L, a 2-by-n location matrix,
% use block coordinate descent to compute:
% min_{W_A,H_A} ||A-W_A*H_A||_F^2 + alpha||H_A - H_A circ H_L||_F^2
% subject to W_A, H_A>=0
function [H_A, W_A, nmf_ids, j_ids] = scHybridNMF(A, L, k, alpha, tol, num_iters)
    % Get optional arguments.
    if isempty(num_iters)
        num_iters = 500;
    end
    if isempty(alpha)
        alpha = 0.2;
    end
    if isempty(tol)
        tol = 0.01;
    end

    %initializes norm differences
    proj_norms = zeros(num_iters,1);

    %get m and n
    [m,n] = size(A);

    %initialize W_A (m-by-k), H_A (k-by-n) by nmf
    beta = m / norm(A,'fro')^2;
    gamma = n / norm(A,'fro')^2;
    [~,H_A] = nmf(A,k,'type','sparse','alpha',beta,'beta',gamma);
    [~,nmf_ids] = max(H_A);

    %initialize W_L (2-by-k), H_L (k-by-n) by kmeans
    [j_ids,W_L] = kmeans(L.',k);
    H_L = double(bsxfun(@eq, j_ids(:), 1:max(j_ids)));
    H_L = H_L.';
    W_L = W_L.';
    [~,j_ids] = max(H_L);
    H_L = confidence_hl(L.',W_L,H_L);
    C = ones(size(H_L)) - H_L;
    
    grad_H_A = zeros(k,n);
    W_A = rand(m,k);
    for i = 1 : num_iters
        %bcd for W_A
        left_W_A = H_A * H_A.';
        right_W_A = H_A * A.';
        W_A = nnlsm_blockpivot(left_W_A, right_W_A, 1, W_A.');
        W_A = W_A.';
        
        %bcd for H_A
        left_H_A = W_A.' * W_A;
        right_H_A = W_A.' * A;
        for j = 1:n
            left_H_A_j = left_H_A + alpha * diag(C(:,j)).^2;
            [H_A(:,j),grad_H_A(:,j)] = nnlsm_blockpivot(left_H_A_j, right_H_A(:,j), 1, H_A(:,j));
        end
        
        left_W_A = H_A * H_A.';
        right_W_A = H_A * A.';
        grad = [grad_H_A, left_W_A * W_A.' - right_W_A];
        proj_norms(i) = norm(grad(grad<=0|[H_A,W_A.']>0));
        if proj_norms(i) < tol * proj_norms(1)
            break;
        end
    end
end

function [H_L] = confidence_hl(L,W_L,H_L)
    [~,n] = size(L);
    for i = 1:n
        dist = pdist2(L(:,i).', W_L.');
        [dist,idx] = mink(dist,2);
        comp_dist = dist(1) + dist(2);
        H_L(idx(1),i) = dist(2) / comp_dist;
        H_L(idx(2),i) = dist(1) / comp_dist;
    end
end
