% Given A, an m-by-n data matrix, and L, a 2-by-n location matrix,
% use block coordinate descent to compute
%
% min_{W_A,H_A} ||A-W_A*H_A||_F^2 + alpha*||H_A - H_A circ H_L||_F^2 
% + beta*||W_A||_F^2 + gamma*sum_{i=1}^n ||H_A(:,i)||_F^2
% subject to W_A, H_A>=0
function [H_A, W_A, nmf_labels, km_labels] = scHybridNMF(A, L, k, H_A, H_L, W_A, W_L, alpha, beta, gamma)
    % Get optional arguments.
    num_iters = 50;
    tol = 0.5;

    %initializes norm differences
    proj_norms = zeros(num_iters,1);

    %get m and n
    [m,n] = size(A);

    %initialize W_A (m-by-k), H_A (k-by-n) by nmf
    if isempty(W_A) && isempty(H_A)
        [W_A,H_A] = nmf(A,k,'type','sparse','alpha',beta,'beta',gamma);
        [~,nmf_labels] = max(H_A);
        while numel(unique(nmf_labels)) < k
            [W_A,H_A] = nmf(A,k,'type','sparse','alpha',beta,'beta',gamma);
            [~,nmf_labels] = max(H_A);
        end
    end
    [~,nmf_labels] = max(H_A);

    %initialize W_L (2-by-k), H_L (k-by-n) by kmeans
    if isempty(H_L)
        start_centroids = zeros(k,2);
        for k_i = 1:k
            k_group = find(nmf_labels == k_i);
            start_centroids(k_i,:) = mean(L(k_group,:));
        end
        [km_labels,W_L] = kmeans(L,k,'Start',start_centroids);
        H_L = double(bsxfun(@eq, km_labels(:), 1:max(km_labels)));
        H_L = H_L.';
        W_L = W_L.';
    else
        [~,km_labels] = max(H_L);
    end
    H_L = confidence_hl(L.',W_L,H_L);
    C = ones(size(H_L)) - H_L;

    %set parameters
    if isempty(alpha)
        alpha = norm(A,'fro')^2 / norm(H_A,'fro')^2;
    end
    if isempty(beta)
        beta = mean(A(:)) * m / norm(A,'fro')^2;
    end
    if isempty(gamma)
        gamma = n / norm(A,'fro')^2;
    end

    grad_H_A = zeros(k,n);
    for i = 1 : num_iters
        %bcd for H_A
        left_H_A = W_A.' * W_A + gamma * ones(k);
        right_H_A = W_A.' * A;
        for j = 1:n
            left_H_A_j = left_H_A + alpha * diag(C(:,j)).^2;
            [H_A(:,j),grad_H_A(:,j)] = nnlsm_blockpivot(left_H_A_j, right_H_A(:,j), 1, H_A(:,j));
        end
        
        %bcd for W_A
        left_W_A = H_A * H_A.' + beta * eye(k);
        right_W_A = H_A * A.';
        W_A = nnlsm_blockpivot(left_W_A, right_W_A, 1, W_A.');
        W_A = W_A.';
        
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
