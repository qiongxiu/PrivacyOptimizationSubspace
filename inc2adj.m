function mAdj = inc2adj(mInc)
if ~issparse(mInc)
    mInc = sparse(mInc);  
end
if ~all(ismember(mInc(:), [-1, 0, 1]))
    error('inc2adj:wrongMatrixInput', 'Matrix must contain only {-1,0,1}');    
end
if ~all(any(mInc, 2))
    error('inc2adj:wrongMatrixInput', 'Invalid incidence matrix - each edge must be connected to at least one node.');
end
mInc = mInc.';
if any(mInc(:) == -1)   % directed graph
    iN_nodes = size(mInc,1);       % columns must be vertices!!!
    
    [vNodes1, dummy] = find(mInc == 1);    % since MATLAB 2009b 'dummy' can be replaced by '~'
    [vNodes2, dummy] = find(mInc == -1);   % since MATLAB 2009b 'dummy' can be replaced by '~'
    
    mAdj = sparse(vNodes1, vNodes2, 1, iN_nodes, iN_nodes);
            
else    % undirected graph
    
    L    = mInc*mInc.';        % using Laplacian
    mAdj = L - diag(diag(L));
    
end
if any(mAdj(:) > 1)
    warning('inc2adj:wrongMatrixInput', 'Multi-edge detected!');
end
end