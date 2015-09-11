function Atilde = OMNICONT_construct_Atilde(C,NODES)
%--------------------------------------------------------------------------
% Constructing sparse matrix Atilde (Section 3.7)
% Atilde is built in blocks which are held in Ccell and Icell
% Ccell contains lambda+1 blocks (lambda defined in thesis)
% Blocks l=1,...,lambda look like
% [ D^(l)      U^(l);
% [ (V^(l))^*   0   ],
% while the final block holds S.
%
% Icell has lambda blocks which look like
% [0  -I;
%  -I  0]
%
% Inputs:  NODES  : Binary tree containing all skeletons, and block indices
% Outputs: Atilde : Sparse system matrix
%--------------------------------------------------------------------------

% BLOCKS FROM FIRST LEVEL
%--------------------------
% This first section constructs D^(1), U^(1),(V^(1))^*=(U^(1))^*,
% corresponding to level 1:
% [ D^(1)    U^(1);
%  (V^(1))^* 0    ]

% D^(1) is built from the self-interaction of diagonal blocks

N   = size(C,2);
nboxes = size(NODES,2);
num_levels=log2(nboxes+1);
Dcell=cell(1,2^(num_levels-1));  % holds diagonal blocks, D^(1)
Ucell=cell(1,2^(num_levels-1));  % holds U^(1)
Ccell=cell(1,num_levels);        % holds the blocks of Atilde

dcount=2^(num_levels-1);         % counter for number of leaves on finest level
sumk=0;                          % sum of ranks- determines dimension of U

for ibox = nboxes:(-1):2   % counting blocks from the bottom of the tree up
    if (NODES{5,ibox}==0)  % if block in tree (ibox) has no children (finest level)
        ind = NODES{6,ibox} - 1 + (1:NODES{7,ibox}); % indices of points in block
        % building D:
        Dcell{1,dcount} = sparse(OMNICONT_construct_A_diag(C,ind));
        % building U:
        k=NODES{9,ibox};
        sumk=sumk+k;
        J=NODES{13, ibox};
        Ui=zeros(length(ind),k);
        Ui(J(1:k),1:k)=speye(k);
        Ui(J(k+1:end),1:k)=NODES{12,ibox}';
        Ucell{1,dcount}=sparse(Ui);
        dcount=dcount-1;
    end
end
U=sparse(blkdiag(Ucell{:}));
Ccell{1,1}=sparse([sparse(blkdiag(Dcell{:})) U; U' zeros(sumk)]);

% figure; spy(Ccell{1,1}) (for debugging)

% BLOCKS FROM REMAINING LEVELS
%-------------------------------
% This section constructs the remaining blocks D^(l), U^(l), (U^(l))^*,
% corresponding to the remaining levels
% The blocks D^(l) are built from sibling interactions from the previous
% (finer) level
% For each level, D^(l) is built first from children on the previous level,
% followed by U^(l) and (U^(l))^*
% Then the identity matrices are built: [0  -I;
%                                        -I  0]
% Then the block [D^(l)     U^(l); is built
%                 (V^(l))^*  0   ]
%

level=NODES{2,nboxes};      % finest level (counts opposite to thesis)
Icell=cell(1,level);        % holds identity matrices for each level in Atilde
Dcell=cell(1,2^(level-1));  % holds diagonal blocks, D^(l) from each level
Ucell=cell(1,2^(level-1));  % holds blocks U^(l) on each level
count=2^(level-1);          % number of subblocks in D^(l) and U^(l)
ccount=2;                   % counting index for Ccell
icount=1;                   % counting index for Icell
sumk=0;
sumkD=0;
for ibox = nboxes:(-1):1    % counting from the bottom of the tree up
    if (NODES{5,ibox}==2)   % if block in tree has 2 children
        % obtaining ibox (block index) of children:
        ison1 = NODES{4,ibox}(1);
        ison2 = NODES{4,ibox}(2);
        
        % skipping building identity matrices and Ccell until D^(l), U^(l)
        % are built (all children are processed on current level)
        if NODES{2,ison1}~=level  % if child is not on the level we are
                                  % processing then we've reached a new
                                  % level and processed all children.
                                  % We can then build D^(l), U^(l), I for
                                  % the current level:
            
            % building Icell
            I=speye(sumkD);
            Icell{1,icount}=sparse([zeros(sumkD) -I; -I zeros(sumkD)]);
            icount=icount+1;
            
            %building U^(l)
            U=sparse(blkdiag(Ucell{:}));
            
            % building block of Atilde
            Ccell{1,ccount}=sparse([sparse(blkdiag(Dcell{:})) U; U' zeros(sumk)]);
            
            ccount=ccount+1;
            level=level-1;   % INCREMENTING LEVEL (opposite to thesis)
            sumk=0;
            sumkD=0;
            % reset Dcell and Ucell
            Dcell=cell(1,2^(level-1));
            Ucell=cell(1,2^(level-1));
            count=2^(level-1); % keeps track of number of subblocks in D, U
        end
        
        % building D^(l), U^(l):
        k1=NODES{9,ison1};
        k2=NODES{9,ison2};
        sum=k1+k2;
        sumkD=sumkD+sum;
        % building D^(l)
        Dcell{1,count}=sparse([zeros(k1) OMNICONT_construct_A_offd(C,NODES{11,ison1},NODES{11,ison2});...
            OMNICONT_construct_A_offd(C,NODES{11,ison2},NODES{11,ison1}) zeros(k2)]);
        
        if ibox~=1  % if we haven't reached the top of the tree
            k=NODES{9,ibox};
            sumk=sumk+k;
            Ik=speye(k);
            J=NODES{13,ibox};
            Ui=zeros(sum, k);
            Ui(J(1:k),1:k)=Ik;
            Ui(J(k+1:end),1:k)=NODES{12,ibox}';
            Ucell{1,count}=sparse(Ui);
            
        else % if we have reached the top of the tree, then build [0 -I; -I S]
            S=sparse(blkdiag(Dcell{:}));
            Ccell{1,ccount}=S;
            I=speye(size(S,1));
            Icell{1,icount}=sparse([zeros(size(I)) -I; -I zeros(size(I))]);
        end
        count=count-1;
    end
end

Atilde=sparse(blkdiag(Ccell{:}));   % blocks in Atilde holding D^(l), U^(l)
Atilde_I=sparse(blkdiag(Icell{:})); % blocks in Atilde holding I
Atilde(N+1:end,N+1:end)=Atilde(N+1:end,N+1:end)+Atilde_I;

%figure; spy(Atilde) (for debugging)

return