% this section must be modified with the values adapted for the user's
% study case

C = 20; % chromosomes number
LowRes = 10000000; % enter the chosen low resolution (to compute whole-genome reference structure)
TADsRes = 2000000; % enter the chosen TADs-resolution (to identify TAD boundaries)
MedRes = 1000000; % enter the chosen medium resolution (to compute chromosomes structure)
HiRes = 100000; % enter the chosen high resolution (to compute TADs structure)

ResRatioT = TADsRes/HiRes; % scaling factor between the resolution used to identify TAD boundaries and the high resolution (used to compute TADs structure)
ResRatioM = MedRes/HiRes; % scaling factor between the medium resolution (used to compute chromosomes structure) and the high resolution (used to compute TADs structure)

%% 
% note: chromosome X is called '20'

% load:
% M_LowRes = whole-genome low resolution contact map (inter- and intra-chromosomal)
% XYZ_LowRes_gen = low resolution 3D structure of the whole genome (computed with the custom 3D reconstruction algorithm); N*3 matrix (x,y,z coordinates of the N points of the structure)
% L_LowRes = array of C components (C = number of chromosomes); each component is the number of bins of each chromosome in the low resolution contact map
% M_TADsRes_#c = intra-chromosomal contact map of each chromosome at the resolution used to compute TAD boundaries (TADsRes); the #c symbol must be replaced with the chromosome number; a total of C matrixes (number of chromosomes)
% TB_TADsRes_#c = array of TAD boundaries position (array of the bins of the TADs-resolution matrix in which there is a TAD boundary); the #c symbol must be replaced with the chromosome number; the length of the array is equal to the number of TAD boundaries identified in the chromosome number #c; a total of C arrays (number of chromosomes)
% M_MedRes_#c = intra-chromosomal contact map of each chromosome at medium resolution (MedRes); the #c symbol must be replaced with the chromosome number; a total of C matrixes (number of chromosomes)
% XYZ_MedRes_chr#c = medium resolution 3D structure of each chromosome (computed with the custom 3D reconstruction algorithm from medium resolution intra-chromosomal contact maps resulted from centromere region removal, here called 'S.MedRes.mapp_regs.chr_#c'); N*3 matrix (x,y,z coordinates of the N points of the structure); the #c symbol must be replaced with the chromosome number; a total of C matrixes (number of chromosomes)
% M_HiRes_#c = intra-chromosomal contact map of each chromosome at high resolution (HiRes); the #c symbol must be replaced with the chromosome number; a total of C matrixes (number of chromosomes)
% XYZ_HiRes_chr#c_TAD#t = high resolution 3D structure of each TAD (computed with the custom 3D reconstruction algorithm from high resolution TADs matrix blocks resulted from centromere region removal, here called 'S.HiRes.mapp_regs.chr_#c.block_#t); N*3 matrix (x,y,z coordinates of the N points of the structure); the #c symbol must be replaced with the chromosome number; the #t symbol must be replaced with the TAD number; a total number of matrixes equal to C * Tc (Tc = number of TADs of each chromosome)

S = [];
S.LowRes.Map = M_LowRes;
S.LowRes.XYZ_gen = XYZ_LowRes_gen; 

b = 0;
for i = 1:C    % for each chromosome ('i' index will refer to chromosome number)
    
    % LowRes
    eval(['S.LowRes.length_chrs.length_chr_' num2str(i) '= L_LowRes(i); ' ]); % store in S the lengths of each low resolution chromosome
    a = 1 + b;
    b = a + eval(['S.LowRes.length_chrs.length_chr_' num2str(i) ]) - 1;
    eval(['S.LowRes.XYZ_chrs.XYZ_' num2str(i) '= S.LowRes.XYZ_gen( a:b , :);']); % separate the whole-genome low resolution 3D structure in the low resolution 3D structures corresponding to each chromosome, and store them in S
    
    % TADsRes
    eval(['S.TADsRes.Maps_chrs.M_' num2str(i) '= M_TADsRes_' num2str(i) ';']); 
    eval(['S.TADsRes.TB_chrs.TB_' num2str(i) '= TB_TADsRes_' num2str(i) ';']); 
    %idx_TADsRes = sum(eval(['S.TADsRes.Maps_chrs.M_' num2str(i) ';'])) == 0; % logical array with 1 in bins corresponding to centromere region (zero bins in the TADs-resolution intra-chromosomal contact map)
    if eval(['S.TADsRes.TB_chrs.TB_' num2str(i) '(1) == 1']) || eval(['S.TADsRes.TB_chrs.TB_' num2str(i) '(1) == 2']) % removal of the first TAD boundary if it is positioned in the first or second bin
            eval(['S.TADsRes.TB_chrs.TB_' num2str(i) '(1) = [];']);
    end
    if eval(['S.TADsRes.TB_chrs.TB_' num2str(i) '(end) == length(S.TADsRes.Maps_chrs.M_' num2str(i) ')']) || eval(['S.TADsRes.TB_chrs.TB_' num2str(i) '(end) == length(S.TADsRes.Maps_chrs.M_' num2str(i) ')-1']) % removal of the last TAD boundary if it is positioned in the last or penultimate bin 
        eval(['S.TADsRes.TB_chrs.TB_' num2str(i) '(end) = [];']);
    end
    %map_TADsRes = find(idx_TADsRes == 0); % array wich maps the original bin numbers (comprising the centromere; corresponding to the array coefficients) to the new bin numbers (not comprising the centromere; corresponding to the array coefficient indexes)
    %eval(['S.TADsRes.TB_chrs.TB_' num2str(i) '= map_TADsRes(S.TADsRes.TB_chrs.TB_' num2str(i) ');']); % maps the TAD boundaries positions computed on contact maps comprising the centromere on the maps which don't comprise the centromere
    
    % medium resolution
    eval(['S.MedRes.Maps_chrs.chr_' num2str(i) '= M_MedRes_' num2str(i) ';']);
    eval(['S.MedRes.length_chrs.chr_' num2str(i) '= length(S.MedRes.Maps_chrs.chr_' num2str(i) ') ;' ]); % store in S the lengths of each medium resolution chromosome
    eval(['S.MedRes.idx_chrs.chr_' num2str(i) '= sum(S.MedRes.Maps_chrs.chr_' num2str(i) ')== 0;']); % logical array with 1 in bins corresponding to centromere region (zero bins in the medium resolution intra-chromosomal contact map)
    
    eval(['S.MedRes.mapp_regs.chr_' num2str(i) '=  S.MedRes.Maps_chrs.chr_' num2str(i) ';' ]); % removal of zero bins (corresponding to centromere) from each medium resolution intra-chromosomal contact map
    eval(['S.MedRes.mapp_regs.chr_' num2str(i) '(S.MedRes.idx_chrs.chr_' num2str(i) ', :) = []; ' ]);
    eval(['S.MedRes.mapp_regs.chr_' num2str(i) '(:, S.MedRes.idx_chrs.chr_' num2str(i) ') = []; ' ]);
    
    eval(['l = length(S.MedRes.mapp_regs.chr_' num2str(i) ');' ]); % chromosomal length not considering centromere region
    
    eval(['S.MedRes.XYZ_mapp_regs.chr_' num2str(i) '= XYZ_MedRes_chr' num2str(i) ';']); 
    
    eval(['S.MedRes.XYZ_chrs.XYZ_' num2str(i) '= zeros(S.MedRes.length_chrs.chr_' num2str(i) ', 3);']); % reindexing of each chromosome 3D structure in order to take into account centromere region (corresponding to zero bins of the contact map)
    eval(['map_MedRes = find(S.MedRes.idx_chrs.chr_' num2str(i) '== 0) ;' ]); % array wich maps the original bin numbers (comprising the centromere; corresponding to the array coefficients) to the new bin numbers (not comprising the centromere; corresponding to the array coefficient indexes)
    for k = 1:length(map_MedRes)
        eval(['S.MedRes.XYZ_chrs.XYZ_' num2str(i) '(map_MedRes(k), :) = S.MedRes.XYZ_mapp_regs.chr_' num2str(i) '(k, :) ;' ]); % chromosome 3D structure with zero coordinates given to the centromere region
    end
    
    % high resolution
    eval(['S.HiRes.Maps_chrs.M_' num2str(i) '= M_HiRes_' num2str(i) ';']); 
    eval(['S.HiRes.TB_chrs.TB_' num2str(i) '= S.TADsRes.TB_chrs.TB_' num2str(i) '* ResRatioT - ResRatioT/2;']); % mapping of TAD boundaries from their position in TADsRes matrixes to their position in HiRes matrixes
    l_HiRes = length( eval(['S.HiRes.Maps_chrs.M_' num2str(i) ';']) ); % high resolution chromosome length 
    t = length(eval(['S.HiRes.TB_chrs.TB_' num2str(i) ';'])); % number of TAD boundaries in the chromosome
    start = 0;
    for j = 1:t+1  % for each TAD ('j' index will refer to TAD number, inside each chromosome)
        
        if j ~= t+1     % for all TADs of the considered chromosome with the exception of the last one (the last TAD needs specific code, see as follows)
            
            % splitting of high resolution intra-chromosomal contact maps in blocks corresponding to TADs
            final = eval(['S.HiRes.TB_chrs.TB_' num2str(i) '(' num2str(j) ');' ]); 
            eval(['S.HiRes.blocks_chrs.chr_' num2str(i) '.block_' num2str(j) '= S.HiRes.Maps_chrs.M_' num2str(i) '(start+1 : final, start+1 : final);']); 
            
            eval(['S.HiRes.TADs_dim.chr_' num2str(i) '.block_' num2str(j) '= final - start ;']); % length of each TAD (number of bins of each block)
            eval(['S.HiRes.TADs_start.chr_' num2str(i) '.block_' num2str(j) '= start+1 ;' ]); % start of each TAD (first bin of each block in the high resolution intra-chromosomal contact map)
            eval(['S.HiRes.TADs_end.chr_' num2str(i) '.block_' num2str(j) '= final ;' ]); % end of each TAD (last bin of each block in the high resolution intra-chromosomal contact map)
            start = final;
            
            % removal of bins corresponding to the centromere from each TAD high resolution contact map 
            eval(['idx = sum(S.HiRes.blocks_chrs.chr_' num2str(i) '.block_' num2str(j) ')== 0;' ]); % logical array with 1 in bins corresponding to zero bins in the high resolution blocks (corresponding to TADs)
            idx_0 = centromere(idx, 20); % logical array with 1 only in bins corresponding to the centromere (considered as at least 20 consecutive zero bins)
            eval(['S.HiRes.idx_blocks_chrs.chr_' num2str(i) '.idx_block_' num2str(j) '= idx_0;' ]); % stores idx_0 for each high resolution block (bins of the centromere)
            eval(['S.HiRes.mapp_regs.chr_' num2str(i) '.block_' num2str(j) '=  S.HiRes.blocks_chrs.chr_' num2str(i) '.block_' num2str(j) ';' ]); 
            eval(['S.HiRes.mapp_regs.chr_' num2str(i) '.block_' num2str(j) '(idx_0, :) = []; ' ]);
            eval(['S.HiRes.mapp_regs.chr_' num2str(i) '.block_' num2str(j) '(:, idx_0) = []; ' ]);
            
            if eval(['isempty(S.HiRes.mapp_regs.chr_' num2str(i) '.block_' num2str(j) ') == 1'])
                eval(['S.HiRes.XYZ_mapp_regs.chr_' num2str(i) '.block_' num2str(j) '= [];' ]);
            else
                eval(['S.HiRes.XYZ_mapp_regs.chr_' num2str(i) '.block_' num2str(j) '= XYZ_HiRes_chr' num2str(i) '_TAD' num2str(j) ';']); % stores in S the high resolution 3D structure of each TAD ('XYZ_HiRes_chr#c_TAD#t')
            end
            
            % reindexing of each TAD 3D structure in order to comprise centromere regions (corresponding to zero bins)
            eval(['S.HiRes.XYZ_TADs_chrs.chr_' num2str(i) '.block_' num2str(j) '= zeros(S.HiRes.TADs_dim.chr_' num2str(i) '.block_' num2str(j) ', 3);']); 
            map_HiRes = find(idx_0 == 0); % array wich maps the original bin numbers (comprising the centromere; corresponding to the array coefficients) to the new bin numbers (not comprising the centromere; corresponding to the array coefficient indexes)
            for k = 1:length(map_HiRes)
                eval(['S.HiRes.XYZ_TADs_chrs.chr_' num2str(i) '.block_' num2str(j) '(map_HiRes(k), :) = S.HiRes.XYZ_mapp_regs.chr_' num2str(i) '.block_' num2str(j) '(k, :) ;' ]); % TAD 3D structure with zero coordinates given to the centromere region
            end
        
            % mapping of TADs_dim, TADs_start, TADs_end from HiRes to MedRes
            eval(['S.MedRes.TADs_dim.chr_' num2str(i) '.block_' num2str(j) '= fix(S.HiRes.TADs_dim.chr_' num2str(i) '.block_' num2str(j) '/ResRatioM);']); % stores in S the length of each TAD at medium resolution 
            eval(['S.MedRes.TADs_start.chr_' num2str(i) '.block_' num2str(j) '= 1 + fix(S.HiRes.TADs_start.chr_' num2str(i) '.block_' num2str(j) '/ResRatioM);' ]); % stores in S the start position of each TAD at medium resolution (first bin of the TAD block in the medium resolution intra-chromosomal contact map)
            eval(['S.MedRes.TADs_end.chr_' num2str(i) '.block_' num2str(j) '= fix(S.HiRes.TADs_end.chr_' num2str(i) '.block_' num2str(j) '/ResRatioM);' ]); % stores in S the end position of each TAD at medium resolution (last bin of the TAD block in the medium resolution intra-chromosomal contact map)
            eval(['S.MedRes.XYZ_TADs.chr_' num2str(i) '.block_' num2str(j) '= S.MedRes.XYZ_chrs.XYZ_' num2str(i) '( S.MedRes.TADs_start.chr_' num2str(i) '.block_' num2str(j) ': S.MedRes.TADs_end.chr_' num2str(i) '.block_' num2str(j) ', :);' ]); % splitting of each medium resolution chromosome 3D structure in medium resolution 3D structures corresponding to TADs identified in the chromosome
        
            % removal of centromere coordinates from TADs medium resolution 3D structures
            XYZ = eval(['S.MedRes.XYZ_TADs.chr_' num2str(i) '.block_' num2str(j) ';' ]); 
            z = zeros(1,length(XYZ));
            for h = 1:length(XYZ)
                if sum(XYZ(h,:)) == 0
                    z(h) = 1;
                end
            end
            z=logical(z);
            XYZ(z,:)=[]; 
            eval(['S.MedRes.XYZ_TADs_mapp.chr_' num2str(i) '.block_' num2str(j) '= XYZ;' ]); % stores in S the medium resolution 3D structures corresponding to TADs resulted from centromere removal
            
    
        else    % repetition of all operations on the last TAD block 
            
            eval(['S.HiRes.blocks_chrs.chr_' num2str(i) '.block_' num2str(t+1) '= S.HiRes.Maps_chrs.M_' num2str(i) '(start+1 : l_HiRes, start+1 : l_HiRes);']); 
            
            eval(['S.HiRes.TADs_dim.chr_' num2str(i) '.block_' num2str(t+1) '= l_HiRes - start ;']);
            eval(['S.HiRes.TADs_start.chr_' num2str(i) '.block_' num2str(t+1) '= start+1 ;' ]);
            eval(['S.HiRes.TADs_end.chr_' num2str(i) '.block_' num2str(t+1) '= l_HiRes ;' ]);
            
            eval(['idx = sum(S.HiRes.blocks_chrs.chr_' num2str(i) '.block_' num2str(t+1) ')== 0;' ]);
            idx_0 = centromere(idx, 20);
            eval(['S.HiRes.idx_blocks_chrs.chr_' num2str(i) '.idx_block_' num2str(t+1) '= idx_0;' ]);
            eval(['S.HiRes.mapp_regs.chr_' num2str(i) '.block_' num2str(t+1) '=  S.HiRes.blocks_chrs.chr_' num2str(i) '.block_' num2str(t+1) ';' ]);
            eval(['S.HiRes.mapp_regs.chr_' num2str(i) '.block_' num2str(t+1) '(S.HiRes.idx_blocks_chrs.chr_' num2str(i) '.idx_block_'  num2str(t+1) ', :) = []; ' ]);
            eval(['S.HiRes.mapp_regs.chr_' num2str(i) '.block_' num2str(t+1) '(:, S.HiRes.idx_blocks_chrs.chr_' num2str(i) '.idx_block_'  num2str(t+1) ') = []; ' ]);
            
            if eval(['isempty(S.HiRes.mapp_regs.chr_' num2str(i) '.block_' num2str(t+1) ') == 1'])
                eval(['S.HiRes.XYZ_mapp_regs.chr_' num2str(i) '.block_' num2str(t+1) '= [];' ]);
            else
                eval(['S.HiRes.XYZ_mapp_regs.chr_' num2str(i) '.block_' num2str(t+1) '= XYZ_HiRes_chr' num2str(i) '_TAD' num2str(t+1) ';']); % stores in S the high resolution 3D structure of each TAD ('XYZ_HiRes_chr#c_TAD#t')
            end
            
            eval(['S.HiRes.XYZ_mapp_regs.chr_' num2str(i) '.block_' num2str(t+1) '= XYZ_HiRes_chr' num2str(i) '_TAD' num2str(t+1) ';']);
            
            eval(['S.HiRes.XYZ_TADs_chrs.chr_' num2str(i) '.block_' num2str(t+1) '= zeros(S.HiRes.TADs_dim.chr_' num2str(i) '.block_' num2str(t+1) ', 3);']);
            map_HiRes = find(idx_0 == 0);
            for k = 1:length(map_HiRes)
                eval(['S.HiRes.XYZ_TADs_chrs.chr_' num2str(i) '.block_' num2str(t+1) '(map_HiRes(k), :) = S.HiRes.XYZ_mapp_regs.chr_' num2str(i) '.block_' num2str(t+1) '(k, :) ;' ]);
            end
            
            eval(['S.MedRes.TADs_dim.chr_' num2str(i) '.block_' num2str(t+1) '= fix(S.HiRes.TADs_dim.chr_' num2str(i) '.block_' num2str(t+1) '/ResRatioM);']);
            eval(['S.MedRes.TADs_start.chr_' num2str(i) '.block_' num2str(t+1) '= 1 + fix(S.HiRes.TADs_start.chr_' num2str(i) '.block_' num2str(t+1) '/ResRatioM);' ]);
            eval(['S.MedRes.TADs_end.chr_' num2str(i) '.block_' num2str(t+1) '= fix(S.HiRes.TADs_end.chr_' num2str(i) '.block_' num2str(t+1) '/ResRatioM);' ]);
            eval(['S.MedRes.XYZ_TADs.chr_' num2str(i) '.block_' num2str(t+1) '= S.MedRes.XYZ_chrs.XYZ_' num2str(i) '(S.MedRes.TADs_start.chr_' num2str(i) '.block_' num2str(t+1) ': length(S.MedRes.XYZ_chrs.XYZ_' num2str(i) '), :);' ]);
            
            XYZ = eval(['S.MedRes.XYZ_TADs.chr_' num2str(i) '.block_' num2str(t+1) ';' ]);
            z = zeros(1,length(XYZ));
            for h = 1:length(XYZ)
                if sum(XYZ(h,:)) == 0
                    z(h) = 1;
                end
            end
            z=logical(z);
            XYZ(z,:)=[]; 
            eval(['S.MedRes.XYZ_TADs_mapp.chr_' num2str(i) '.block_' num2str(t+1) '= XYZ;' ]);    
            
        end
        

    end
    
    % stores in S the number of TADs identified in each chromosome
    eval(['S.HiRes.num_TADs_chrs.chr_' num2str(i) '= t+1;']); 
end

%% composing MedRes-chromosomes with HiRes-TADs  ->  building HiRes chromosomes 3D structures

for i = 1:C     % for each chromosome
    t = eval(['S.HiRes.num_TADs_chrs.chr_' num2str(i) ';']);
    for j = 1:t     % each TAD of each chromosome
        if eval(['isempty(S.HiRes.XYZ_mapp_regs.chr_' num2str(i) '.block_' num2str(j) ') == 0'])
            
            eval(['XYZ_Med = S.MedRes.XYZ_TADs_mapp.chr_' num2str(i) '.block_' num2str(j) ';']);    % splitting of the MedRes chromosome structure in the MedRes structure of the considered TAD
            eval(['XYZ_Hi = S.HiRes.XYZ_mapp_regs.chr_' num2str(i) '.block_' num2str(j) ';']);      % HiRes structure of the considered TAD
            n_Med = size(XYZ_Med, 1);   % number of points of the MedRes TAD structure
            n_Hi= size(XYZ_Hi, 1);  % number of points of the HiRes TAD structure
        
            % subsampling of the HiRes TAD structure 
            XYZ_sub = zeros(n_Med, 3);
            XYZ_sub(1,:) = XYZ_Hi(1,:);
            XYZ_sub(n_Med,:) = XYZ_Hi(n_Hi,:);
            delta = fix((n_Hi)/(n_Med));
            k = 1 + delta;      % index of the next considered point on the HiRes structure
            h = 2;      % index of the next considered point on the subsampled structure
            XYZ_sub(h,:) = XYZ_Hi(k,:);
            while k<n_Hi & h<n_Med & XYZ_sub(h+1,:) == [0 0 0]
                h = h+1;
                k = k+delta;
                XYZ_sub(h,:) = XYZ_Hi(k,:); 
            end
        
            % finding Proctustes transform between the MedRes structure corresponding to the considered TAD and the subsampled HiRes TAD structure
            [~, ~, transform] = procrustes(XYZ_Med, XYZ_sub);
            R = transform.T;
            b = transform.b;
            c = repmat(transform.c(1,:)', 1, n_Hi)';
            XYZ_Hi_aligned = b * XYZ_Hi * R + c;
            eval(['S.MedRes.XYZ_Hi_aligned.chr_' num2str(i) '.block_' num2str(j) ' = XYZ_Hi_aligned;' ]);   % stores in S the high resolution TAD structures, aligned to its chromosome
        
        else
            eval(['S.MedRes.XYZ_Hi_aligned.chr_' num2str(i) '.block_' num2str(j) ' = zeros(S.HiRes.TADs_dim.chr_' num2str(i) '.block_' num2str(j) ', 3);' ]);   % if the considered matrix block is completely zero, we assign zero coordinates to the correspondign structure (matlab will discard them)
        end
    end
end

% composing high resolution chromosomes structures and storing them in S.MedRes.XYZ_Hi_chrs.chr_#c
for i = 1:C
t = eval(['S.HiRes.num_TADs_chrs.chr_' num2str(i) ';']);
X = [];
    for j = 1:t
        X = cat(1, X, eval(['S.MedRes.XYZ_Hi_aligned.chr_' num2str(i) '.block_' num2str(j) ]));
    end
    eval(['S.MedRes.XYZ_Hi_chrs.chr_' num2str(i) ' = X; ' ]); 
end

%% procrustes tra S.MedRes.XYZ_mapp_regs e S.LowRes.XYZ_chrs
for i = 1:C
    eval(['XYZ_Med = S.MedRes.XYZ_mapp_regs.chr_' num2str(i) ';']);
    eval(['XYZ_Low = S.LowRes.XYZ_chrs.XYZ_' num2str(i) ';']);
    n_Med = size(XYZ_Med, 1);
    n_Low = size(XYZ_Low, 1);
        
    % subsampling of the MedRes chromosome structure
    XYZ_sub = zeros(n_Low, 3);
    XYZ_sub(1,:) = XYZ_Med(1,:);
    XYZ_sub(n_Low,:) = XYZ_Med(n_Med,:);
    delta = fix((n_Med)/(n_Low));
    k = 1 + delta;  % index of the next considered point on the MedRes structure
    h = 2;   % index of the next considered point on the subsampled structure
    XYZ_sub(h,:) = XYZ_Med(k,:);
    while k<n_Med & h<n_Low & XYZ_sub(h+1,:) == [0 0 0]
        h = h+1;
        k = k+delta;
        XYZ_sub(h,:) = XYZ_Med(k,:); 
    end
    
    % finding Proctustes transform between the LowRes structure corresponding to the considered chromosome and the subsampled MedRes chromosome structure
    [~, ~, transform] = procrustes(XYZ_Low, XYZ_sub);
    R = transform.T;
    b = transform.b;
    c = repmat(transform.c(1,:)', 1, n_Med)';
    XYZ_Med_aligned = b * XYZ_Med * R + c;
    eval(['S.LowRes.XYZ_Med_aligned.chr_' num2str(i) ' = XYZ_Med_aligned;' ]);  % stores in S the medium resolution chromosome structures, aligned to the corresponding low resolution genome portion
    
    % alignment of high resolution chromosomes structures to the whole genome
    eval(['W = S.MedRes.XYZ_Hi_chrs.chr_' num2str(i) ';' ]);   % high resolution chromosomes structures
    n_Hi = size(W, 1);
    w = repmat(transform.c(1,:)', 1, n_Hi)';
    XYZ_Hi_aligned = b * W * R + w; 
    eval(['S.LowRes.XYZ_Hi_aligned.chr_' num2str(i) ' = XYZ_Hi_aligned;' ]);
end

% composing HIGH RESOLUTION WHOLE-GENOME STRUCTURE (stored in S.LowRes.XYZ_GENOME)
X = [];
for i = 1:C
    X = cat(1, X, eval(['S.LowRes.XYZ_Hi_aligned.chr_' num2str(i) ]));
end
eval('S.LowRes.XYZ_GENOME = X; ');







