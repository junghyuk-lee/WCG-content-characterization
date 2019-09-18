%% WCG image content characterization
% To run this code, you need to download cssim code in 
% https://telin.ugent.be/~bortiz/color

D = content_characterization('test.tif');

function D = content_characterization(imname)
% Pre-defined gamut conversion matrix of Rec.2020, P3, Rec.709, and Toy
load conversion_matrix.mat

N = 3;
D = zeros(N,1);

M_0 = M_2020;
M = {M_P3,M_709,M_Toy};

% I_0: WCG source (Rec.2020) image
I_0 = im2double(imread(imname));
for n=1:N
    % Generate gamut-reduced image I_n with all colors in G_n
    I_n = gamut_clipping(I_0, M_0, M{n});
    
    % Calculate d_n = PD(I_n,I_0)
    d_n = PD(I_n,I_0);
    
    % Obtain vector D = [d_1,...,d_N]'
    D(n) = d_n;
end
end
%% WCG database characterization
function [C, C_total, U, U_total] = database_characterization(Ds)
N = size(Ds,2); % number of reduced gamuts

% Coverage
maxD = 2;
C = (max(Ds,[],1)-min(Ds,[],1))/maxD;

% Total coverage
D_norm = Ds/maxD;
k = convhull(D_norm(:,1),D_norm(:,2));
C_total = sqrt(polyarea(D_norm(k,1),D_norm(k,2)));

% Uniformity
nbins = 10;
U = zeros(1,N);
for i=1:N
    [H, ~]= histcounts(Ds,nbins);
    H = H/size(Ds,1);
    U(i) = -nansum(H.*log10(H));
end

% Total uniformity
Xedges = 0:0.2:2;
Yedges = 0:0.2:2;
H= histcounts2(C(:,1),C(:,2),Xedges,Yedges,'Normalization','probability');
H = H(:);
U_total = -nansum(H.*log10(H)/2);

end
%% WCG source content selection
function selected_contents = content_selection(image_set, Ds)
k = 5; % number of clusters
n = 1; % number of selected contents for each cluster
idx = kmeans(Ds,k);

selected_idx = zeros(n*k,1);
for i=1:k
    cluster_i = find(idx==i);
    s = randperm(length(cluster_i),k);
    selected_idx((i-1)*n+1:i*n) = cluster_i(s);
end
selected_contents = image_set(selected_idx);
end

%% etc.
function out = gamut_clipping(img, M1, M2)
    % Gamut clipping function
    % M1: Gamut conversion matrix of input img's gamut
    % M2: Gamut conversion matrix of target gamut
    
    [h,w,~] = size(img);
    out = reshape(transpose(M2\M1*transpose(reshape(img,[],3))),h,w,3);
    out = max(min(out,1),0);
end

function d_n = PD(I1,I2)
    % PD prediction function using cssim metric

    cssim = cd04_color_SSIM(I1,I2);
    d_n = 2/(1+10^(-3.5*(1.9-cssim)));
end
