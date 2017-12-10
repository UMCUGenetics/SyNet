function [Net_Adj, Net_Threshold] = SelectTopFromNet(Net_Adj, SelectionType, MAX_N)

%% Initialization
% Net_Adj = randi(10, 100);
% SelectionType = 'p';
% MAX_N = 70;
n_gene = size(Net_Adj,1);
Net_Adj(1:n_gene+1:end) = 0; % Set diagonal to zero

%% Sorting links
if min(nonzeros(Net_Adj(:)))<1e-9, error('Not implemented for smaller thab 1e-9 values.'); end
Net_Eps = triu(Net_Adj + rand(n_gene)*1e-10, 1); % Add small variation to make sure links with same weight do not exists
[scr_val, scr_ind] = sort(Net_Eps(:), 'Descend');
n_lnk = numel(scr_ind);

%% Selecting top items
if strcmpi(SelectionType, 'P')
	clear scr_ind
	fprintf('Selecting top %d interactions.\n', MAX_N);
	Net_Threshold = scr_val(MAX_N);
	Net_Adj(Net_Eps < Net_Threshold) = 0;
	clear Net_Eps scr_val
else
	clear Net_Eps scr_val
	MAX_N_GENE = min([n_gene MAX_N]);
    fprintf('Selecting top %d genes.\n', MAX_N_GENE);
    top_pair = zeros(n_lnk, 2, 'uint16');
    [top_pair(:,1), top_pair(:,2)] = ind2sub([n_gene n_gene], scr_ind);
    if any(top_pair(:)>=65535), error('Not implemented for such a large number of genes!'); end
    [~, top_ind] = unique(top_pair', 'Stable'); % (1:1e6,:)
    Pivot_Index = ceil(top_ind(MAX_N_GENE)/2);
    Net_Threshold = Net_Adj(scr_ind(Pivot_Index));
    Net_Adj(scr_ind(Pivot_Index+1:end)) = 0;
    clear top_ind scr_ind
end

%% Adjustments
Net_Adj = max(Net_Adj, Net_Adj');
fprintf('[%d] genes and [%d] links are left in the network.\n', sum(sum(Net_Adj)~=0), numel(nonzeros(triu(Net_Adj))));
end