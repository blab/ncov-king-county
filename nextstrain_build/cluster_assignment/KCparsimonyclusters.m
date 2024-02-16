% get the clade mapping 
exlude = {'USA/WA1/2020'};
sample_cutoff = '2022-06-03';
%% distance based clustering
date_cutoff = datenum(sample_cutoff);

tree = phytreeread('../results/other_tree.nwk');
% get all WA sequences
f = fopen('../results/other_sub_subsampled_metadata.tsv');
line = strsplit(fgets(f), '\t');
div_id = find(ismember(line,'division'));
date_id = find(ismember(line,'date'));
originating_id = find(ismember(line,'originating_lab'));
submitting_id = find(ismember(line,'submitting_lab'));
location_id = find(ismember(line,'location'));
strain_id = find(ismember(line,'strain'));
c=1;
id = cell(0,0);
clade_membership = cell(0,0);
division = cell(0,0);
lab = cell(0,0);
sub = cell(0,0);
date = cell(0,0);
date_val = zeros(0,0);
strain = cell(0,0);
while ~feof(f)
    line = strsplit(fgets(f), '\t','CollapseDelimiters', false);
    id{c,1} = line{1};
    strain{c,1} = line{strain_id};
    date{c,1} = line{date_id};
    lab{c,1} = line{originating_id};
    sub{c,1} = line{submitting_id};


    if ~isempty(find(ismember(id{c,1}, exlude)))
        division{c,1} = 'NA';
    elseif ~contains(line{date_id}, 'X') && sum(date{c,1}=='-')==2 && datenum(date{c,1})<=date_cutoff
        date_val(c,1) = datenum(date{c,1});
        division{c,1} = line{location_id};
    else
        division{c,1} = 'NA';
    end
    c=c+1;
end
fclose(f);

% get all leafnames
leafs = get(tree, 'leafnames');

%%
is_kc = find(ismember(division, 'King County WA'));

is_kc_leaf = false(length(leafs), 1);
for i = 1 : length(leafs)
    if sum(ismember(id(is_kc), leafs{i}))==1
        is_kc_leaf(i) = true;
    end
end

nodenames = get(tree, 'nodenames');    
    
% get the connectivity matrix and the distances between nodes
mat = getmatrix(tree);
%dist = get(tree, 'Distances');
    
%% initialize the location vector
location = cell(length(nodenames),1);
visited = false(length(nodenames),1);

for j = 1 : length(is_kc_leaf)           
    visited(j) = true;
    if is_kc_leaf(j)
        location{j} = "King County";         
    else
        location{j} = "Outside";
    end
end

% upwards path of the parsimony calling    
not_visited = find(~visited);
while ~isempty(not_visited)
    disp(length(not_visited))
    not_visited = find(~visited);
    for j = 1 : length(not_visited)
        
        children = find(mat(not_visited(j),:));
     %   distances = dist(not_visited(j), children);
        if ~isempty(location{children(1)}) && ~isempty(location{children(2)})
  %          if sum(distances>0)==2
                % get the distance of the new sample from a basel sequence
                int = intersect(location{children(1)}, location{children(2)});
                if isempty(int)
                    location{not_visited(j)} = [location{children(1)}, location{children(2)}];
                else
                    location{not_visited(j)} = int;
                end
%            else
 %              location{not_visited(j)} = [location{children(1)}, location{children(2)}];
      %      end
            visited(not_visited(j)) = true; 
        end
    end
end





f = fopen('node_locations.tsv', 'w');
for i = 1 : length(nodenames)
    if contains(nodenames{i},'NODE')
        if strcmp(location{i},'King County')
            fprintf(f, '%s\tKing County\n',nodenames{i});
        else
            fprintf(f, '%s\tNot King County\n',nodenames{i});
        end
    end
end
fclose(f);

%% downwards calling 
visited(length(is_kc_leaf)+1:end) = false;
visited(end) = true;

not_visited = find(~visited);
while ~isempty(not_visited)
    not_visited = find(~visited);
    for j = length(not_visited) : -1 : 1
        parent = find(mat(:,not_visited(j)));    
    %    parent_dist = dist(parent,not_visited(j));
        
        uni_loc = unique(location{not_visited(j)});

        if sum(ismember(not_visited, parent))==0
       %     if parent_dist == 0
        %        location{not_visited(j)} = unique(location{parent});
            if length(uni_loc) < length(location{not_visited(j)})
                combined_locs = [location{not_visited(j)} location{parent}];
                uni_comb = unique(combined_locs);        
                freqs = zeros(length(uni_comb),1);
                for k = 1 : length(uni_comb)
                    freqs(k) = sum(combined_locs==uni_comb{k});
                end
                ind_max = find(freqs==max(freqs));
                location{not_visited(j)} = uni_comb(ind_max);
            else
                int = intersect(location{not_visited(j)}, location{parent});
                if ~isempty(int)
                    location{not_visited(j)} = int;
                end
            end
            if length(location{not_visited(j)})>1
                location{not_visited(j)};
            end
            visited(not_visited(j)) = true;
        end
    end
end


%%

% get all nodes that are only in KC
onlyKC = false(length(visited),1);   
uncertainCOunt = 0;
for j = 1 : length(onlyKC)
    if length(location{j})==1 && location{j}=="King County"
        onlyKC(j) = true;
    elseif length(location{j})==2
       % disp(location{j})
        onlyKC(j) = true;
        uncertainCOunt = uncertainCOunt+1;
    end
end




isInKC = find(onlyKC(1:length(is_kc_leaf)));

kc_names = str2double(nodenames(isInKC));    

clustered_node = false(size(nodenames));
clustered_node(isInKC) = true;


% for each isInWA leaf, get all parent nodes in WA
KCParent = cell(length(isInKC),1);
for a = 1 : length(isInKC)
    parent = find(mat(:,isInKC(a)));
    while onlyKC(parent)
        clustered_node(parent) = true;
        KCParent{a} = [KCParent{a} parent];
        parent = find(mat(:,parent));
    end
end

in_cluster_mat = zeros(length(isInKC), length(isInKC));

new_mat = zeros(size(in_cluster_mat));
% get for each pair of leaves if they are in the same cluster
for a = 1 : length(isInKC)
    %disp(a)
    if ~isempty(KCParent{a})
        for b = a+1 : length(isInKC)
            if ~isempty(KCParent{b})
                int = intersect(KCParent{a}, KCParent{b});
                if ~isempty(int)
                    new_mat(a,b) = 1;
                end
            end
        end     
    end
end

%%
kc_leafs = leafs(is_kc_leaf);

connected = new_mat;

stop=false;
while ~stop
    [a,b] = find(connected==1,1,'first');
    if isempty(a)
        stop = true;
    else
        % combine a and b
        kc_leafs{a} = [kc_leafs{a} ',' kc_leafs{b}];
        % combine the rows of a and b
        for i = 1 : length(kc_leafs)
            connected(a,i) = max([connected(a,i),connected(b,i)]);
        end
        % remove b
        connected(b,:) = 0;
        connected(:,b) = 0;
        kc_leafs{b} = 'NA';
    end
end

% get the cluster
cl_ind = find(~ismember(kc_leafs, 'NA'));
kc_clusters = kc_leafs(cl_ind);


f = fopen('other_new_cluster_assignment.tsv', 'w');
fprintf(f,'strain\tcluster\n');
for a = 1 : length(kc_clusters)
    seqs = strsplit(kc_clusters{a}, ',');
    for b = 1 : length(seqs)
        if length(seqs{b})>2
            fprintf(f, '%s\t%d\n', seqs{b}, a);
        end
    end
end
fclose(f);

