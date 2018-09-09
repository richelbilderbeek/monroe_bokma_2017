function [PhT, E] = Evo(Tree,sa2,sc2)
%returns phenotypes PhT of species evolved on phylogenetic Tree with rates
%of anagenetic and cladogenetic evolution given by variance-covariance 
%matrices sa2 and sc2, respectively.
%E is an array with indexes of the extant species.

nt = length(sa2); %number of traits to be simulated
B = get(Tree,'Pointers'); %list of child nodes
D = get(Tree,'Distances'); %list of branch lenghts
pi = max(max(B))+1;%deepest node
PhT(pi,:) = zeros(nt,1); %allocate memory
E(pi) = 0;%set deepest node to extinct

for j = size(B,1):-1:1 %start at root
    dc = round(rand(1,length(B))); %randomly assign which species is novel
    dc(2,:) = 1-dc; %corresponding sisters
    for si = 1:2 %for each sister
        d = randn(1,nt).*sqrt(sc2)*dc(si,j); %cladogenetic change
        d = d + randn(1,nt).*sqrt(D(B(j,si))*sa2); %anagenetic change
        PhT(B(j,si),:) = PhT(pi,:) + d; %update phenotype of this node
        E(B(j,si)) = E(pi) + D(B(j,si)); %update time of this node
    end %si
    pi = pi-1;
end %j
E = find(abs(E-max(E))<0.00001); %find indices that are very close to max ARBITRARY
