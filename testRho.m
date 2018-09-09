function [ns,robs,p]=testRho(T,MassNames,Masses)

[LeafNames,B,nsi,nb]=get(T,'LeafNames','Pointers','NumLeaves','NumBranches');%retrieve number of leaves from tree
[nm,T]=OrderMasses(T,LeafNames,MassNames,Masses,nsi);%put masses in order of species on the tree
if (length(nm)<3) | (length(nm)<0.5*nsi)
    ns=length(nm); robs=NaN; p=NaN;
    return
end
[LeafNames,B,ns,nb]=get(T,'LeafNames','Pointers','NumLeaves','NumBranches');%retrieve number of leaves from tree
[nm,T]=OrderMasses(T,LeafNames,MassNames,Masses,ns);%put masses in order of species on the tree

robs=GetRhoMax(B,nm,nb,ns);%calculate rmax with observed masses
%get p-value
sa2=0.1*ones(1,1000); sc2=zeros(1,1000);
[PhT, E] = Evo(T,sa2,sc2); PhT=PhT(E,:);
for j=1:1000
    r(j)=GetRhoMax(B,PhT(:,j),nb,ns);
end
p=(sum(r>=robs)+1)/(1000+1);

end %APC

    function rb=GetRhoMax(B,nm,nb,ns)
        rb=GetRho(B,nm,nb,ns);
        for j=1:ns-1
            B(j,:)=fliplr(B(j,:));
            ra=GetRho(B,nm,nb,ns);
            if ra<rb
                B(j,:)=fliplr(B(j,:));
            else
                rb=ra;
            end
        end
    end %GetRhoMax
 
function r=GetRho(B,nm,nb,ns)%function with imput variables B=pointers, m=masses, nb=number of branches, ns=number of species
    order=B(end,:);%from the end of the tree to the begining
    for j=ns+nb-1:-1:ns+1;%using the last branch to calculate the position of the species
         k=find(order==j);%find the species order from the branch
         order=[order(1:k-1) B(j-ns,:) order(k+1:end)];%put the species in order 
    end
    r=corr(nm(order),[1:ns]','type','Kendall');
end%GetRho

function [nm,T]=OrderMasses(T,LeafNames,MassNames,m,ns)
    nm=zeros(ns,1); nm(:)=NaN;
    for j=1:ns%look up species in array of masses
        k = find(strcmp(LeafNames{j}, MassNames)); %locates and gives the position in the order of the tree to the masses
        if ~isempty(k); %only does the step above if there is a match
            nm(j)=m(k(1));%reorders the masses so that the masses are in the order of the tree
        end
    end %j
    nm(nm==-999)=NaN;
    toprune = find(isnan(nm));%get species lacking mass data
    T = prune(T, toprune);%remove these species from the tree 
    nm=nm(~isnan(nm));
end%OrderMasses


