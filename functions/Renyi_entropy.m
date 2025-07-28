function [h] = Renyi_entropy(TF,alpha)
[row,col]=size(TF);
TF_abs=abs(TF);
%TF_abs_alpha=abs(TF).^(alpha);
%{
dec=0;
for i=1:row
    for j=1:col
        %summation=summation+TF_abs_alpha(i,j);
        dec=dec+TF_abs(i,j);       
    end
end
%}
dec = sum(sum(sum(TF_abs)));
TF_norm=TF./dec;
TF_norm_alpha=abs(TF_norm).^alpha;
%{
summation=0;
for i=1:row
    for j=1:col
        summation=summation+TF_norm_alpha(i,j);      
    end
end
%}
summation = sum(sum(sum(TF_norm_alpha)));
summation=log2(summation);
summation=(summation)/(1-alpha);
h=summation;
end

