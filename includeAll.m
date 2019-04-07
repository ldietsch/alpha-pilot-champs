function Xref = includeAll(xref)
k = size(xref,2);
Xref = zeros(3, 4*k + 1);
j = 1;
for i = 1:k
    Xref(:, j) = xref(:,i);
    j = j + 4;
end
Xref = Xref(:);

end

