figure;
sigs =[];
for i = [1,2,4:23]
    for j =1:2
        ind = HOT_sig(:,1) == i & HOT_sig(:,2) == j; 
        x= HOT_sig(ind,4);
        y=HOT_sig(ind,5);
        z=HOT_sig(ind,3);
        %plot(x,y)
        %hold on
        sigs = [sigs; i, j, z(find(y==min(y))),x(find(y==min(y))),min(y)];
    end
end

    