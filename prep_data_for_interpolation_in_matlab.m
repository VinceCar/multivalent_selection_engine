load(Source_kD_file)
X1 = unique(X(:,1));
X2 = unique(X(:,2));
X3 = unique(X(:,3));
X4 = unique(X(:,4));
X5 = unique(X(:,5));
[X_g_1,X_g_2,X_g_3,X_g_4,X_g_5] = ndgrid(X1,X2,X3,X4,X5);
y_1_g = zeros(size(X_g_1));
y_2_g = zeros(size(X_g_1));
y_3_g = zeros(size(X_g_1));

for ii = 1:length(y_1_g(:))
    wel = find((abs(X_g_1(ii)-X(:,1))<eps).*(abs(X_g_2(ii)-X(:,2))<eps).*(abs(X_g_3(ii)-X(:,3))<eps).*(abs(X_g_4(ii)-X(:,4))<eps).*(abs(X_g_5(ii)-X(:,5))<eps));
    if ~isempty(wel)
        y_1_g(ii)=y(wel,1);
        y_2_g(ii)=y(wel,2);
        y_3_g(ii)=y(wel,3);
    else
        y_1_g(ii)=NaN;
        y_2_g(ii)=NaN;
        y_3_g(ii)=NaN;
        
    end
end