% Computes Sampson distance (a fast way)
function [SampDist, AlgebDist_R] = FMDistances(FM,X)

X1 = X(1:2,:)';  
X2 = X(4:5,:)';  
    rFhat2 = [X2 ones(length(X2),1)]* FM;   % [x' y' 1] * F
    AlgebDist_R = rFhat2(:,1) .* X1(:,1) + rFhat2(:,2) .* X1(:,2) + rFhat2(:,3);
    drx =  FM(1,1)*X2(:,1) + FM(2,1)*X2(:,2)+FM(3,1); % F11 . x' + F21 . y' + F31
    dry =  FM(1,2)*X2(:,1) + FM(2,2)*X2(:,2)+FM(3,2);
    drxp = FM(1,1)*X1(:,1) + FM(1,2)*X1(:,2)+FM(1,3); % F11 . x + F12 . y + F13
    dryp = FM(2,1)*X1(:,1) + FM(2,2)*X1(:,2)+FM(2,3);
     
    SampDist = AlgebDist_R ./ sqrt(drx .^ 2 + dry .^ 2 + drxp .^ 2 + dryp .^ 2);
end
