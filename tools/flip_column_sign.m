function U = flip_column_sign(U, sta)
% flip signs of columns to match that of sta
% U - columns vectors
% sta - a row vector

sign_of_inner_products = (sta*U)>=0;
% find negatives
idx_neg = find(sign_of_inner_products<0);

U(:,idx_neg) = -U(:,idx_neg);


end
