function Us = flip_column_sign(U, sta)
% flip signs of columns to match that of sta
% U - columns vectors
% sta - a row vector

sign_of_inner_products = sign(U*sta);


Us = bsxfun(@times, U, sign_of_inner_products);


end
