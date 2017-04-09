function [J grad] = myTestFunction(p)
#
# f = @(p) myTestFunction(p)
# [X, fX, i] = fmincg(f, p)
J = p'*p + 5;
grad = 2*p;
end;