function val = ExponentialShape(a, xi, der, p, h)
if a == 1
    if der == 0
        val = ( exp(p*h*(xi - 1) / 2) - exp(-p*h*(xi - 1) / 2) ) / (exp(- p * h) - exp( p * h));
    elseif der == 1
        val = ( exp(p*h*(xi - 1) / 2) + exp(-p*h*(xi - 1) / 2) ) / (exp(- p * h) - exp( p * h))* p * (h / 2);
    end
elseif a == 2
    if der == 0
        val = (-exp(p*h*(xi + 1) / 2) + exp(-p*h*(xi + 1) / 2) ) / (exp(- p * h) - exp( p * h));
    elseif der == 1
        val = (-exp(p*h*(xi + 1) / 2) - exp(-p*h*(xi + 1) / 2) ) / (exp(- p * h) - exp( p * h)) * p * (h / 2);
    end 
end