function val = HermiteShape(a, xi, der, h)
if a == 1
    if der == 0
        val = 0.25 * (xi - 1) * (xi - 1) * (xi + 2);
    elseif der == 1
        val = 0.75 * (xi * xi - 1);
    elseif der == 2
        val = 1.5 * xi;
    end
elseif a == 2
    if der == 0
        val = 0.125 * h * (xi + 1) * (xi - 1) * (xi - 1);
    elseif der == 1
        val = 0.125 * h * (xi - 1) * (3*xi + 1);
    elseif der == 2
        val = 0.25 * h * (3*xi - 1);
    end
elseif a == 3
    if der == 0
        val = 0.25 * (xi + 1) * (xi + 1) * (2 - xi);
    elseif der == 1
        val = 0.75 * (1 - xi * xi );
    elseif der == 2
        val = -1.5 * xi;
    end
elseif a == 4
    if der == 0
        val = 0.125 * h * (xi + 1) * (xi + 1) * (xi - 1);
    elseif der == 1
        val = 0.125 * h * (xi + 1) * (3*xi - 1);
    elseif der == 2
        val = 0.25 * h * (3*xi + 1);
    end
end






